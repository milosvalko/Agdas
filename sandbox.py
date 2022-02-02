# import pyqtgraph as pg
# import pyqtgraph.exporters
from classes import Graph
import sqlite3 as sql
import numpy as np
from scipy.signal import butter, filtfilt
from functions import rssq, movingAverage
import matplotlib.pyplot as plt
from time import time
from math import atan2


def computeLST(A, z, frmin, frmax):


    x=np.linalg.lstsq(A,z,rcond=None) # solution of LST
    covar = np.matrix(np.dot(A.T, A)).I # covariance matrix
    m02=x[1]/(frmax-frmin-A.shape[1]) # m02=res*res'/(frmax-frmin-k)
    stdX=np.sqrt(np.multiply(np.diag(covar),m02)) # standart devivations of x
    std=np.dot(stdX,stdX)
    res=np.subtract(z,np.matmul(A,x[0])) # res=z-A*X
    m0=np.sqrt(np.dot(res[frmin:frmax], res[frmin:frmax])/(frmax-frmin)*(covar[2,2]))

    return x, covar, m02, std, stdX, res, m0

nset = 15
Lmin = 3
Lmax = 16
frmin=150
frmax=9400
frmaxplot = 10400
hefm = 114.483921664619
gamma = 3.626e-06
tfit = np.loadtxt('tfit.csv', delimiter = ';')
resm = np.loadtxt('resm.csv', delimiter = ';')
v0m = np.loadtxt('v0m.csv', delimiter = ';')
g0m = np.loadtxt('g0m.csv', delimiter = ';')
zzh = np.loadtxt('zzh.csv', delimiter = ';')

ressets2 = np.zeros((frmaxplot, nset))
A0 = np.ndarray((frmax-frmin+1, 3))
AH = np.ndarray((frmax-frmin+1, 3))
AL1 = np.ndarray((frmax-frmin+1, 3))
LA = np.ndarray((frmax-frmin+1, 3))
harmonic_res = np.ndarray((10*Lmax-10*Lmin+1, nset))
T0_m = np.ndarray((10*Lmax-10*Lmin+1, nset))
A0_m = np.ndarray((10*Lmax-10*Lmin+1, nset))
B0_m = np.ndarray((10*Lmax-10*Lmin+1, nset))

xg=[]
xvgg=[]
mx0=[]
mxL1=[]
dg=[]
mdg=[]
vgg=[]
mvgg=[]
# minimum_index_array = []
# minimum_array=[]
LLmin=[] #wave lengths m
T0min=[] #zeros member
A0min=[] #cosinus member
B0min=[] #sinus member
harm_amp=[] #amplitude of harmonic
harm_phase=[] #phase of harmonic
t=time()
for n in range(nset):

    #==========================================================#
    A0[:, 0] = 1
    A0[:, 1] = tfit[frmin-1: frmax]
    A0[:, 2] = 0.5*tfit[frmin-1: frmax]*tfit[frmin-1: frmax]
    b0=self.meanResSets[n, frmin-1:frmax]*1e-9

    x, covar, m02, std, stdX, res, m0=computeLST(A0, b0, frmin-1, frmax)
    mx0.append(m0*1e8)
    xg.append(x[0][2])
    #==========================================================#


    ressets2[:, n] = (self.meanResSets[n, :frmaxplot] - x[0][0] - x[0][1]*tfit[:frmaxplot] - 0.5*xg[n]*tfit[:frmaxplot]*tfit[:frmaxplot])*1e-9

    #==========================================================#

    AL1[:, 0] = 1
    AL1[:, 1] = tfit[frmin-1: frmax]
    AL1[:, 2] = self.v0m_bysets[n]/6*tfit[frmin-1: frmax]**3 + self.g0m_bysets[n]/24*tfit[frmin-1: frmax]**4 - hefm/2*tfit[frmin-1: frmax]*tfit[frmin-1: frmax] - hefm*gamma/24*tfit[frmin-1: frmax]**4
    bL1 = ressets2[frmin-1:frmax, n]

    x1, covar, m02, std, stdX, res, m0=computeLST(AL1, bL1, frmin-1, frmax)
    xvgg.append((gamma + x1[0][2])*1e6)
    mxL1.append(m0*1e6)
    #==========================================================#

    #calculating residuals and their std in range 3-16 cm
    LLL=[]
    for i in range(10*Lmin, 10*Lmax+1):
        LL=i/1000
        LLL.append(LL)

        LA[:, 0] = 1
        LA[:, 1] = np.sin(2*np.pi*zzh[frmin-1: frmax, n]/LL)
        LA[:, 2] = np.cos(2*np.pi*zzh[frmin-1: frmax, n]/LL)

        b = ressets2[frmin-1: frmax, n]

        x2, covar, m02, std, stdX, res, m0=computeLST(LA, b, frmin-1, frmax)

        T0_m[i-10*Lmin, n] = x2[0][0]*1e9
        A0_m[i-10*Lmin, n] = x2[0][1]*1e9
        B0_m[i-10*Lmin, n] = x2[0][2]*1e9
        harmonic_res[i-10*Lmin, n] = np.std(res)*1e9

    #searching of minimal residuals std
    minimum = min(harmonic_res[:, n])
    minimum_index = np.argmin(harmonic_res[:, n])

    T0min.append(T0_m[minimum_index, n])
    A0min.append(A0_m[minimum_index, n])
    B0min.append(B0_m[minimum_index, n])
    LLmin.append(Lmin/100 - 0.001 + minimum_index/1000)
    # print(B0_m[minimum_index, n])

    amplitude = np.sqrt(A0min[-1]*A0min[-1] + B0min[-1]*B0min[-1])
    harm_amp.append(amplitude)

    phase = np.arctan2(B0min[-1], A0min[-1])*(180/np.pi)
    if phase < 0:
        harm_phase.append(phase + 360)
    else:
        harm_phase.append(phase)

    #==========================================================#

    bH = ressets2[frmin-1:frmax, n] - (T0min[n] - A0min[n]*np.sin(2*np.pi*zzh[frmin-1:frmax, n]/LLmin[n]) - B0min[n]*np.cos(2*np.pi*zzh[frmin-1:frmax, n]/LLmin[n]))/1e9
    x1, covar, m02, std, stdX, res, m0=computeLST(A0, bH, frmin-1, frmax)
    dg.append(x1[0][2]*1e8)
    mdg.append(m0*1e8)

    #==========================================================#
    x1, covar, m02, std, stdX, res, m0=computeLST(AL1, res, frmin-1, frmax)
    vgg.append((gamma + x1[0][2])*1e6)
    mvgg.append(m0*1e6)








# np.savetxt('A0.csv', A0min, delimiter = ',')




# d = sql.connect(r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\res1\data.db')
# c = d.cursor()
# c.execute('select avg(EffHeight + CorToEffHeight) from results where Accepted = 1')
# res = c.fetchall()
# print(res[0][0])
# res = [i[0] for i in res]
#
#
#
# g=Graph(path='finals', name='residuals_sifted', project = 'xlb', show=True, x_label='Time /s', y_label = 'Shifted Residuals /nm', title='Set residuals', winsize=(15,10))
# # g.error_bar(x, grad, m0, 'r', ms=5, capsize=5)
# for i in range(harmonic_res.shape[1]):
#
#     g.plotXY(x=[LLL], y=[harmonic_res[:,i]], mark=['-'], columns_name=['pokus'], lw=[1], legend = ['Set {}'.format(i+1)])
# # g.plotXY(x=XX, y=YY, mark=markk, columns_name=col_name, lw=lww)
# # g.text(x=[x[frmin], x[frmax]], y=[0.3, 0.3], t=['Start fringe', 'Final fringe'], c=['b', 'b'])
# # g.text(x=text_x, y=text_y, t=col_name, c=text_color)
# # g.saveSourceData()
# # g.histogram(hist_data = res, fit=True)
# g.save()
