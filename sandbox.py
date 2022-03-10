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
from functions import allan

def allanGraph(a, tau, path):

    p = plt
    p.loglog(tau[:len(a)], [i[0] for i in a], '.r', ms=20)
    p.loglog(tau[:len(a)], [a[0][0]/np.sqrt(tau[i]) for i in range(len(a))], '-r')
    p.plot([tau[:len(a)], tau[:len(a)]], [[a[i][0]-a[i][0]/np.sqrt(a[i][2]) for i in range(len(a))], [a[i][0]+a[i][0]/np.sqrt(a[i][2]) for i in range(len(a))]], '-k', lw = 2)
    p.plot([[i*0.95 for i in tau[:len(a)]], [i*1.05 for i in tau[:len(a)]]], [[a[i][0]+a[i][0]/np.sqrt(a[i][2]) for i in range(len(a))], [a[i][0]+a[i][0]/np.sqrt(a[i][2]) for i in range(len(a))]], '-k', lw = 2)
    p.plot([[i*0.95 for i in tau[:len(a)]], [i*1.05 for i in tau[:len(a)]]], [[a[i][0]-a[i][0]/np.sqrt(a[i][2]) for i in range(len(a))], [a[i][0]-a[i][0]/np.sqrt(a[i][2]) for i in range(len(a))]], '-k', lw = 2)
    p.legend(['Mean values', 'White noise'])
    p.title('Allan deviation - normalized data')
    p.xlabel('Drop number (n)')
    p.ylabel(r'σ (n)/nm.s$^-2$')
    # # p.errorbar(tau, a1s, '*k', tau, sme2, '-')
    p.savefig(path)
    plt.close()

db = sql.connect(r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\res\data.db')
c = db.cursor()
c.execute('select gTopCor from results where Accepted = 1')
res = c.fetchall()
tau=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500]

res = [i[0] for i in res]

a = allan(res, tau)


allanGraph(a, tau)


# rate = 1
# tmstep = 1/rate
#
# f = open(r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\res\Files\X153_glb_allan.csv').read().splitlines()
#
# # tau = np.zeros((len(f)-1, 1))
#
# tau = []
# a1 = []
# a1s=[]
# a2 = []
# a2s=[]
# a3 = []
# a3s=[]
# smn1=[]
# smn2=[]
# smn3=[]
# sme1 = []
# sme2 = []
# sme3 = []
# for i in range(1, len(f)-1):
#     l = f[i].split(';')
#
#     tau.append(float(l[0]))
#     a1.append(float(l[1]))
#     a2.append(float(l[3]))
#     a3.append(float(l[5]))
#
#     a1s.append(float(l[2]))
#     a2s.append(float(l[4]))
#     a3s.append(float(l[6]))
#
#     smn1.append(a1[0]/np.sqrt(tau[i-1]))
#     smn2.append(a2[0]/np.sqrt(tau[i-1]))
#     smn3.append(a3[0]/np.sqrt(tau[i-1]))
#
#     sme1.append(a1s[-1]/np.sqrt(tau[-1]))
#     sme2.append(a2[-1]/np.sqrt(tau[-1]))
#     sme3.append(a3s[-1]/np.sqrt(tau[-1]))
