from classes import Graph
import sqlite3 as sql
import numpy as np


x=range(8900, 10400+1)
print(len(x))
# dgrm = np.loadtxt('dgrm.csv', delimiter = ';')
# dgr = np.loadtxt('dgr.csv', delimiter = ';')
# # tttt = np.loadtxt('tttt.csv', delimiter = ';')
# tttt=range(8900, 10400+1)
# frmax=9400
#
# Y=[]
# X=[]
# l=[]
# cn=[]
# m=[]
# lw=[]
# g=Graph(path='finals', name='set_std', project = 'x', show=True, x_label='Set #', y_label = 'Set standart deviation /nm.s^(-2)', title='Set gravity at top of the drop')
# # g.plotXY(x=[tttt], y=[dgrm], mark=['k-'], columns_name=['Mean'], legend =['Mean'], lw=1)
# for i in range(len(dgr)):
#     # g.plotXY(x=[tttt], y=[dgr[i,:]], mark=['C'+str((i)%10)+ '-'], columns_name=['Set ' + str(i+1)], legend =['Set ' + str(i+1)])
#     X.append(tttt)
#     Y.append(dgr[i,:])
#     l.append('Set ' + str(i+1))
#     cn.append('Set ' + str(i+1))
#     m.append('C'+str((i)%10)+ '-')
#     lw.append(0.3)
#
# X.append(tttt)
# Y.append(dgrm)
# l.append('Mean')
# cn.append('Mean')
# m.append('k-')
# lw.append(1)
#
#
# g.plotXY(x=[tttt], y=[[0 for i in range(len(tttt))]], mark=['b-'], columns_name='zero', legend ='', lw=[0.3])
# g.plotXY(x=[[frmax, frmax]], y=[[-10, 10]], mark=['b-'], columns_name='zero', legend ='', lw=[0.3])
# g.plotXY(x=X, y=Y, mark=m, columns_name=cn, legend =l, lw=lw)
#
#
#
#
# # g.error_bar(range(1, self.nset + 1), self.stodch, self.stodchs, 'r')
# # g.saveSourceData()
# g.save()
