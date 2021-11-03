from time import time as tm
()
# import matplotlib.pyplot as plt
import sqlite3 as sql
# import os, csv
# from CONFIG import statistic, matrDatabase
# from numpy import sqrt, pi, sin, cos
# from datetime import datetime, timedelta
# from classes import dropFile
# import glob
from functions import graph
import numpy as np
from scipy.stats import t
# import py2exe
from glob import glob
from statistics import mean
from CONFIG import matrDatabase
from math import floor
import inspect
import sys

from PyQt5.QtCore import *
from PyQt5.QtGui import *

from PyQt5.QtWidgets import QApplication, QDialog, QMainWindow, QPushButton, QSplashScreen, QProgressBar, QTextBrowser, QLineEdit, QVBoxLayout

import queue

q1=queue.Queue(5)

q1.put(1)
q1.put(2)
q1.put(3)
q1.put(4)
q1.put(5)



for i in range(6,17):
    q1.get()
    q1.put(i)
    print(q1.queue)




















# def roundList(l, index):
#     for j in index:
#         a='{:.'+str(j[1])+ 'f}'
#         l[j[0]]=a.format(l[j[0]])
#
#     return l























# print(np.arange(0,10,2.33))
# a=np.array([[1,2],[3,4]])
# b=np.array([[1,2, 4],[3,4, 8]])
# np.matmul(b,a)
# np.matrix(np.dot(a, a)).I
# np.savetxt('mat.txt', a, delimiter=';')
# # for dropfile in glob.glob(os.getcwd()+'/data/*.drop.txt'):
# #     d=dropFile(dropfile)
# #     columns_dropfile=d.dropHeader4()
# #     lines=d.dropLines()
# #
# # tides=[]
# # time_gr=[]
# # id=0
# # for i in lines:
# #     tide=i.split()[8]
# #     tides.append(float(tide))
# #
# #     time1=i.split()[2].split(':')
# #     time1=int(time1[2])/3600+int(time1[1])/60+int(time1[0])
# #     if id>0:
# #         if time1<time_gr[-1]:
# #             time1+=24
# #     time_gr.append(time1)
# #     id+=1
# #
# # graph(x=[time_gr], y=[tides], xLabel='Time /h', yLabel='Tides /μGal', title='Tidal acceleration',path= os.getcwd()+'/res/Graphs', name='tides.png', hist=False, mark=['b+'])
#
# #
# # xpole=0.0480
# # ypole= 0.3532
# # fi= 49.91381 *pi/180
# # lam=   14.78559*pi/180
# #
# # dg=-19.139*sin(2*fi)*(xpole*cos(lam)-ypole*sin(lam))
# # print(dg)
# # print(0.82)
#
#
#
#
# import sqlite3 as sql
# import numpy as np
# from math import ceil
#
# matr_db=sql.connect(r'c:\Users\Jakub\Desktop\pecny\pyAgdasGui\zaloha\data.db')
# cursor=matr_db.cursor()
# cursor.execute('select Set1, Drop1, Date,  round(g0_Gr,2), "STD", round(CorrToTop,2), round(Tide,2), round(Load, 2), round(Baro, 2), round(Polar, 2), round(gTopCor, 2), round(g0, 2), round(EffHeight,3), round(CorToEffHeight, 3), Accepted  from results')
# res=cursor.fetchall()
#
# nset=max([x[0] for x in res])
#
# Std=[]
# for i in range(nset):
#     cursor.execute('select gTopCor from results where Accepted =1 and Set1 == {}'.format(i+1))
#     v=cursor.fetchall()
#
#     Std.append(np.std(v)*10/np.sqrt(len(v)))
#
#
#
#
# # for i in range(len(res)):
#
#
#
# for i in res:
#     line=list(i[0:4])
#     line.append(ceil(Std[i[0]-1]*100)/100)
#     line.extend(list(i[5:]))
#     print(line)
