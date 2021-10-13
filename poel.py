import os, sys
import csv
from datetime import datetime, timedelta
import glob
from classes import dropFile
from math import pi, sin, cos
import sqlite3 as sql
import numpy as np

ps=10
sens_tn=1
sens_tx=150*ps
frmax=940*ps
frmin=15*ps
frmaxplot=1040*ps
nforfft=4501


tt=np.loadtxt('tt.txt')
resmm=np.loadtxt('resmm.txt', delimiter = ';')
resm=np.loadtxt('ress.csv', delimiter=';')


matr_db=sql.connect(os.getcwd()+'/zaloha/data.db')
cursor=matr_db.cursor()
# print(type(cursor))
cursor.execute('select Set1, Res from results where Accepted = 1')
r=cursor.fetchall()


#
cursor.execute('''select count(*) from results
where Accepted = 1
group by Set1''')
count=cursor.fetchall()

countAcc=0
for c in count:
    countAcc+=c[0]


#
# cursor.execute('select Accepted, Set1, Drop1 from results')
# d=cursor.fetchall()

tin=np.linspace(tt[frmin-1], tt[frmax-1], nforfft)
tinc=np.linspace(tt[0], tt[frmaxplot-1], nforfft)
ttx2=tt[0:frmaxplot]

# resxxx=np.interp(tin, ttx2, resmm)
ttx=tt[frmin-1:frmax]

x=int((nforfft-1)/2)

yfd=np.zeros((len(r),nforfft), dtype=complex)
yfdMeanBySet=np.zeros((15,x))
yfdMean=np.zeros((1,x))

it=0
for set, res in r:
    ress=res.split(',')
    ress=np.array([float(x) for x in ress])
    ress=ress[frmin-1:frmax]

    resd=np.interp(tin, ttx, ress)
    resd=resd-np.mean(resd)

    fft= 2/nforfft*np.fft.fft(resd)

    yfd[it, :] =fft

    l=np.absolute(fft[0:x])/count[set-1][0]
    yfdMeanBySet[set-1,:]+=np.real(l)

    l=np.absolute(fft[0:x]/countAcc)
    yfdMean[0,:]+=np.real(l)


    it+=1

yfs=np.zeros((15, 4501), dtype = complex)
yfsa=np.zeros((15,x))
for i in range(0,15):
    ress=np.interp(tin, ttx, resm[i,frmin-1:frmax])

    ressm=ress-np.mean(ress)

    fft = 2/nforfft*np.fft.fft(ressm)

    yfs[i, :] = fft

    yfsa[i, :] = np.real(np.absolute(fft[0:x]))







np.savetxt('pokus.csv',yfdMean , delimiter=';')

















# d=int(len(r[0][1].split(',')))
# res=np.zeros((15,d))
#
#
#
# res_mean=np.zeros((1,d))
#
# # print(r[1][0].split(','))
# for i in range(0,len(r)):
#     set=r[i][0]
#     ress=r[i][1:][0].split(',')
#     ress=np.array([float(x) for x in ress])
#
#     # for j in range(len(ress)):
#     #     c=float(ress[j])
#     #
#     #     if c > 1:
#     #         ress[j]=0
#     #     else:
#     #         ress[j]=c
#
#
#     res[set-1,:]=res[set-1,:]+ress
#     # print(ress[3])
#     res_mean=res_mean+ress
#
#
# # for i in range(15):
# #     print(res[i,:])
#
# dgl=np.zeros((15,1500))
# dgr=np.zeros((15,1500))
# dglm=np.zeros((1,1500))
# dgrm=np.zeros((1,1500))
# c=0
# for i in range(len(res)):
#     res[i,:]=res[i,:]/count[i][0]
#     c+=count[i][0]
# np.savetxt('resmm.csv', res_mean/c,delimiter=';')
#
#     for j in range(sens_tn, sens_tx+1):
#
#         x=tt[j-1:frmax]
#         y=res[i, j-1:frmax]
#
#         koef=np.polyfit(x,y,deg=2)
#
#         dgl[i,j-1]=koef[0]*2
#
#         x=tt[frmin-1:j+8899]
#         y=res[i, frmin-1:j+8899]
#
#         koef=np.polyfit(x, y, deg=2)
#
#         dgr[i, j-1]=koef[0]*2
#
#     dglm = dglm + dgl[i, :]
#     dgrm = dgrm + dgr[i, :]
#
# dglm=dglm/15
# dgrm=dgrm/15
#
# sensa_tn=3*ps
# sensa_tx=80*ps
# sensa_bn=900*ps
# sens_bn=frmaxplot-150*ps
# sensa_bx=1000*ps
#
#
# rozd=sensa_tn-sens_tn
# celk=sensa_tx-sensa_tn
#
# dglc=dgl[:, rozd:rozd+1+celk]
#
# dglrms=np.sqrt(np.sum(np.square(dglc.transpose()), axis=0))/np.sqrt(celk+1)
# print(dglrms)
#
# rozd=sensa_bn-sens_bn
# celk=sensa_bx-sensa_bn
# dgrc=dgr[:, rozd:rozd+1+celk]
# dgrrms=np.sqrt(np.sum(np.square(dgrc.transpose()), axis=0))/np.sqrt(celk+1)
# print(dgrrms)

# np.savetxt('dglrms.csv', dglrms, delimiter = ';')
# np.savetxt('dgrm.csv', dgrm, delimiter = ';')
#
#
#
#
# #
# #
# np.savetxt("dgr.csv", dgr, delimiter=";")


# for i in range(0,len(res)):
#     print(res[i,2])
#
#
# res_mean=res_mean/c
# for i in range(0, len(res_mean)):
#     print(res_mean[i])




#
# for dropfile in glob.glob(os.getcwd()+'/data\*.drop.txt'):
#     d=dropFile(dropfile)
#     columns_dropfile=d.dropHeader4()
#     lines=d.dropLines()
#
#
# # open and load file from IERS
# file=open( os.getcwd()+'/finals/finals2000A.all.csv', 'r')
# # reader = csv.DictReader(file, delimiter=';')
# reader = csv.reader(file, delimiter=';')
# rows=list(reader)
# file.close()
#
# # date of first day
# a=datetime(int( lines[2].split()[4]),1,1) + timedelta(int( lines[2].split()[3]))
#
# month=(str(a.month))
# if len(month)==1:
#     month='0'+month
#
# day=(str(a.day))
# if len(day)==1:
#     day='0'+day
#
# date=[str(a.year), month,  day]
#
# # get unique DOY
# doys=[]
# years=[]
# for l in  lines:
#     doy=int(l.split()[3])
#     if doy not in doys:
#         doys.append(doy)
#         years.append(int(l.split()[4]))
#
# # get index of date in finals file
# i=0
# for row in reversed(rows):
#     if row[1:4]==date:
#         today=row
#         break
#     i+=1
#
# # coordinates of pole from finals to interpolation
# x=[]
# y=[]
# for j in range(-i,-i+len(doys)+1):
#     x.append(rows[j][5])
#     y.append(rows[j][7])
#
# fi= 49.91381 *pi/180
# lam=   14.78559*pi/180
#
# print(x)
#
# # compute pole corrections
# dg=[]
# for l in  lines:
#     line=l.split()
#
#     doy=line[3]
#
#     xInt=[float(x[doys.index(int(doy))]), float(x[doys.index(int(doy))+1])]
#     yInt=[float(y[doys.index(int(doy))]), float(y[doys.index(int(doy))+1])]
#
#     time=(line[2].split(':'))
#     time=float(time[2])/3600+float(time[1])/60+float(time[0])
#
#     ypole=time*(yInt[1]-yInt[0])/24+yInt[0]
#     xpole=time*(xInt[1]-xInt[0])/24+xInt[0]
#
#     dg.append(-19.139*sin(2*fi)*(xpole*cos(lam)-ypole*sin(lam)))
#     break
#
# # print(xInt)
# print(dg)
