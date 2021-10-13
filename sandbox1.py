from allantools import adev, oadev, mdev, totdev
import sqlite3 as sql
import numpy as np
from math import floor

matr_db=sql.connect('zaloha/data.db')
cursor=matr_db.cursor()
cursor.execute('select max(Set1) from results')
r=cursor.fetchall()
cursor.execute('select avg(gTopCor) from results where Accepted = 1 group by Set1')
g=cursor.fetchall()
# r=[i[0] for i in r]
r=r[0][0]
nset = (r)

ksmooth=3
stodch=[]
stdodchpadu=[]
count=0
weight=[]
gfinal=0
sumweight=0
for i in range(r):
    cursor.execute('select gTopCor from results where Accepted = 1 and Set1 = {}'.format(i+1))
    d=cursor.fetchall()
    count+=len(d)
    stodch.append(np.std(d)/np.sqrt(len(d)))
    stodchmod = stodch[-1]*stodch[-1]+ksmooth*ksmooth
    stdodchpadu.append(stodch[-1]*np.sqrt(len(d)))
    weight.append(100/stodchmod)
    sumweight+=weight[-1]
    gfinal+=g[i][0]*weight[-1]
    # print(weight[-1])
#
#
gfinal=gfinal/sumweight

vv=[]
mm=0
for i in range(r):
    vv.append((g[i][0]-gfinal)*np.sqrt(weight[i]))
    mm+=vv[-1]*vv[-1]
    # print(weight[i])

mm=np.sqrt(mm/(nset-1))
gstd=np.std(vv)/np.sqrt(sumweight)

cursor.execute('select gTopCor, Set1 from results where Accepted = 1')
gtop=cursor.fetchall()

normres=[]
for i in gtop:
    normres.append((i[0]-gfinal)/stdodchpadu[i[1]-1]*gstd*np.sqrt(count))
    # print(normres[-1])

# cursor.execute('select Gradient from results where Accepted = 1')
# r=cursor.fetchall()
# r=[i[0]*1000 for i in r]
# np.loadtxt('gradient.csv')
# file=open('gradient.txt', 'r')
# r=file.read().splitlines()
# r[0]=r[0][3:]
#
# r=[float(i) for i in r]



f1='1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500'
f1=f1.split()
f1=[int(i) for i in f1]
print(f1)
#
#
def allan(data, tau):

    res=[]
    for f in tau:
        k=0
        int=[]
        s=0
        it=0
        while k<len(data):
            se=data[k:k+f]

            if len(se)==f:
                int.append(np.mean(se))

                if it>0:
                    s+=(int[-1]-int[-2])**2

            k+=f
            it+=1

        if len(int)-1>0:
            v=np.sqrt(1/(2*(len(int)-1))*s)
            err=v/np.sqrt(floor(len(data)/f))
            res.append([f, v, err])

    return res

a=allan(data=normres, tau=f1)
for i in a:
    print(i)
