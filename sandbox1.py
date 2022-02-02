import sqlite3 as sql
import numpy as np
from time import time
import matplotlib.pyplot
import matplotlib as mpb
mpb.use('ps')

frmax = 9400
frmin = 150
siz = 15
frmaxplot = 10400

d = sql.connect(r'c:\Users\Jakub\Desktop\pecny\data\x153_glb\res\data.db')
c = d.cursor()
c.execute('select Accepted from results')
r = c.fetchall()

res = np.loadtxt('resm.csv', delimiter = ';')

t=time()
p = mpb.pyplot
p.rcParams['figure.figsize']=(20,13)
p.ylim([-siz, siz])
p.plot([frmin, frmin], [-siz, siz], '-b', lw = 1)
p.plot([frmax, frmax], [-siz, siz], '-b', lw = 1)
p.title('Residuals for all drops')
p.ylabel('Residuals /nm')
p.xlabel('Fringe #')
# x=range(frmaxplot)
# for i in range(len(r)):
#     acc = r[i][0]
#
#     if acc:
#         # p.plot(x, res[i, :], '-k', lw=0.2)
#         maxi = np.max(res)
#     # else:
#     #     p.plot(x, res[i, :], '-', color = "0.5", lw = 0.2)
max_env = []
min_env = []
for i in range(res.shape[1]):

    max_env.append(np.max(res[:, i]))
    min_env.append(np.min(res[:, i]))

p.plot(max_env, '.k')
p.plot(min_env, '.k')
p.savefig('pokus.png')
p.close()
print(time()-t)
