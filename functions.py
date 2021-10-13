import json
import matplotlib.pyplot as plt
import numpy as np
from math import floor

def allan(data, tau):
    """
    Function return allans deviation
    """

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
            res.append([v, err])

    return res

def printDict(dict):
    """
    Print dictionary like json format
    """
    print(json.dumps(dict,indent=5))


def graph(x, y, xLabel, yLabel, title, path, name, hist, mark, project, columns_name):
    """
    Create graph and save him to 'path' folder
    """
    gr=plt
    if hist == True:
        gr.hist(y,edgecolor='black')

    if hist == False:
        for i in range(len(x)):
            gr.plot(x[i],y[i], mark[i])


    gr.xlabel(xLabel)
    gr.ylabel(yLabel)
    gr.title(title)
    gr.savefig(path+'/'+project+'_'+name)
    gr.close()

    d=open(path+'/'+project+'_'+name+'.csv', 'w')
    h=''
    for i in columns_name:
        h+=i+'_x'
        h+=';'
        h+=i+'_y'
        h+=';'
    d.write(h+'\n')
    n=max([len(i) for i in x])

    for i in range(n):
        l=[]
        for j in range(len(x)):

            try:
                l.append(x[j][i])
            except IndexError:
                x[j].append('-')
                l.append(x[j][i])

            try:
                l.append(y[j][i])
            except IndexError:
                y[j].append('-')
                l.append(y[j][i])

        ll=''
        for k in l:
            ll+=str(k)
            ll+=';'
        ll+='\n'
        d.write(ll)

    d.close()

def roundList(l, index):
    for j in index:
        a='{:.'+str(j[1])+ 'f}'
        l[j[0]]=a.format(l[j[0]])

    return l

# def allan(data, tau):
#
#
#     res=[]
#     for f in tau:
#         k=0
#         int=[]
#         s=0
#         it=0
#         while k<len(data):
#             se=data[k:k+f]
#
#             if len(se)==f:
#                 int.append(np.mean(se))
#
#                 if it>0:
#                     s+=(int[-1]-int[-2])**2
#
#             k+=f
#             it+=1
#
#         if len(int)-1>0:
#             v=np.sqrt(1/(2*(len(int)-1))*s)
#             err=v/np.sqrt(floor(len(data)/f))
#             res.append([f, v, err])
#
#     return res
