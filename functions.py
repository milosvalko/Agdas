import json
import matplotlib.pyplot as plt
import numpy as np
from math import floor

def allan(data, tau):
    """Short summary.

    Parameters
    ----------
    data : type
        Description of parameter `data`.
        - data vector
    tau : type
        Description of parameter `tau`.
        - vector with counts of intervals

    Returns
    -------
    type
        Description of returned object.

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
    """Short summary.
        - this function create and save graph to path+project+name+.png file
        - in path+project+name+.csv are source data of graph

    Parameters
    ----------
    x : type
        Description of parameter `x`.
        - vector of x dates x = [x1, ..., xn]
    y : type
        Description of parameter `y`.
        - vector of y dates y = [y1, ..., yn]
    xLabel : type
        Description of parameter `xLabel`.
        - x labels, xLabel=[x1Label, ..., x1Label]
    yLabel : type
        Description of parameter `yLabel`.
        - y labels, yLabel=[y1Label, ..., y1Label]
    title : type
        Description of parameter `title`.
        - title of graph
    path : type
        Description of parameter `path`.
        - path to save graph
    name : type
        Description of parameter `name`.
        - name of graph
    hist : type
        Description of parameter `hist`.
        - boolean hist, create histogram?
    mark : type
        Description of parameter `mark`.
        - list of string, mark = [mark1, ..., markn]
        - e. g. mark = ['b+', 'r*']
    project : type
        Description of parameter `project`.
        - string, prefix of graph name
    columns_name : type
        Description of parameter `columns_name`.
        - list of strings for csv file

    Returns
    -------
    type
        Description of returned object.


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

    #Write graph source data to file
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
    """Short summary.
        - this function round line with writing outputs

    Parameters
    ----------
    l : type
        Description of parameter `l`.
        - line of output files
    index : type
        Description of parameter `index`.
        - list od 2D vectors
        - index = [[i, k], ...]
            - i, position rounded number
            - k, number of decimal places

    Returns
    -------
    type
        Description of returned object.

    """
    for j in index:
        a='{:.'+str(j[1])+ 'f}'
        if l[j[0]]=='-':
            l[j[0]]='-'
        else:
            l[j[0]]=a.format(l[j[0]])

    return l
