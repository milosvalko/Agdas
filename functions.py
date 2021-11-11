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
