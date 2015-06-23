__author__ = 'siteshindra'

import numpy as np
import random
import operator
import collections
from operator import add
import matplotlib.pyplot as plt
import time



def readFile():
    f=open("Ex",'r')
    data=[]
    for line in f:
        lineContent = line.strip('\n')
        lineContent = lineContent.split(",")
        #print lineContent
        data.append([float(i) for i in lineContent[:-1]])
    return data

if __name__ == "__main__":
    dataList = readFile();
    print dataList
    x_list = [x for x in range(0,201,1)]
    print x_list
    #plt.ion()
    #plt.show()

    for steps in dataList:
        plt.plot(x_list,steps)
        plt.pause(0.01)
        plt.clf();
    #plt.show()
    #plt.show()


