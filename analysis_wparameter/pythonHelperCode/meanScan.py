import sys
import os
import pandas as pd
import numpy as np
import ROOT as rt
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate
from scipy import optimize
from scipy.stats import poisson
from scipy.interpolate import UnivariateSpline
from scipy import constants as spc
matplotlib.style.use('ggplot')

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

file = str(sys.argv[1])
data = pd.read_csv(file, sep=" ", header=None)

cols = len(data.columns)

xdata = data.iloc[:,0]
ydata = data.iloc[:,1:1025]

xlines = np.linspace(0.988,1.0,1000)

print xlines
print xdata
print ydata
ydata = data.iloc[:,1:1025].values
xlin = np.linspace(np.min(xdata),np.max(xdata),100)
fig = plt.figure()

plt.plot(xdata,ydata)
#plt.show(block=False)

fig2 = plt.figure()
for i in range(0,1024):
    mean = np.float64(np.max(ydata[:,i]))
    sig = (np.float64(np.max(xdata) - np.min(xdata)))
    if (mean > 1.0E-20 and sig!=0):
        print mean,sig
        print ydata[:,i]

    ##plt.plot(xdata,ydata[:,i],'.')
    ##plt.show()
    x = np.array(xdata.values)
    x = x.astype(np.float64)
    y = np.array(ydata[:,i])
    y = y.astype(np.float64)
    #s = UnivariateSpline(x,y,s=0)
    plt.plot(xlines,s(xlines))
    #g1 = rt.TGraph(len(x),x,y)
    #f = g1.Fit("gaus","QS")
    #a= f.Value(0)
    #x0 = f.Value(1)
    #sigma = f.Value(2)
    #plt.plot(xlin, gaus(xlin,a,x0,sigma))


plt.show()
plt.savefig(file.replace(".txt",".pdf"))
plt.close()


