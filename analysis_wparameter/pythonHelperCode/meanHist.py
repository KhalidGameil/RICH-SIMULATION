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

file = str(sys.argv[1])
compFile = str(sys.argv[2])
data = pd.read_csv(file, sep=" ", header=None)

true_data = pd.read_csv(compFile,sep=" ", header=None)
print true_data

cols = len(data.columns)

xdata = data.iloc[:,0]
ydata = data.iloc[:,1:1025]
ydata = data.iloc[:,1:1025].values
true_ydata = true_data.iloc[:,1:1025].values
path = file.replace(".txt","")+"_HISTOGRAMS"
if not os.path.exists(path):
        os.mkdir(path)

for i in range(0,1024):
    if np.mean(ydata[:,i])>5E-20:
        fig = plt.figure()
        plt.hist(ydata[:,i])
        plt.axvline(true_ydata[0,i])
        plt.savefig(path+"/"+file.replace(".txt","_")+str(i)+".pdf")
        plt.close()
        print np.mean(ydata[:,i]), true_ydata[0,i],np.std(ydata[:,i]), (np.mean(ydata[:,i])- true_ydata[0,i])/np.std(ydata[:,i])


