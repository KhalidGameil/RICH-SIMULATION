import sys
import os
import pandas as pd
import numpy as np
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
data = pd.read_csv(file, sep=" ", header=None, nrows=35)

cols = len(data.columns)

xdata = data.iloc[0:34,0]
ydata = data.iloc[0:34,1:1025]
ydata = data.iloc[0:34,1:1025].values

print xdata

fig = plt.figure()
plt.plot(xdata,ydata)
plt.xlabel("Beta")
plt.ylabel("Mean Number of Photons")
plt.savefig(file.replace(".txt",".pdf"))
plt.close()


