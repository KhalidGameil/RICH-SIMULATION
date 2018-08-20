import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate
from scipy.stats import poisson
from scipy.interpolate import UnivariateSpline
from scipy import constants as spc
import sys
matplotlib.style.use('ggplot')

def chunkyMeanPlot(data,temp):
    data = data.sort('beta')
    beta = np.array(data['beta'])
    ll = np.array(data['ll'])
    dupStop = False

    betaMean = []
    betaMean.append(beta[0])
    llMean = []
    llMean.append(ll[0])
    llm = ll[0]
    m = []
    meanCount = 1
    m.append(meanCount)
    for i in range(1,len(ll)):
        if beta[i-1]!=beta[i]:
            dupStop = True
            betaMean.append(beta[i])
            llMean.append(llm/meanCount)
            m.append(meanCount)
            llm =ll[i];
            meanCount =1;
        else:
            llm+=ll[i]
            meanCount+=1

    fig = plt.figure()
    betaMean = np.array(betaMean)
    llMean = np.array(llMean)
    print betaMean
    print llMean
    print m
    plt.plot(betaMean,llMean)
    fig.savefig(temp+"BETAMEAN.pdf")
    plt.close()


# Get the total number of args passed to the demo.py
total = len(sys.argv)
 
print total
files = [f for f in os.listdir('.') if os.path.isfile(f)]
pionMass =  137.2735 #MeV
if total == 2:
    files = [str(sys.argv[1])] #get second argument

for temp in files:
    if ".py" not in temp and ".pdf" not in temp and "betaMeanScan.txt" not in temp:
        print temp
        path = temp.replace(".txt",'')
 	path = temp.replace("/","_")
        path = temp.replace(".","")
	os.mkdir(path)
        file = open(temp)
        line = file.readlines()
        line = line[0]
        parameterArray = np.fromstring(line, dtype=float, sep=" ")
        parameterArray[0] = parameterArray[0]/np.sqrt(parameterArray[0]**2 + pionMass**2)
        print parameterArray
        data = pd.read_csv(temp, sep=" ",skiprows=1, header=None)
        data.columns = ['i','x','errX','y','errY','theta','thetaErr','phi','phiErr','beta','betaErr','msg1','msg2','msg3','msg4','ll','throw']
        data = data.sort('beta')
        data = data.drop(labels='throw', axis=1)
        chunkyMeanPlot(data,path+"/"+path)
        parameterIndex = [1,2,3,4,0]
        warning = data['msg1']+data['msg2']+data['msg3']+data['msg4']
        for i in range(0,5):
            fig, ax = plt.subplots(1)
            d = (np.array(data[[2*i+1]])-parameterArray[parameterIndex[i]])/parameterArray[parameterIndex[i]]
            if (total ==1):
                d = d[warning == 3]
            elif (total == 2):
                d = d[d>-0.00005]
            betaError = data['betaErr']
            print "error Correct", betaError[warning==3]
            print "error InCorrect",betaError[warning!=3]
            #print np.array(data[[2*i+1]]),parameterArray[parameterIndex[i]],np.array(data[[2*i+2]]), d
            plt.hist(d,50)
            plt.xlabel(data.columns[2*i+1])
            mu = d.mean()
            sigma = d.std()
            #plt.title(r'$\mathrm{Histogram\ of\ Parameter\ Fractional\ Error')
            textstr = '$\mu=%.6f$\n$\sigma=%.6f$\n$Entries=%d$'%(mu, sigma,len(d))
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
                    verticalalignment='top', bbox=props)
            plt.savefig(path+"/"+path+data.columns[2*i+1]+".pdf")
            plt.close()

        beta = np.array(data['beta'])
        ll = np.array(data['ll'])
        if (total ==1):
            beta = beta[warning == 3] #[np.where(data['msg1']==0 and data['msg2']==0 and data['msg3']==0 and data['msg4']==3)]
            ll = ll[warning == 3]#[np.where(data['msg1']==0 and data['msg2']==0 and data['msg3']==0 and data['msg4']==3)]

        fig, ax = plt.subplots(1)
        d = np.array(beta)
        plt.hist(d,50)
        mu = d.mean(dtype=np.float64)
        print mu
        sigma = d.std(dtype=np.float64)
        print sigma
        textstr = '$\mu_{True}=%.6f$\n$\mu=%.6f$\n$\sigma=%g$\n$Entries=%d$'%(parameterArray[0],mu, sigma,len(beta))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        plt.savefig(path+"/"+path+"BETA_MINIMIZED_VALUES.pdf")

        fig, ax = plt.subplots(1)
        plt.plot(beta,ll,'o')
        plt.xlabel('Beta')
        plt.ylabel('ll')
        plt.savefig(path+"/"+path+"LL_v_Beta.pdf")
        plt.close()


