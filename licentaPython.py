import numpy.fft as fft
import numpy as np
import math 
from scipy import signal
from matplotlib import pyplot as plt
from obspy import read
from obspy.signal.filter import envelope

def LoadData(file,delimiter):
    data = np.genfromtxt(file,delimiter=delimiter)
    acc_x = data[:, 0]
    acc_y = data[:, 1]
    acc_z = data[:, 2]
    
    return acc_x, acc_y, acc_z

def UtilsData():
    for i in range(samples):
        if (x[i] < -260 or x[i] > -240) and t[i] > 15.000:
            margine_inf = math.floor(t[i]/0.01)
            break
            
    for i in range(samples):
        if (x[i] < -250 or x[i] > -260) and t[i] > 25.000:
            margine_sup = math.floor(t[i]/0.01)
            break  
    
    dateUtile_x = np.zeros(margine_sup - margine_inf + 1)
    dateUtile_y = np.zeros(margine_sup - margine_inf + 1)
    dateUtile_z = np.zeros(margine_sup - margine_inf + 1)
    j = 0
    for i in range(margine_inf, margine_sup):
        dateUtile_x[j] = x[i]
        dateUtile_y[j] = y[i]
        dateUtile_z[j] = z[i]
        j = j + 1
        
    return dateUtile_x, dateUtile_y, dateUtile_z
  
def FindSamples(x):
    samples = len(x);
    t = np.zeros(samples)
    dF = Fs/samples;
    for i in range(samples):
        t[i] = i*0.01;
    return t
    
def CenteredData(semnal):
    R = 0.9
    samples = len(semnal)
    semnalCentrat = np.zeros(samples);
    for i in range(1, samples):
        semnalCentrat[i] = semnal[i] - semnal[i-1] + R*semnalCentrat[i-1]
      
    return semnalCentrat
    
def filtrareAnvelopa(semnal, dimFereastra):
    stop = len(semnal) - dimFereastra
    semnalFiltrat = semnal
    for i in range(0,stop):
        R = 0
        semnalDeAjutor = semnal[i:i+dimFereastra]
        for j in range(0,dimFereastra):
            R = R + semnalDeAjutor[j]**2
        R = math.sqrt(R/dimFereastra)
        semnalFiltrat[i] = R
        
    for i in range(len(semnal)-dimFereastra,len(semnal)-1):
        R = 0
        dimFereastra = dimFereastra - 1
        semnalDeAjutor = semnalFiltrat[i:i+dimFereastra]
        for j in range(0,dimFereastra):
            R = R + semnalDeAjutor[j]**2
        R = math.sqrt(R/dimFereastra)
        semnalFiltrat[len(semnal)-dimFereastra-1] = R
        
    return semnalFiltrat
    
def fftFunction(semnalFiltrat):
    semnalFFT = fft.fftshift(fft.fft(semnalFiltrat))*2/samples
    semnalFFT[round(samples2/2)] = 0
    semnalFFT[round(samples2/2)+1] = 0
    semnalFFT[round(samples2/2)-1] = 0
    return semnalFFT
      
x,y,z = LoadData("./fierastrau_achizitierapida.csv",delimiter=',')

Fs = 0.19
samples = len(x)
t = FindSamples(x)
    
#plt.subplot(3,4,1),plt.plot(t, x)
#plt.subplot(3,4,5),plt.plot(t, y)
#plt.subplot(3,4,9),plt.plot(t, z)
    
x_util,y_util,z_util = UtilsData()

samples2 = len(x_util)
t2 = FindSamples(x_util)    

#plt.subplot(3,4,2),plt.plot(t2, x_util)
#plt.subplot(3,4,6),plt.plot(t2, y_util)
#plt.subplot(3,4,10),plt.plot(t2, z_util)

semnalCentrat_x = CenteredData(x_util)
semnalCentrat_y = CenteredData(y_util)
semnalCentrat_z = CenteredData(z_util)

#plt.subplot(3,4,3),plt.plot(t2, semnalCentrat_x)
#plt.subplot(3,4,7),plt.plot(t2, semnalCentrat_y)
#plt.subplot(3,4,11),plt.plot(t2, semnalCentrat_z)


upper_x = envelope(semnalCentrat_x)
upper_y = envelope(semnalCentrat_y)
upper_z = envelope(semnalCentrat_z)
#plt.subplot(3,4,4),plt.plot(t2,upper_x)


filteredSignal = filtrareAnvelopa(upper_x, 15)
#plt.subplot(3,4,8), plt.plot(t2, filteredSignal)


acc_x_fft = fftFunction(filteredSignal)
dF2 = Fs/samples2
f_x = np.arange(-Fs/2, Fs/2, dF2)

plt.plot(f_x, abs(acc_x_fft))
plt.show()