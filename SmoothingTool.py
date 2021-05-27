# Smoothing tool for solar neutrino spectra
# Davide Basilico - Dec 2020 - davide.basilico@mi.infn.it

import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import matplotlib.pylab as pylab
from scipy import signal
from scipy import interpolate
from scipy import ndimage
import sys
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import CubicSpline
from sklearn.linear_model import LinearRegression
from patsy import cr
from scipy.fft import fft, fftfreq
from scipy.signal.signaltools import wiener

def Convolution(x, window_len, window):

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
       print("Window is none of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'. Please retry")
       quit()

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]

def SavitzkyGolay(yy,Window,Pol):
	yy_smoothed = signal.savgol_filter(yy,Window,Pol)
	return yy_smoothed

def GaussianFilter(yy,Sigma):
	yy_smoothed = ndimage.gaussian_filter1d(yy, Sigma)	
	return yy_smoothed	

def ButterFilter(yy,Window):
	yy_smoothed = scipy.signal.butter(yy,0.4,mysize=Window)
	return yy_smoothed	
	
def SplineFilterOld(xx,yy):
	Spline = UnivariateSpline(xx, yy)
	Spline.set_smoothing_factor(0.5)	
	yy_smoothed = np.array(Spline(xx))
	return yy_smoothed	
	
def SplineFilter(xx,yy):
    x_basis = cr(xx, df=200, constraints="center")
    model = LinearRegression().fit(x_basis, yy)

    # Get estimates

    yy_smoothed = model.predict(x_basis)
    yy_smoothed = np.array(yy_smoothed)	
    print(yy_smoothed)
    return(yy_smoothed)
	
if __name__=='__main__':
	
	# filename is the text input file. It must contain the energy bins as the first column.
	# The other columns (you can have as many columns as you want) contains the events for each of the bin spectra energies (first column).
	

	if(len(sys.argv)<5 or (sys.argv[4]!='WF' and sys.argv[4]!='US' and sys.argv[4]!='SG' and sys.argv[4]!='GF' and sys.argv[4]!='CV')):
		print('python3 SmoothingTool.py [InputFilename] [ColumnToBeSmoothed] [OutputName] [WhichFilter] [Additional]\n')
		print('- SavitzkyGolay Filter   => [WhichFilter] = SG 	2 additional par: [WindowLength] [PolDegree]')
		print('- Gaussian Filter 	 => [WhichFilter] = GF 	1 additional par: [Sigma]')
		print('- Convolution		 => [WhichFilter] = CV 	2 additional par: [WindowLength] [Shape=]')
		print('WindowLength (for SG and CV): the dimension of the smoothing window. It must be odd')
		print('')
		quit()
  		
	filename = sys.argv[1]						# input filename		
	i = sys.argv[2]						# which column (= which species spectrum) do you want to smooth
	spectrum = np.loadtxt(filename,usecols=(0,int(i)))
	OutName = sys.argv[3] + '.data'
	OutNamePDF = sys.argv[3] + '.pdf'				
	Mode = sys.argv[4] 						# which filter
  	
	xx = spectrum[:,0]						# energy binning
	yy = spectrum[:,1]						# counts

#	yy = np.exp(-0.5*(xx-700)*(xx-700)/(100*100)) + 0.005*xx*np.exp(-xx/400)+ 0.02*np.random.randn(xx.size)
#	yy = np.sin(2*3.14/30.*xx) + 0.2*np.random.randn(xx.size)
#	yy = 0.1*np.random.normal(0,1,8001)

	# Plotting section
	
	yy_smoothed = np.empty(len(xx))
	
	if(Mode=="SG"):
		if (int(sys.argv[5]) % 2) == 0 :
			print("Smoothing window must be odd. Please retry!")
			quit()
		yy_smoothed = SavitzkyGolay(yy,int(sys.argv[5]),int(sys.argv[6]))
		
		l = ['Original signal', 'Filtered ' + Mode+'_Window' + sys.argv[5] + '_Pol' + sys.argv[6]]	
	
	if(Mode=="GF"):
		yy_smoothed = GaussianFilter(yy,int(sys.argv[5]))
		l = ['Original signal', 'Filtered ' + Mode+'_Sigma' + sys.argv[5]]			

	if(Mode=="WF"):
		yy_smoothed = WienerFilter(yy,int(sys.argv[5]))
		l = ['Original signal', 'Filtered ' + Mode+'_Window' + sys.argv[5]]			
	
	if(Mode=="CV"):
		if (int(sys.argv[5]) % 2) == 0 :
			print("Smoothing window must be odd. Please retry!")
			quit()
		yy_smoothed = Convolution(yy,int(sys.argv[5]),sys.argv[6])

		l = ['Original signal', 'Filtered ' + Mode+'_Window' + sys.argv[5] + '_Shape-' + sys.argv[6]]
		
	if(Mode=="US"):
		yy_smoothed = SplineFilter(xx,yy)
		l = ['Original signal', 'Filtered ' + Mode+'_k' + sys.argv[5]]								

	params = {'legend.fontsize': 'x-large',
          'figure.figsize': (7, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
	pylab.rcParams.update(params)
	print(type(yy_smoothed))
	
	MovingRMSArrayPostSmoothing = []
	MovingRMSArrayPreSmoothing = []

	for i in range(0,int(xx.size),1):
		#print(np.std(yy_smoothed[i:i+5]))
		MovingRMSArrayPreSmoothing.append(np.std(yy[i:i+5]))
		MovingRMSArrayPostSmoothing.append(np.std(yy_smoothed[i:i+5]))
		#print(MovingRMSArray[i])
	
	MovingRMSArrayPreSmoothing=np.array(MovingRMSArrayPreSmoothing)
	MovingRMSArrayPostSmoothing=np.array(MovingRMSArrayPostSmoothing)
	
	fig, axs = plt.subplots(3,figsize=(8,12))
	
	
	axs[0].plot(xx,yy,color='xkcd:medium blue',lw=0,marker='o',ms=0.5)			# plotting original un-smoothed histo
	axs[0].plot(xx,yy_smoothed,color='xkcd:vermillion',lw=1)			# plotting smoothed histo
	#plt.yscale('log')

#	Residuals = np.abs((yy_smoothed-yy)/(np.sqrt(yy)))
	Residuals = np.abs(yy_smoothed-yy)
	
	axs[1].plot(xx,Residuals,color='xkcd:jade',lw=1)
	axs[1].set_ylim(0,0.0001)
	axs[2].set_ylim(1e-9,1)
		
	#plt.yscale('log')
	
	#axs[3].plot(xx,MovingRMSArrayPreSmoothing,color='xkcd:medium blue',lw=1)	
	#axs[3].plot(xx,MovingRMSArrayPostSmoothing,color='xkcd:vermillion',lw=1)	
	
	#x_Frequency = fftfreq(EndpointPlots+1,EndpointPlots/10.)
	
	samplingFrequency = 1
	
	samplingInterval = 1/samplingFrequency
	
	beginTime = 0
	endTime=len(xx)
	
	time = np.arange(beginTime,endTime,samplingInterval)
	
	y_FrequencyPreSmoothing = np.fft.fft(yy)
	y_FrequencyPostSmoothing = np.fft.fft(yy_smoothed)
	Frequency = np.fft.fftfreq(time.shape[-1])

	#for i in range(5,65,4):
	#	yy_smoothed_2 = GaussianFilter(yy,int(i))
	#	axs[0].plot(xx,yy_smoothed_2,lw=0.5)			# plotting original un-smoothed histo
	

	for i in range(0,2):
		axs[i].set_xlim(0,len(xx))
		axs[i].grid()
		axs[i].set_xlim(0,3000)

	#axs[0].legend(l,fontsize=18)
	axs[2].legend(l,fontsize=18)

	#axs[2].plot(x_Frequency,y_FrequencyPreSmoothing,color='xkcd:medium blue',lw=1)
	#axs[2].plot(x_Frequency,y_FrequencyPostSmoothing,color='xkcd:vermillion',lw=1)

	arr1inds = Frequency.argsort()
	Frequency = Frequency[arr1inds[::-1]]
	y_FrequencyPostSmoothing = y_FrequencyPostSmoothing[arr1inds[::-1]]
	y_FrequencyPreSmoothing = y_FrequencyPreSmoothing[arr1inds[::-1]]

	
	axs[2].plot(Frequency,np.abs(y_FrequencyPreSmoothing),color='xkcd:medium blue',lw=1)
	axs[2].plot(Frequency,np.abs(y_FrequencyPostSmoothing),color='xkcd:vermillion',lw=1)
	#axs[3].set_xlim(0,3000)
	#axs[3].set_yscale('log')
	
	axs[0].set_ylabel('PDFs')
	axs[1].set_ylabel('Difference')
	axs[2].set_ylabel('Amplitude')

	axs[0].set_xlabel('Energy [p.e.]')
	axs[1].set_xlabel('Energy [p.e.]')
	axs[2].set_xlabel('Frequency [1/p.e.]')
	
	#Residuals = np.isfinite(Residuals)
	
	print("High frequency average")
	HighFrequencyAvg = np.log10(np.average(np.abs(y_FrequencyPostSmoothing[0:200])/np.abs(y_FrequencyPreSmoothing[0:200])))
	print(HighFrequencyAvg)

	print("Residuals average")
	ResidualsAverage = np.sum(Residuals[np.isfinite(Residuals)])	
	print(ResidualsAverage)
#	axs[2].plot(abs(y_FrequencyPostSmoothing),color='xkcd:vermillion',lw=1)

#	axs[0].set_yscale('log')

	axs[2].set_yscale('log')

#	axs[1].set_ylim(-0.01,0.01)

	
	#for i in Residuals:
    	#	print(i, end = '\n')
	#axs[2].set_xlim(-100,100)
	axs[2].grid()
	
	#axs[2].set_xlim(-100,100)
	
	plt.savefig(OutNamePDF,bbox_inches='tight')

	# Textfile output
	#np.savetxt(OutName, yy_smoothed)
	#np.savetxt(OutName, sys.argv[5] + ' ' + str(HighFrequencyAvg) + ' ' + str(ResidualsAverage))

	#plt.show()

	with open(OutName, 'a') as f:
		sys.stdout = f # Change the standard output to the file we created.
		print(sys.argv[5] + ' ' + str(HighFrequencyAvg) + ' ' + str(ResidualsAverage))

