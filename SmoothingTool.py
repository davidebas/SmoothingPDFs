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

if __name__=='__main__':

	# filename is the text input file. It must contain the energy bins as the first column.
	# The other columns (you can have as many columns as you want) contains the events for each of the bin spectra energies (first column).


	if(len(sys.argv)<5 or (sys.argv[4]!='SG' and sys.argv[4]!='GF' and sys.argv[4]!='CV')):
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
	Spare = 12345
	
	# Plotting section

	yy_smoothed = np.empty(len(xx))
	res = np.empty(len(xx))
	res_nozero = np.empty(len(xx))

	if(Mode=="SG"):
		if (int(sys.argv[5]) % 2) == 0 :
			print("Smoothing window must be odd. Please retry!")
			quit()
		yy_smoothed = SavitzkyGolay(yy,int(sys.argv[5]),int(sys.argv[6]))
		l = ['Original signal', Mode+'_Window' + sys.argv[5] + '_Pol' + sys.argv[6]]

	if(Mode=="GF"):
		yy_smoothed = GaussianFilter(yy,int(sys.argv[5]))
		l = ['Original signal', Mode+'_Sigma' + sys.argv[5]]

	if(Mode=="CV"):
		if (int(sys.argv[5]) % 2) == 0 :
			print("Smoothing window must be odd. Please retry!")
			quit()
		yy_smoothed = Convolution(yy,int(sys.argv[5]),sys.argv[6])
		l = ['Original signal', Mode+'_Window' + sys.argv[5] + '_Shape-' + sys.argv[6]]


	fig, axs = plt.subplots(3, 1, gridspec_kw={'height_ratios': [3, 1,1]},figsize=(8,10))
	
	
	
	for i in range (0,len(res)):
		if yy[i]==0:
			res[i] = 0
			res_nozero[i]=Spare
		else:
			res[i] = (yy_smoothed[i]-yy[i])/np.sqrt(yy[i])
			res_nozero[i]=res[i]
	
	res_nozero = res_nozero[res_nozero != Spare]
	
	params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
	pylab.rcParams.update(params)
	
	axs[0].plot(xx,yy,color='lightskyblue',lw=1)			# plotting original un-smoothed histo
	axs[0].plot(xx,yy_smoothed,color='red',lw=2)			# plotting smoothed histo
	axs[0].legend(l,fontsize=14)
	axs[0].set(xlabel="Reconstructed energy [photoelectrons]", ylabel="PDF")
	axs[0].grid()

	axs[1].plot(xx,res,color='black',lw=0.001, ms=10)
	axs[1].set(xlabel="Reconstructed energy [photoelectrons]", ylabel="Residuals")
	axs[1].set_ylim(-np.max(np.absolute(res)),np.max(np.absolute(res)))
	axs[1].grid()

	axs[2].hist(res_nozero,color='black',lw=0.001,bins=100)
	axs[2].set(xlabel="Residuals", ylabel="Entries")
	axs[2].set_xlim(-np.max(np.absolute(res)),np.max(np.absolute(res)))
	axs[2].grid()
	l2 = [r'$\mu=$' + str(round(np.mean(res_nozero),4)) + '\n$\sigma=$' + str(round(np.std(res_nozero),4))]
	axs[2].legend(l2,fontsize=14)
	
	print('Smoothing performances: residual histogram')
	print('mu = ' + str(round(np.mean(res_nozero),4)) + '\nsigma = ' + str(round(np.std(res_nozero),4)))

	plt.savefig(OutNamePDF,bbox_inches='tight')

	# Textfile output
	np.savetxt(OutName, yy_smoothed)

	#plt.show()
