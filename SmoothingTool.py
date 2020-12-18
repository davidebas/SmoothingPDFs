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

	# Plotting section
	
	yy_smoothed = np.empty(len(xx))
	
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

	params = {'legend.fontsize': 'x-large',
          'figure.figsize': (8, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
	pylab.rcParams.update(params)

	plt.plot(xx,yy,color='lightskyblue',lw=1)			# plotting original un-smoothed histo
	plt.plot(xx,yy_smoothed,color='red',lw=2)			# plotting smoothed histo
	plt.xlabel('Reconstructed energy [photoelectrons]')
	plt.ylabel('PDF')
	plt.grid()
	plt.legend(l,fontsize=18)
	plt.savefig(OutNamePDF,bbox_inches='tight')

	# Textfile output
	np.savetxt(OutName, yy_smoothed)

	#plt.show()


