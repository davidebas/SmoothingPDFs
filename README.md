# SmoothingPDFs

SmoothingTool is a python3 standalone which include algorithms which smooth (or filter) histograms. The main goal is the removing of statistical fluctuations in the histograms, without distorting the signal tendency or peculiarities.

This reveals important in data analysis for a wide range of astroparticle experiments, when the Monte Carlo simulations statistics is lower, or comparable with, the expected statistics on data.

The three implemented algorithms are
- Savitzky Golay filter: a particular type of low-pass filter, well adapted for data smoothing. The core of the method is the successive sub-sets of adjacent data points with a low-degree polynomial by the method of linear least squares. When the data points are equally spaced, an analytical solution to the least-squares equations is applied.
- Gaussian filter: convolution with a Gaussian function
- Convolution method: convolution of a scaled window with the signal. The type of window constitutes the kernel of the convolution method and can be 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'. A flat window produces a simple moving average smoothing.

#### HOW TO RUN:

python3 SmoothingTool.py InputFilename ColumnToBeSmoothed OutputName WhichFilter Additional

InputFilename is the text input file. It must contain the energy bins as the first column. The other columns contains the events for each of the bin spectra energies (first column). You can have as many columns as you want. Look at the ExampleSpectra_JUNO_MC.txt file, which contains spectra from JUNO MC simulations: the first column is the reconstructed energy expressed in photoelectrons; the other columns are the rec. energy spectra from Be7 neutrinos, pep neutrinos, Bi210 decay, Po210 decay, U238 chain decays.

The original and smoothed histograms are plotted in OutputName.pdf , while the latter content is printed also in OutputName.data

1) SavitzkyGolay Filter   => WhichFilter = SG 	2 additional par: WindowLength PolDegree
2) Gaussian Filter 	      => WhichFilter = GF 	1 additional par: Sigma
3) Convolution		        => WhichFilter = CV 	2 additional par: WindowLength Shape

Examples:
1) python3 SmoothingTool.py ExampleSpectra_JUNO_MC.txt 2 Output SG 81 3
2) python3 SmoothingTool.py ExampleSpectra_JUNO_MC.txt 3 Output GF 21 
3) python3 SmoothingTool.py ExampleSpectra_JUNO_MC.txt 4 Output CV 11 flat

WindowLength (SG and CV) must be odd. 
PolDegree (SG only) is the degree of interpolating polynomial.
Sigma (GF only) is the width of the gaussian kernel.

Instructions: python3 SmoothingTool.py

#### MISCELLANEOUS

Convolution method lines have been widely inspired by https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html

Original Savitzky-Golay paper:  "Smoothing and Differentiation of Data by Simplified Least Squares Procedures" 1964, Analytical Chemistry. 36 (8): 1627–39 https://pubs.acs.org/doi/abs/10.1021/ac60214a047
