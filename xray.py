import scipy.odr as odr
import scipy.stats
from scipy import constants as cst
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statistics as st
import math

# #############################################################################
# 1 mA, 21kV to 35kV, from 2.5 to 20 degrees
# #############################################################################

# # Calculate SDOM of the experimental data and plot it

U = np.arange(21,37,2)
x = 1/U

data_1 = [57.5*10**(-9), 56.3*10**(-9), 56.9*10**(-9), 57.3*10**(-9)] # 4 experimental repetitions for n=1
data_2 = [52.9*10**(-9), 52.9*10**(-9), 51.4*10**(-9), 50.9*10**(-9)] # 4 experimental repetitions for n=2
data_3 = [47.2*10**(-9), 48.3*10**(-9), 48.5*10**(-9), 47.9*10**(-9)] # ...
data_4 = [44.9*10**(-9), 43.7*10**(-9), 43.1*10**(-9), 44.5*10**(-9)]
data_5 = [41.4*10**(-9), 40.3*10**(-9), 41.8*10**(-9), 40.4*10**(-9)]
data_6 = [38.0*10**(-9), 38.2*10**(-9), 37.4*10**(-9), 37.9*10**(-9)]
data_7 = [35.7*10**(-9), 36.8*10**(-9), 36.0*10**(-9), 35.3*10**(-9)]
data_8 = [34.5*10**(-9), 34.5*10**(-9), 35.7*10**(-9), 35.1*10**(-9)]
  
def Mean(lst):
    """
    Returns the mean of a list of numbers.
    """
    return st.mean(lst)

def SDOM (lst):
    """
    Returns the standard deviation of the mean (SDOM).
    """
    return scipy.stats.sem(lst)
    
wav_min=[Mean(data_1),Mean(data_2),Mean(data_3),Mean(data_4),Mean(data_5),Mean(data_6),Mean(data_7),Mean(data_8)]
wav_sdom=[SDOM(data_1)+2.5*10**(-10),SDOM(data_2)+2.5*10**(-10),SDOM(data_3)+2.5*10**(-10),SDOM(data_4)+2.5*10**(-10),SDOM(data_5)+2.5*10**(-10),SDOM(data_6)+2.5*10**(-10),SDOM(data_7)+2.5*10**(-10),SDOM(data_8)+2.5*10**(-10)]


plt.errorbar(x ,wav_min, wav_sdom,fmt='.k')
plt.xlim(0.025,0.05)
plt.ylim(3.4*(10**(-8)),6.0*(10**(-8)))
plt.title('Graph of the experimental data and SDOM')
plt.xlabel('1/Accelerating Potential')
plt.ylabel('Minimum Wavelength')


zoomxmin = 0.0315; zoomxmax = 0.0330; zoomymin = 3.6*10**(-8); zoomymax = 4*10**(-8);
zoomwidth = zoomxmax - zoomxmin; zoomheight = zoomymax - zoomymin
plt.gca().add_patch(mpl.patches.Rectangle((zoomxmin,zoomymin),zoomwidth,zoomheight, fill=False, linestyle='dashed'))
plt.axes([0.25,0.6, 0.20,0.20])
plt.errorbar(x,wav_min,wav_sdom,fmt='k.')

plt.xlim(zoomxmin,zoomxmax)
plt.ylim(zoomymin,zoomymax)
plt.show()

# Fit the data using f = m*x + b and add theoretical data
xplot=np.arange(0.025,0.050,0.025/14)


wav_theoretical = ((cst.h * cst.c)/(cst.e))*(xplot)
plt.plot(xplot, wav_theoretical, label='Theoretical predictions')
plt.title('X-ray experiment for the Duane-Hunt relation')
plt.xlabel('1/Accelerating Voltage (1/kV)')
plt.ylabel('Minimum Wavelength (m)')
plt.legend()
plt.show

y = np.array(wav_min)
sig_y = np.array(wav_sdom)
err_x = np.array([0.0002,0.00018,0.00016,0.00014,0.00012,0.00010,0.00008,0.00008])

A_start=0.5*10**-4
B_start=1.5*10**-9

def f(B, x):
    return B[0] + B[1]*x

odr_model = odr.Model(f)
odr_data = odr.RealData(x,y)
odr_obj = odr.ODR(odr_data,odr_model,beta0=[A_start,B_start])
odr_res = odr_obj.run()
par_best = odr_res.beta
par_sig_ext = odr_res.sd_beta
par_cov = odr_res.cov_beta 
print(" The (INTERNAL!) covariance matrix  = \n", par_cov)
 
chi2 = odr_res.sum_square
print("\n Chi-squared = ", chi2)
chi2red = odr_res.res_var
print(" Reduced chi-squared = ", chi2red, "\n")

odr_res.pprint()

plt.errorbar(x,y,xerr=err_x,yerr=sig_y,fmt='.k',label='Experimental data and SDOM')
plt.plot(xplot,par_best[0] + par_best[1]*xplot,'r-',label='Line of best fit')
plt.xlim(0.025,0.05)
plt.ylim(3.4*(10**(-8)),6.0*(10**(-8)))
plt.legend(loc='lower right')
plt.show()
