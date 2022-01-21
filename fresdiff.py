# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 11:03:41 2021

@author: Haadi
"""

import scipy.odr as odr
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import statistics as st

# #############################################################################
# 635nm laser, 0.04 to 0.36 mm slits, both m = 1 and m = -1 used
# #############################################################################

# # Calculate SDOM of the experimental data and plot it

#initialize array with various slith widths represented by variable b

b = np.array([0.04, 0.08, 0.12, 0.16, 0.24, 0.36])
x = 1/b

data_1 = [0.014114, 0.019496, 0.018988, 0.019453] # 4 experimental repetitions for b = 0.04mm
data_2 = [0.009316, 0.008525, 0.008620, 0.009128, 0.007948, 0.008740, 0.009369, 0.008930] # 4 experimental repetitions for b = 0.08mm
data_3 = [0.006372, 0.005597, 0.006932, 0.005666, 0.006674, 0.006329, 0.006691, 0.005425] # 4 experimental repetitions for b = 0.12mm
data_4 = [0.004409, 0.003832, 0.004366, 0.004039, 0.003927, 0.004478, 0.004099, 0.004263] # 4 experimental repetitions for b = 0.16mm
data_5 = [0.002781, 0.002678, 0.003160, 0.002876, 0.002833, 0.002687, 0.002471, 0.002635] # 4 experimental repetitions for b = 0.24mm
data_6 = [0.002773, 0.000741, 0.000749, 0.000741, 0.001361, 0.001438, 0.001498, 0.001386] # 4 experimental repetitions for b = 0.36mm
  
#for data_2 through data_6, the diffraction agle is calculated for m=1 and m=-1, hence creating 8 data points
#this could not be done for data_1 as the intensity distribution was too wide to identify m = 1

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
    
theta_mean= np.array([Mean(data_1),Mean(data_2),Mean(data_3),Mean(data_4),Mean(data_5),Mean(data_6)])
theta_sdom= np.array([SDOM(data_1), SDOM(data_2), SDOM(data_3), SDOM(data_4), SDOM(data_5), SDOM(data_6)])

plt.errorbar(x, theta_mean, theta_sdom,fmt='.k')
plt.title('Graph of the experimental data and SDOM')
plt.xlabel('1/Slit Width (1/mm)')
plt.ylabel(' θ (Radians)')
plt.show()

# zoomxmin = ; zoomxmax = ; zoomymin = ; zoomymax = ;
# zoomwidth = zoomxmax - zoomxmin; zoomheight = zoomymax - zoomymin
# plt.gca().add_patch(mpl.patches.Rectangle((zoomxmin,zoomymin),zoomwidth,zoomheight, fill=False, linestyle='dashed'))
# plt.axes([0.25,0.6, 0.20,0.20])
# plt.errorbar(x,wav_min,wav_sdom,fmt='k.')

# plt.xlim(zoomxmin,zoomxmax)
# plt.ylim(zoomymin,zoomymax)
# plt.show()

# Fit the data using f = m*x + b and add theoretical model 

xplot=np.arange(0,35, 35/6)

lamb_theoretical = 635*10**(-6)*xplot
plt.plot(xplot, lamb_theoretical, label='Theoretical prediction')
plt.xlim(2, 27)
plt.ylim(0.0025, 0.02)
plt.title('Diffraction experiment for the Fresnel Relation')
plt.xlabel('1/Slit Width (1/mm)')
plt.ylabel(' θ (Radians) ')
plt.legend()
plt.show

y = np.array(theta_mean)
sig_y = np.array(theta_sdom)
err_x = 1/2 * abs( (1/(b+0.005)) - (1/b-0.005))

A_start=0.5*10**-4
B_start=1.5*10**-9

def f(B, x):
    print(B[1])
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
plt.plot(xplot ,par_best[0] + par_best[1]*xplot,'r-',label='Line of best fit')
plt.xlim(0, 27)
plt.ylim(0, 0.02)
plt.legend(loc='lower right')
plt.show()

# Uncertainty Propagation:
    
theta_uncert = 1/2* abs( (theta_mean + theta_sdom)*b - (theta_mean - theta_sdom)*b )
b_uncert = 1/2* abs( theta_mean*(b + 0.005) - theta_mean*(b - 0.005) )

lamb_uncert = ((theta_uncert)**2 + (b_uncert)**2)**(0.5)

#print(lamb_uncert)
