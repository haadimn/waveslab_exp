# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 14:09:12 2021

@author: Haadi
"""

import scipy.odr as odr
import scipy.stats 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statistics as st
import math

# #############################################################################
# # Copper 0,4 mm: 150g (4 repetitions of the experiment)
# #############################################################################

# # Calculate SDOM of the experimental data and plot it

x = np.array([1,2,3,4,5,6,7,8])

data_1 = 15.96,17.05,12.97,11.05 # 4 experimental repetitions for n=1
data_2 = 32.80,30.20,33.98,31.23 # 4 experimental repetitions for n=2
data_3 = 48.98,44.21,48.95,45.01 # ...
data_4 = 63.16,65.09,61.29,60.01
data_5 = 80.44,79.43,83.09,78.96
data_6 = 94.64,92.56,91.67,91.05
data_7 = 113.8,115.01,110.03,111.96
data_8 = 129.2,128.79,132.67,127.01
  
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
    
f_exp150=[Mean(data_1),Mean(data_2),Mean(data_3),Mean(data_4),Mean(data_5),Mean(data_6),Mean(data_7),Mean(data_8)]
f_exp150_sdom=[SDOM(data_1),SDOM(data_2),SDOM(data_3),SDOM(data_4),SDOM(data_5),SDOM(data_6),SDOM(data_7),SDOM(data_8)]

plt.errorbar(x,f_exp150,f_exp150_sdom,fmt='k.')
plt.title('Graph of the experimental data and SDOM (r=0.44mm, m=150g)')
plt.xlabel('Harmonic')
plt.ylabel('Experimental frequency')


zoomxmin = 2.85; zoomxmax = 3.15; zoomymin = 45; zoomymax = 52;
zoomwidth = zoomxmax - zoomxmin; zoomheight = zoomymax - zoomymin
plt.gca().add_patch(mpl.patches.Rectangle((zoomxmin,zoomymin),zoomwidth,zoomheight, fill=False, linestyle='dashed'))
plt.axes([0.2,0.6, 0.20,0.20])
plt.errorbar(x,f_exp150,f_exp150_sdom,fmt='k.')
plt.xlim(zoomxmin,zoomxmax)
plt.ylim(zoomymin,zoomymax)
plt.show()

# Fit the data using f = m*x + b and add theoretical data

f_theoretical_150=(1/(2*1.482))*(math.sqrt((0.150*9.81/0.00136)))*x
plt.plot(x,f_theoretical_150,label='Theoretical predictions')

plt.errorbar(x,f_exp150,f_exp150_sdom,fmt='k.', label='Experimental data and SDOM')
plt.title('Melde experiment for lower orders of n ignoring dispersion')
plt.xlabel('Harmonic')
plt.ylabel('Experimental frequency')


y = np.array(f_exp150)
sig_y = np.array(f_exp150_sdom)
A_start=0.5
B_start=1.5

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

xplot=np.arange(1,9,1)
plt.errorbar(x,y,yerr=sig_y,fmt='k.')
plt.plot(xplot,par_best[0] + par_best[1]*xplot,'r-',label='Line of best fit')
plt.legend(loc='lower right')


zoomxmin = 2.85; zoomxmax = 3.15; zoomymin = 45; zoomymax = 52;
zoomwidth = zoomxmax - zoomxmin; zoomheight = zoomymax - zoomymin
plt.gca().add_patch(mpl.patches.Rectangle((zoomxmin,zoomymin),zoomwidth,zoomheight, fill=False, linestyle='dashed'))
plt.axes([0.2,0.6, 0.20,0.20])
plt.errorbar(x,f_exp150,f_exp150_sdom,fmt='k.')
plt.xlim(zoomxmin,zoomxmax)
plt.ylim(zoomymin,zoomymax)

plt.errorbar(x,y,yerr=sig_y,fmt='k.')
plt.plot(xplot,par_best[0] + par_best[1]*xplot,'r-',label='Line of best fit')

plt.show()

# #############################################################################
# # Copper 0,4 mm: 200g (2 repetitions of the experiment)
# #############################################################################

data_200_1 = 13.90,10.55 # 2 experimental repetitions for n=1
data_200_2 = 29.87,26.01 # 2 experimental repetitions for n=2
data_200_3 = 43.70,39.92 # ...
data_200_4 = 53.71,49.83
data_200_5 = 68.47,65.01
data_200_6 = 78.21,82.29
data_200_7 = 93.92,97.25
data_200_8 = 106.8,109.8
data_200_9 = 124.6,128.2
data_200_10 = 141.8,138.0
data_200_11 = 161.9,165.0

x = np.array([1,2,3,4,5,6,7,8,9,10,11])

f_exp200=[Mean(data_200_1),Mean(data_200_2),Mean(data_200_3),Mean(data_200_4),Mean(data_200_5),Mean(data_200_6),Mean(data_200_7),Mean(data_200_8),Mean(data_200_9),Mean(data_200_10),Mean(data_200_11)]
f_exp200_sdom=[SDOM(data_200_1),SDOM(data_200_2),SDOM(data_200_3),SDOM(data_200_4),SDOM(data_200_5),SDOM(data_200_6),SDOM(data_200_7),SDOM(data_200_8),SDOM(data_200_9),SDOM(data_200_10),SDOM(data_200_11)]

plt.errorbar(x,f_exp200,f_exp200_sdom,fmt='k.')
plt.title('Graph of the experimental data and SDOM (r=0.44mm, m=200g)')
plt.xlabel('Harmonic')
plt.ylabel('Experimental frequency')

zoomxmin = 2.85; zoomxmax = 3.15; zoomymin = 37; zoomymax = 48;
zoomwidth = zoomxmax - zoomxmin; zoomheight = zoomymax - zoomymin
plt.gca().add_patch(mpl.patches.Rectangle((zoomxmin,zoomymin),zoomwidth,zoomheight, fill=False, linestyle='dashed'))
plt.axes([0.2,0.6, 0.20,0.20])
plt.errorbar(x,f_exp200,f_exp200_sdom,fmt='k.')
plt.xlim(zoomxmin,zoomxmax)
plt.ylim(zoomymin,zoomymax)
plt.show()

# Fit the data using f = m*x + b and add theoretical data

f_theoretical_200=((x/(2*1.482))*(math.sqrt((0.200*9.81/0.00136))))*(1+(x**2)*((np.pi**3*0.00044**4*130*10**9)/(8*9.81*1.482**2*0.2)))
plt.plot(x,f_theoretical_200,label='Theoretical predictions')

plt.errorbar(x,f_exp200,f_exp200_sdom,fmt='k.', label='Experimental data and SDOM')
plt.title('Melde experiment for high orders of n')
plt.xlabel('Harmonic')
plt.ylabel('Experimental frequency')

y = np.array(f_exp200)
sig_y = np.array(f_exp200_sdom)
A_start=0.5
B_start=1.5

def f(B, x):
    return 0.038*B[0]*x**3 + B[1]*x

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

xplot=np.arange(1,12,1)
plt.errorbar(x,y,yerr=sig_y,fmt='k.')
plt.plot(xplot,par_best[0] + par_best[1]*xplot,'r-',label='Line of best fit')

plt.legend(loc='upper left')
plt.show()

#############################################################################
# Copper 0,22 mm: 200g (2 repetitions of the experiment)
#############################################################################

data_200_1_sm = 26.47,24.32 # 2 experimental repetitions for n=1
data_200_2_sm = 55.41,57.22 # 2 experimental repetitions for n=2
data_200_3_sm = 83.24,86.91 # ...
data_200_4_sm = 110.1,107.3
data_200_5_sm = 134.4,137.2
data_200_6_sm = 166.4,163.3
data_200_7_sm = 193.7,190.8
data_200_8_sm = 211.2,209.0

x = np.array([1,2,3,4,5,6,7,8])

f_exp200_sm=[Mean(data_200_1_sm),Mean(data_200_2_sm),Mean(data_200_3_sm),Mean(data_200_4_sm),Mean(data_200_5_sm),Mean(data_200_6_sm),Mean(data_200_7_sm),Mean(data_200_8_sm)]
f_exp200_sdom_sm=[SDOM(data_200_1_sm),SDOM(data_200_2_sm),SDOM(data_200_3_sm),SDOM(data_200_4_sm),SDOM(data_200_5_sm),SDOM(data_200_6_sm),SDOM(data_200_7_sm),SDOM(data_200_8_sm)]

plt.errorbar(x,f_exp200_sm,f_exp200_sdom_sm,fmt='k.')
plt.title('Graph of the experimental data and SDOM (r=0.22mm, m=200g)')
plt.xlabel('Harmonic')
plt.ylabel('Experimental frequency')

zoomxmin = 2.85; zoomxmax = 3.15; zoomymin = 80; zoomymax = 90
zoomwidth = zoomxmax - zoomxmin; zoomheight = zoomymax - zoomymin
plt.gca().add_patch(mpl.patches.Rectangle((zoomxmin,zoomymin),zoomwidth,zoomheight, fill=False, linestyle='dashed'))
plt.axes([0.2,0.6, 0.20,0.20])
plt.errorbar(x,f_exp200_sm,f_exp200_sdom_sm,fmt='k.')
plt.xlim(zoomxmin,zoomxmax)
plt.ylim(zoomymin,zoomymax)
plt.show()

# Fit the data and add theoretical predictions

f_theoretical_200_sm=((x/(2*1.482))*(math.sqrt((0.200*9.81/(np.pi*0.00022**2*8960)))))*(1+(x**2)*((np.pi**3*0.00022**4*130*10**9)/(8*9.81*1.482**2*0.2)))
plt.plot(x,f_theoretical_200_sm,label='Theoretical predictions')

plt.errorbar(x,f_exp200_sm,f_exp200_sdom_sm,fmt='k.', label='Experimental data and SDOM')
plt.title('Melde experiment for high orders of n (r=0.22mm)')
plt.xlabel('Harmonic')
plt.ylabel('Experimental frequency')

y = np.array(f_exp200_sm)
sig_y = np.array(f_exp200_sdom_sm)
A_start= 1
B_start=10

def f(B, x):
    return 0.002*B[0]*x**3 + B[1]*x

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

xplot=np.arange(1,9,1)
plt.errorbar(x,y,yerr=sig_y,fmt='k.')
plt.plot(xplot,par_best[0] + par_best[1]*xplot,'r-',label='Line of best fit')

plt.legend(loc='upper left')
plt.show()

#############################################################################
# Copper 0,2 mm: 250g (2 repetitions of the experiment)
#############################################################################

data_250_1_sm = 35.61,33.21 # 2 experimental repetitions for n=1
data_250_2_sm = 62.29,59.82 # 2 experimental repetitions for n=2
data_250_3_sm = 94.04,96.28 # ...
data_250_4_sm = 124.2,121.3
data_250_5_sm = 156.4,152.7
data_250_6_sm = 191.2,189.0
data_250_7_sm = 220.4,223.4
data_250_8_sm = 254.1,251.9
data_250_9_sm = 283.5,281.4

x = np.array([1,2,3,4,5,6,7,8,9])

f_exp250_sm=[Mean(data_250_1_sm),Mean(data_250_2_sm),Mean(data_250_3_sm),Mean(data_250_4_sm),Mean(data_250_5_sm),Mean(data_250_6_sm),Mean(data_250_7_sm),Mean(data_250_8_sm),Mean(data_250_9_sm)]
f_exp250_sdom_sm=[SDOM(data_250_1_sm),SDOM(data_250_2_sm),SDOM(data_250_3_sm),SDOM(data_250_4_sm),SDOM(data_250_5_sm),SDOM(data_250_6_sm),SDOM(data_250_7_sm),SDOM(data_250_8_sm),SDOM(data_250_9_sm)]

plt.errorbar(x,f_exp250_sm,f_exp250_sdom_sm,fmt='k.')
plt.title('Graph of the experimental data and SDOM (r=0.22mm, m=250g)')
plt.xlabel('Harmonic')
plt.ylabel('Experimental frequency')

zoomxmin = 2.85; zoomxmax = 3.15; zoomymin = 90; zoomymax = 99;
zoomwidth = zoomxmax - zoomxmin; zoomheight = zoomymax - zoomymin
plt.gca().add_patch(mpl.patches.Rectangle((zoomxmin,zoomymin),zoomwidth,zoomheight, fill=False, linestyle='dashed'))
plt.axes([0.2,0.6, 0.20,0.20])
plt.errorbar(x,f_exp250_sm,f_exp250_sdom_sm,fmt='k.')
plt.xlim(zoomxmin,zoomxmax)
plt.ylim(zoomymin,zoomymax)
plt.show()

# Fit the data and add theoretical predictions

f_theoretical_250_sm=((x/(2*1.482))*(math.sqrt((0.250*9.81/(np.pi*0.00022**2*8960)))))*(1+(x**2)*((np.pi**3*0.00022**4*130*10**9)/(8*9.81*1.482**2*0.25)))
plt.plot(x,f_theoretical_250_sm,label='Theoretical predictions')

plt.errorbar(x,f_exp250_sm,f_exp250_sdom_sm,fmt='k.', label='Experimental data and SDOM')
plt.title('Melde experiment for high orders of n (r=0.2mm)')
plt.xlabel('Harmonic')
plt.ylabel('Experimental frequency')

y = np.array(f_exp250_sm)
sig_y = np.array(f_exp250_sdom_sm)
A_start= 1
B_start=10

def f(B, x):
    return 0.03*B[0]*x**2 + B[1]*x

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

xplot=np.arange(1,10,1)
plt.errorbar(x,y,yerr=sig_y,fmt='k.')
plt.plot(xplot,par_best[0] + par_best[1]*xplot,'r-',label='Line of best fit')


plt.legend(loc='upper left')
plt.show()
