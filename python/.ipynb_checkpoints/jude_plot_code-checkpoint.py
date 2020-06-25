import pandas as pd
import numpy as np
import math as math
import random as random
import sys
import copy as copy
from scipy.stats import norm
import statsmodels.api as sm
from scipy import stats
from numpy.linalg import inv
import pymc3
import scipy.stats as sst
from sklearn import preprocessing
#import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams.update({'errorbar.capsize': 2})

def do_BiasBoxplot(results):
    fig= plt.figure(1, figsize=(9, 6))
    # Create an axes instance
    ax = fig.add_subplot(111)
    # Create the boxplot
    bp = ax.boxplot(results, showfliers=0)
    plt.rc('ytick',labelsize=18)
    ## add patch_artist=True option to ax.boxplot()
    ## to get fill color
    bp = ax.boxplot(results, patch_artist=True, showfliers=0)
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = '#1b9e77' )
    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)
    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)
    
    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)
    
    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)
    ax.set_xticklabels([r'$g_J$', r'$g_Y$', r'$T_{opt}$', r'$s$', r'$R_0$', r'$d_J$', r'$m_J$', r'$m_Y$', r'$m_A$',  r'$K$'], fontsize=18)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(width = 2, direction = "out")
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    ax.set_ylabel('Bias in Parameter Estimation', fontsize=18)
    ax.set_xlabel('Parameters',fontsize=18)
    # Save the figure
    return fig, ax
def do_forcast( N1, temperatures, time, filename):
    yerr=np.empty((2, len(temperatures)))
    N1L=np.percentile(N1, 25, axis =0)
    N1M=np.percentile(N1, 50, axis =0)
    N1U=np.percentile(N1, 75, axis =0)
    yerr[0,:]=N1M-N1L
    yerr[1,:]=N1U-N1M
    fig=plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.errorbar(time, N1M, yerr, fmt='k--', elinewidth=2)
    ax.set_ylabel('Bias in abundance', fontsize=18)
    ax.set_xlabel('Time (years)',fontsize=18)
    fig.savefig(filename)
    plt.close()
def do_realdata(NS, NO, filename):
    NO=pd.Series(NO,index=pd.Series(range(1980,2013)))
    NS=pd.Series(NS,index=pd.Series(range(1980,2013)))
    fig, ax=plt.subplots()
    plt.plot(NO.index, NO, 'k', linewidth=2, label='Observed')
    plt.plot(NS.index, NS, 'r', linewidth=2, label='Simulated')
    plt.legend(loc='best')
    ax.set_ylabel('Abundance', fontsize=18)
    ax.set_xlabel('Time (years)',fontsize=18)
    ax.tick_params(width = 2, direction = "out")
    ax.set_xticks([1980, 1985, 1990, 1995, 2000, 2005, 2013])
    ax.set_xticklabels([1980, 1985, 1990, 1995, 2000, 2005, 2013])
    fig.savefig(filename, bbox_inches='tight')
    plt.close()
    
def do_scatter(NS, NO, filename):
    NO=pd.Series(NO,index=pd.Series(range(1968, 2001)))
    NS=pd.Series(NS,index=pd.Series(range(1968,2001)))
    fig, ax=plt.subplots()
    plt.scatter(NO.index, NO, s=700, c='r', alpha=0.5, label='Observed')
    plt.scatter(NS.index, NS,  s=700, c='c', alpha=0.5, label='Simulated')
    plt.legend(loc='best')
    ax.set_ylabel('Abundance', fontsize=18)
    ax.set_xlabel('Time (years)',fontsize=18)
    ax.tick_params(width = 2, direction = "out")
    ax.set_xticks([1968, 1978, 1988, 1998])
    ax.set_xticklabels([1968, 1978, 1988, 1998])
    fig.savefig(filename, bbox_inches='tight')
    plt.close()

