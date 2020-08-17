# this script simulates stage-structured population sizes in space and time for given parameters and forcin
# functions
from __future__ import print_function, division
import sys

# NB for EAM: this runs on Mushu in pipenv, so call pipenv run python fish_abc.py

# import dependencies
import numpy as np
import random as random
import math as math
from scipy.stats import norm
from sklearn.metrics import mean_squared_error
from random import triangular
import scipy.stats as sst
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
#price = pd.Series(np.random.randn(150).cumsum(),index=pd.date_range('2000-1-1', periods=150, freq='B'))
########################################
##################################################################
#There 16 lat( from 30 to 45), and 33 years (from 1968 -2000). length  at maturity is 68.3
df = pd.read_csv("jude_groundfish_training.csv", usecols=['sppocean','Freq', 'year', 'wtcpue', 'lat', 'bottemp', 'surveyfact', 'LENGTH'])
df=df.loc[(df['sppocean'] =='gadus morhua_Atl')]
df=df.loc[(df['surveyfact'] =='NEFSC_NEUSSpring')]
df.lat = df.lat.round().astype(np.int)
df=df.interpolate(method ='linear', limit_direction ='forward')
df=df.interpolate(method ='linear', limit_direction ='backward')
#print(np.where(pd.isnull(df.SURFTEMP)))
#Temp = pd.isnull(df["SURFTEMP"])
#Pop = pd.isnull(df["ABUNDANCE"])
# filling a missing value with
# previous ones
#df.interpolate(method ='linear', limit_direction ='forward')
#df.interpolate(method ='linear', limit_direction ='backward')

#df.fillna(method ='pad')
# filling  null value using fillna() function
#df.fillna(method ='bfill')
#print(df[Temp])
#print(df[Pop])
#np.where(df.SURFTEMP.applymap(lambda x: x == ''))
#pd.isna(df.SURFTEMP['column_name'])
#sys.exit()
#df1=df.groupby(["LAT", "SURFTEMP", "YEAR", "LENGTH" ])["ABUNDANCE"].sum()
##############################################################################
#print(df[0:5])
#print(df['ABUNDANCE'].max())
#sys.exit()
#sys.exit()
nn=df['lat'].max()-df['lat'].min()#15
#print(df['LENGTH'].max())
#print(df['LENGTH'].min())
#sys.exit()
#print(df['lat'].min())30
#print(df['lat'].max())45
#print(df['year'].min())
#print(df['year'].max())
#print(nn)
#sys.exit()
#print(df['YEAR'].min())
#sys.exit()
D={}
D1={}
#D1={}
#for q in range(1,2):
#   D['J_patch'+ str(q)]=df.loc[(df['LAT'] == 28+ q) & (df['LENGTH']<13)]
#  D['Y_patch'+ str(q)]=df[(df['LAT'] == 28+ q) & (df['LENGTH'].isin(range(13,27)))]
# D['A_patch'+ str(q)]=df[(df['LAT'] == 28+ q) & (df['LENGTH']>26)]
#print(D['J_patch'+ str(q)].ABUNDANCE.values)
#print(len(D['J_patch'+ str(q)].ABUNDANCE.values))
#print(len(D['J_patch'+ str(q)].YEAR.values))
#sys.exit()(1,18):
for q in range(1,nn+2):
    D['J_patch'+ str(q)]=df.loc[(df['lat'] == 29+ q) & (df['LENGTH']<34)]
    n=len(D['J_patch'+ str(q)].year.values)
    m=df['year'].max()-df['year'].min()
    A=np.empty((n, 3))
    Abun_TemJ=np.empty((m+1, 3))
    A[:,0]=D['J_patch'+ str(q)].year.values
    A[:,1]=D['J_patch'+ str(q)].bottemp.values
    A[:,2]=D['J_patch'+ str(q)].Freq.values
    #####################################
    kJ=0
    kY=0
    kA=0
    for i in range (0, m+1):
        Abun_TemJ[i,0]=kJ
        N=np.zeros(len(A[:,0]))
        T=np.zeros(len(A[:,0]))
        for j in range (0, len(A[:,0])):
            if A[j,0] == 1968+i:
               T[j]=A[j,1]
               N[j]=A[j,2]
            else:
               T[j]=0
               N[j]=0
        Abun_TemJ[i,1]=np.nan_to_num(T[T!=0].mean())
        Abun_TemJ[i,2]=np.sum(N)
        kJ=kJ +1
    D1['J_patch'+ str(q)]=Abun_TemJ
#print( D1['J_patch'+ str(q)][:,2].max())
 ##############################
    D['Y_patch'+ str(q)]=df[(df['lat'] == 29+ q) & (df['LENGTH'].isin(range(34,68)))]
    n=len(D['Y_patch'+ str(q)].year.values)
    m=df['year'].max()-df['year'].min()
    B=np.empty((n, 3))
    Abun_TemY=np.empty((m+1, 3))
    B[:,0]=D['Y_patch'+ str(q)].year.values
    B[:,1]=D['Y_patch'+ str(q)].bottemp.values
    B[:,2]=D['Y_patch'+ str(q)].Freq.values
    for i in range (0, m+1):
        Abun_TemY[i,0]=kY
        N=np.zeros(len(B[:,0]))
        T=np.zeros(len(B[:,0]))
        for j in range (0, len(B[:,0])):
            if B[j,0] == 1968+i:
                T[j]=B[j,1]
                N[j]=B[j,2]
            else:
                T[j]=0
                N[j]=0
        Abun_TemY[i,1]=np.nan_to_num(T[T!=0].mean())
        Abun_TemY[i,2]=np.sum(N)
        kY=kY +1
    D1['Y_patch'+ str(q)]=Abun_TemY
# print( D1['Y_patch'+ str(q)][:,2].max())
 #############################
    D['A_patch'+ str(q)]=df[(df['lat'] == 29+ q) & (df['LENGTH']>67)]
    n=len(D['A_patch'+ str(q)].year.values)
    m=df['year'].max()-df['year'].min()
    C=np.empty((n, 3))
    Abun_TemA=np.empty((m+1, 3))
    C[:,0]=D['A_patch'+ str(q)].year.values
    C[:,1]=D['A_patch'+ str(q)].bottemp.values
    C[:,2]=D['A_patch'+ str(q)].Freq.values
    for i in range (0, m+1):
        Abun_TemA[i,0]=kA
        N=np.zeros(len(C[:,0]))
        T=np.zeros(len(C[:,0]))
        for j in range (0, len(C[:,0])):
            if C[j,0] == 1968+i:
                T[j]=C[j,1]
                N[j]=C[j,2]
            else:
                T[j]=0
                N[j]=0
        Abun_TemA[i,1]=np.nan_to_num(T[T!=0].mean())
#Abun_TemA[i,1]=np.amax(T)
        Abun_TemA[i,2]=np.sum(N)
        kA=kA +1
    D1['A_patch'+ str(q)]=Abun_TemA
# print( D1['J_patch'+ str(q)][:,2].max())
 ##############################

#A=np.array([[1963, 1963, 1963, 1968],[1,2,3,4], [7,6,5,4]])
#A=np.transpose(A)
#B=np.empty((10, 3))
#k=0
#for i in range (0, 10):
#    B[i,0]=k
#    C=np.zeros(len(A[:,0]))
#  D=np.zeros(len(A[:,0]))
#   for j in range (0, len(A[:,0])):
#       if A[j,0] == 1963+i:
#         C[j]=A[j,2]
#          D[j]=A[j,1]
#      else:
#           C[j]=0
#           D[j]=0
#   B[i,1]=np.nan_to_num(D[D!=0].mean())
#   B[i,2]=np.sum(C)#  k=k+1
#print(B)
#D.groupby(["LAT", "YEAR", "SURFTEMP"])["ABUNDANCE"].sum()
#df.groupby(["id", "year"])["value"].sum()
#D1=D['J_patch'+ str(q)].loc[D['J_patch'+ str(q)]['YEAR'] == i, ['ABUNDANCE', 'BIOMASS']].sum()
#for i in range(1963,2019):
#df.loc[df['YEAR'] == 1, '['A', 'D']'].sum()
#print(df1[0:5])
#D['A_bias'+ str(q)]=df[(df['LAT'] == 29) & df['LENGTH']>26]
###############################################################################
#print(df1[0:5])
#print(D['J_patch1'].SURFTEMP[:5].values)
#print(D['J_patch1'].YEAR.min())
#print(D['A_patch2'][:5])
#print(D['Y_patch3'][:5])
#rslt_df = dataframe[dataframe['Percentage'] > 80]
#rslt_df = dataframe.loc[dataframe['Percentage'] > 80]
#print('\nResult dataframe :\n', rslt_df)

#options = ['Math', 'Commerce']

# selecting rows based on condition
#rslt_df = dataframe[dataframe['Stream'].isin(options)]

#print('\nResult dataframe :\n', rslt_df)

# selecting rows based on condition
#rslt_df = dataframe.loc[~dataframe['Stream'].isin(options)]
# selecting rows based on condition
#rslt_df = dataframe[(dataframe['Age'] == 21) &
#          dataframe['Stream'].isin(options)]
# selecting rows based on condition
#rslt_df = dataframe.loc[(dataframe['Age'] == 21) &
#                       dataframe['Stream'].isin(options)]

#df.loc[df['a'] == 1, 'b'].sum()
