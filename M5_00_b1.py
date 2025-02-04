# -*- coding: utf-8 -*-
'''
Pull 10-min data from M5
'''
import os
cd=os.path.dirname(__file__)
import utils as utl
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
import glob
import datetime as dt
import pandas as pd
import warnings
plt.close('all')
warnings.filterwarnings("ignore")

#%% Inputs
source='Y:/Wind-data/Public/Projects/Met135/MetData/M5Twr'
storage=os.path.join(cd,'data/nwtc/nwtc.m5.b1/')
date1='2023-04-01'#start date
date2='2023-06-30'#end date
zero_datenum=719529#[days] 1970-01-01 in matlab time

z_sonic=[15,41,61,74,100,119]#[m] sonic heights
z_cup=[3,10,30,38,55,80,87,105,122,130]#[m] cup-vane height
vars_sonic=['Wind_Direction_Sonic_{z}m','Wind_Speed_CupEq_Sonic_{z}m','Ti_CupEq_Sonic_{z}m']
vars_cup=['Wind_Direction_Vane_{z}m_mean', 'Wind_Direction_Vane_{z}m_sdev','Wind_Speed_Cup_{z}m','Ti_Cup_{z}m']

#%% Functions
def _check_keys( dict):
    """
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    """
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


def loadmat(filename):
    """
    this function should be called instead of direct scipy.io .loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

#%% Initialization
os.makedirs(os.path.join(storage),exist_ok=True)

d1=dt.datetime.utcfromtimestamp(utl.datenum(date1,'%Y-%m-%d'))
d2=dt.datetime.utcfromtimestamp(utl.datenum(date2,'%Y-%m-%d'))

days=[] 
d=d1
while d<=d2:
    days.append(d)
    d+=dt.timedelta(1,0,0)
    
#%% Main
for d in days:      
    files=glob.glob(os.path.join(source,d.strftime('%Y/%m/%d'),'summary_data','*.mat'))
    if len(files)>0:
        Data=pd.DataFrame()   

        for f in files:
            D=pd.DataFrame([])
            
            mat = loadmat(f)
            D.index=[(mat['Raw_Sonic_x_15_mean']['date']-zero_datenum)*24*3600]
            for v in vars_sonic:
                for z in z_sonic:
                    x=mat[v.format(z=z)]['val']
                    u=mat[v.format(z=z)]['units']
                    D[v.format(z=z)+' ('+u+')']=x
            for v in vars_cup:
                for z in z_cup:
                    x=mat[v.format(z=z)]['val']
                    u=mat[v.format(z=z)]['units']
                    D[v.format(z=z)+' ('+u+')']=x
                
            Data=pd.concat([Data,D])
            print('\r'+f+' done')
            
        # Output
        Data.to_csv(os.path.join(storage,'nwtc.m5.b1.'+dt.datetime.strftime(d,'%Y%m%d.%H%M%S')+'.csv'))
    
    else:
        print(f'No files on {d.strftime("%Y%m%d")}')
print('Complete')
input()
