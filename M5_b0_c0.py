# -*- coding: utf-8 -*-
import os
cd=os.path.dirname(__file__)
    
import pandas as pd
import numpy as np
import utils as utl
from matplotlib import pyplot as plt
import warnings
import glob
import matplotlib
import matplotlib.dates as mdates

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.size'] = 12

warnings.filterwarnings('ignore')
plt.close('all')

#%% Inputs
source=os.path.join(cd,'data/nwtc/nwtc.m5.b0/*csv')# make sure this are the QC data
source_tilt=os.path.join(cd,'data/20220512.000000-20220809.000000_m5_tilt.csv')
z_S=np.array([41,61,74,119])#[m] sonics' heights
WD_offset=8#[deg] offset with respect to true north, it will be added

#%% Initialization
files=glob.glob(source)
Output=pd.DataFrame()

TILT=pd.read_csv(source_tilt)
TILT=TILT.set_index('z [m AGL]')

#%% Main
for f in files:
    Data=pd.read_csv(f)
    Data=Data.set_index('Time (UTC)')
    
    #tilt correction
    for z in z_S:
        if 'Sonic_x_clean_{z}m [m/s]'.format(z=z) in Data.keys():
        
            Um=Data['Sonic_x_clean_{z}m [m/s]'.format(z=z)]
            Vm=Data['Sonic_y_clean_{z}m [m/s]'.format(z=z)]
            Wm=Data['Sonic_z_clean_{z}m [m/s]'.format(z=z)]

            pitch=TILT['Pitch [deg]'][z]
            roll=TILT['Roll [deg]'][z]
            bias=TILT['Bias [m/s]'][z]
            
            if np.isnan(pitch):
                pitch=0
            if np.isnan(roll):
                roll=0
            if np.isnan(bias):
                bias=0
                
            C=np.array([[1,0,0],[0,utl.cosd(roll),-utl.sind(roll)],[0,utl.sind(roll),utl.cosd(roll)]])
            D=np.array([[utl.cosd(pitch),0,utl.sind(pitch)],[0,1,0],[-utl.sind(pitch),0,utl.cosd(pitch)]])
           
            P=np.matmul(D.T,C.T)
            
            UVWm=np.array([Um,Vm,Wm-bias])
            UVWp=np.matmul(P,UVWm)
            
            Data['Sonic_x_corr_{z}m [m/s]'.format(z=z)]=UVWp[0,:]
            Data['Sonic_y_corr_{z}m [m/s]'.format(z=z)]=UVWp[1,:]
            Data['Sonic_z_corr_{z}m [m/s]'.format(z=z)]=UVWp[2,:]
            
    Data_avg=Data.mean()
    
    #mean wind
    for z in z_S:
        if 'Sonic_x_clean_{z}m [m/s]'.format(z=z) in Data.keys():
            
            #no tilt correctin
            WS,WD0=utl.cart2pol(Data_avg['Sonic_x_clean_{z}m [m/s]'.format(z=z)],Data_avg['Sonic_y_clean_{z}m [m/s]'.format(z=z)])
            WD=(270-WD0)%360
            Data_avg['Sonic_WS_{z}m [m/s]'.format(z=z)]=WS
            Data_avg['Sonic_WD_{z}m [deg]'.format(z=z)]=WD+WD_offset
            
            #tilt correction
            WS,WD0=utl.cart2pol(Data_avg['Sonic_x_corr_{z}m [m/s]'.format(z=z)],Data_avg['Sonic_y_corr_{z}m [m/s]'.format(z=z)])
            WD=(270-WD0)%360
            Data_avg['Sonic_WS_corr_{z}m [m/s]'.format(z=z)]=WS
            Data_avg['Sonic_WD_corr_{z}m [deg]'.format(z=z)]=WD+WD_offset
        
    Output=Output.append(Data_avg,ignore_index=True)

#%% Output
Output['Time (UTC)']=[utl.num_to_dt64(t) for t in Output['Timestamp']]
Output=Output.set_index('Time (UTC)')

name_save=os.path.join(cd,'data/'+utl.datestr(Output['Timestamp'].min(),'%Y%m%d.%H%M%S')+'-'+\
                      utl.datestr(Output['Timestamp'].max(),'%Y%m%d.%H%M%S')+'.'+os.path.basename(f)[:-20].replace('b0','c0')+'.csv')
Output.to_csv(name_save)
    
print('10-min statistics saved as ' + name_save)
    
#%% Plots
plt.figure(figsize=(18,10))
ctr=1
for z in z_S:
    ax=plt.subplot(len(z_S),1,ctr)
    plt.plot(Output.index,Output['Sonic_WS_{z}m [m/s]'.format(z=z)],'-r',label='No tilt correction')
    plt.plot(Output.index,Output['Sonic_WS_corr_{z}m [m/s]'.format(z=z)],'-k',label='Tilt corrected')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
    plt.title(f'$z$={z} m')
    plt.grid()
    plt.ylabel('Wind speed [m s$^{-1}$]')
    if ctr==0:
        plt.legend()
    ctr+=1
plt.xlabel('Time (UTC)')

plt.tight_layout()
plt.savefig(name_save.replace('.csv','.wind_speed.png'))

plt.figure(figsize=(18,10))
ctr=1
for z in z_S:
    ax=plt.subplot(len(z_S),1,ctr)
    plt.plot(Output.index,Output['Sonic_WD_{z}m [deg]'.format(z=z)],'-r',label='No tilt correction')
    plt.plot(Output.index,Output['Sonic_WD_corr_{z}m [deg]'.format(z=z)],'-k',label='Tilt corrected')
    plt.ylim([0,360])
    plt.title(f'$z$={z} m')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
    plt.grid()
    plt.ylabel('Wind direction [$^\circ$]')
    if ctr==0:
        plt.legend()
    ctr+=1
plt.xlabel('Time (UTC)')
plt.tight_layout()
plt.savefig(name_save.replace('.csv','.wind_direction.png'))
