# -*- coding: utf-8 -*-

import os
cd=os.path.dirname(__file__)
import utils as utl
import scipy.io as spio
import matplotlib.pyplot as plt
import glob
import datetime as dt
import pandas as pd
import time
import os
import re
import warnings
plt.close('all')
warnings.filterwarnings("ignore")

#%% Inputs
source='Y:/Wind-data/Public/Projects/Met135/MetData/M5Twr'
storage=os.path.join(cd,'data/nwtc/nwtc.m5.a0/')
sdate='2024-01-01'
edate='2024-01-02'
replace=False

zero_datenum=719529#[days] 1970-01-01 in Inputlab time
variables=['Sonic_x_clean_(\d+)m','Sonic_y_clean_(\d+)m','Sonic_z_clean_(\d+)m']

#%% Initialization
utl.mkdir(storage)
d1=dt.datetime.utcfromtimestamp(utl.datenum(sdate,'%Y-%m-%d'))
d2=dt.datetime.utcfromtimestamp(utl.datenum(edate,'%Y-%m-%d'))

days=[]
d=d1
while d<=d2:
    days.append(d)
    d+=dt.timedelta(1,0,0)

for d in days:      
    files=glob.glob(os.path.join(source,d.strftime('%Y/%m/%d'),'raw_data','*.mat'))
        
    counter=0
    t0=time.time()
        
    #%% Main
    for f in files:
        try:
            Output=pd.DataFrame()
            Input = spio.loadmat(f, struct_as_record=False, squeeze_me=True)
            tnum=(Input['time_UTC'].val-zero_datenum)*24*3600
            filename='nwtc.m5.a0.'+utl.datestr(tnum[0],'%Y%m%d.%H%M%S')+'.csv'

            if replace==False and os.path.exists(os.path.join(storage,filename))==True:
                print(filename+' skipped')
            else:
                for v in variables:
                    heights= [Inputch.group(1) for string in Input.keys() if (Inputch := re.search(v, string))]
                    for h in heights:
                        varname=v.replace('(\d+)',h)
                        Output[varname+' ['+Input[varname].units+']']=Input[varname].val
                
                Output['Timestamp']=tnum
                Output['Time (UTC)']=[utl.num_to_dt64(t) for t in tnum]
                
                Output=Output.set_index('Time (UTC)')
                Output.to_csv(os.path.join(cd,storage,filename))
                counter+=1
                print(f+' done')
        
        except Exception as e:
            print('Error at '+f+': '+str(e))
            pass