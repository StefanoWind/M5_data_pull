# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 18:34:08 2024

@author: sletizia
"""
import numpy as np
from matplotlib import pyplot as plt
import warnings
import matplotlib
import os

warnings.filterwarnings('ignore')
plt.close('all')

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'cm' 
matplotlib.rcParams['font.size'] = 14

#%% System
def mkdir(path):
    '''
    Makes recursively folder from path, no existance error
    '''
    import os
    path=path.replace('\\','/')
    folders=path.split('/')
    upper=''
    for f in folders:
        try:
            os.mkdir(upper+f)           
        except:
            pass
                
        upper+=f+'/'

def create_logger(filename,level='info'):
    import logging
    # Create a logger
    logger = logging.getLogger()
    
    if level=='info':
        logger.setLevel(logging.INFO)
    if level=='debug':
        logger.setLevel(logging.DEBUG)
    
    # Clear existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    # Create a file handler
    handler = logging.FileHandler(filename)
    handler.setLevel(logging.DEBUG)

    # Create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)

    # Add the handler to the logger
    if not logger.hasHandlers():
        logger.addHandler(handler)
    
    return logger, handler

def close_logger(logger, handler):
    # Remove handler
    logger.removeHandler(handler)
    # Close handler
    handler.close()


#%% Error analysis
def plot_lin_fit(x,y,lb=[],ub=[],units=''):
    '''
    Plots linear fit and key error metrics
    '''
      
    from scipy.stats import linregress
    
    reals=~np.isnan(x+y)
    
    if np.sum(reals)<=1:
        print("Insufficent number of valid points for linear regression")
        return
    
    if lb==[]:
        lb=np.nanmin(np.append(x,y))
        
    if ub==[]:
        ub=np.nanmax(np.append(x,y))
        
    lf=np.round(linregress(x[reals],y[reals]),3)
    rho=np.round(np.corrcoef(x[reals],y[reals])[0][1],3)
    bias=np.round(np.nanmean(y[reals]-x[reals]),3)
    err_SD=np.round(np.nanstd(y[reals]-x[reals]),3)

    scatter=plt.plot(x,y,'.k',markersize=5,alpha=0.25)
    
    line1_1=plt.plot([lb,ub],[lb,ub],'--g')
    trendline=plt.plot(np.array([lb,ub]),np.array([lb,ub])*lf[0]+lf[1],'r',linewidth=2)
    
    txt=plt.text(lb+(ub-lb)*0.05,lb+(ub-lb)*0.05,r'$y='+str(lf[1])+r'+'+str(lf[0])+r'~x $'+'\n'+
                 r'$\rho='+str(rho)+r'$'+'\n'+r'$\overline{\epsilon}='+str(bias)+'$ '+units+
             '\n'+r'$\sqrt{\overline{\epsilon^{\prime 2}}}='+str(err_SD)+'$ '+units,color='k',fontsize=15,bbox=dict(facecolor=(1,1,1,0.25), edgecolor='k'))

    axis_equal()
    plt.xticks(plt.gca().get_yticks())
    plt.xlim([lb,ub])
    plt.ylim([lb,ub])
    txt_position=txt.get_window_extent().transformed(plt.gca().transData.inverted())
    txt.set_position((lb+(ub-lb)*0.05,lb+(ub-lb)*0.95-txt_position.height))
    plt.grid()
    
    return lf,rho,bias,err_SD,scatter,line1_1,trendline,txt


def simple_bins(x,y,bins=10,bin_method='equally populated'):
    from scipy import stats
    reals=~np.isnan(x+y)
    if isinstance(bins, int):
        if bin_method=='equally populated':
            bins=np.unique(np.nanpercentile(x, np.linspace(0,100,bins)))
        elif bin_method=='equally spaced':
            bins=np.unique(np.linspace(np.nanpercentile(x,0),np.nanpercentile(x,100),bins))
        else:
            raise ValueError('The bin selection method is unkown')
        
    avg=stats.binned_statistic(x[reals], y[reals],statistic='mean',bins=bins)[0]
    std=stats.binned_statistic(x[reals], y[reals],statistic='std',bins=bins)[0]
    
    scatter=plt.plot(x,y,'.k',markersize=5,alpha=0.25)
    errorbar=plt.errorbar(mid(bins), avg, std,color='r',elinewidth=1,barsabove=True,capsize=1)
    scatter_avg=plt.plot(mid(bins), avg,'.r',markersize=5)
    
    return avg,std,bins,scatter, scatter_avg,errorbar

def simple_bins_2d(x,y,z,bins_x=9,bins_y=10,bin_method='equally spaced'):
    from scipy import stats
    reals=~np.isnan(x+y+z)
    if isinstance(bins_x, int):
        if bin_method=='equally populated':
            bins_x=np.unique(np.nanpercentile(x, np.linspace(0,100,bins_x)))
        elif bin_method=='equally spaced':
            bins_x=np.unique(np.linspace(np.nanpercentile(x,0),np.nanpercentile(x,100),bins_x))
        else:
            raise ValueError('The bin selection method is unkown')
    
    if isinstance(bins_y, int):
        if bin_method=='equally populated':
            bins_y=np.unique(np.nanpercentile(y, np.linspace(0,100,bins_y)))
        elif bin_method=='equally spaced':
            bins_y=np.unique(np.linspace(np.nanpercentile(y,0),np.nanpercentile(y,100),bins_y))
        else:
            raise ValueError('The bin selection method is unkown')
    
    avg=stats.binned_statistic_2d(x[reals], y[reals],z[reals],statistic='mean',bins=[bins_x,bins_y])[0]
    pcolor=plt.pcolor(mid(bins_x),mid(bins_y),avg.T,cmap='coolwarm',vmin=np.nanpercentile(avg,5),vmax=np.nanpercentile(avg,95))
    scatter=plt.scatter(x.ravel(),y.ravel(),s=7,c=z.ravel(),cmap='coolwarm',vmin=np.nanpercentile(avg,5),vmax=np.nanpercentile(avg,95),edgecolors='black')
    
    return avg,pcolor,scatter


def bins_unc(x,y,bins,M_BS=100,perc_lim=[5,95],p_value=0.05,bin_method='equally populated'):
    from scipy import stats
    
    if isinstance(bins, int):
        if bin_method=='equally populated':
            bins=np.unique(np.nanpercentile(x, np.linspace(0,100,bins)))
        elif bin_method=='equally spaced':
            bins=np.unique(np.linspace(np.nanpercentile(x,0),np.nanpercentile(x,100),bins))
        else:
            raise ValueError('The bin selection method is unkown')
    
    avg=stats.binned_statistic(x,y,statistic=lambda x:filt_mean(x,perc_lim),bins=bins)[0]
    low=stats.binned_statistic(x,y,statistic=lambda x:filt_BS_mean(x,perc_lim,p_value/2*100,M_BS),bins=bins)[0]
    top=stats.binned_statistic(x,y,statistic=lambda x:filt_BS_mean(x,perc_lim,(1-p_value/2)*100,M_BS),bins=bins)[0]
    
    return avg,low,top

def filt_mean(x,perc_lim=[5,95]):
    '''
    Mean with percentile filter
    '''
    x[x<np.nanpercentile(x,perc_lim[0])]=np.nan
    x[x>np.nanpercentile(x,perc_lim[1])]=np.nan    
    return np.nanmean(x)

def filt_BS_mean(x,perc_lim,p_value,M_BS,min_N=10):
    '''
    Mean with percentile filter and bootstrap
    '''
    x[x<np.nanpercentile(x,perc_lim[0])]=np.nan
    x[x>np.nanpercentile(x,perc_lim[1])]=np.nan
    x=x[~np.isnan(x)]
    
    if len(x)>=min_N:
        x_BS=bootstrap(x,M_BS)
        mean=np.mean(x_BS,axis=1)
        BS=np.nanpercentile(mean,p_value)
    else:
        BS=np.nan
    return BS

def bootstrap(x,M):
    '''
    Bootstrap sample drawer
    '''
    i=np.random.randint(0,len(x),size=(M,len(x)))
    x_BS=x[i]
    return x_BS

    
#%% Dates
def datenum(string,format="%Y-%m-%d %H:%M:%S.%f"):
    '''
    Turns string date into unix timestamp
    '''
    from datetime import datetime
    num=(datetime.strptime(string, format)-datetime(1970, 1, 1)).total_seconds()
    return num

def datestr(num,format="%Y-%m-%d %H:%M:%S.%f"):
    '''
    Turns Unix timestamp into string
    '''
    from datetime import datetime
    string=datetime.utcfromtimestamp(num).strftime(format)
    return string
    

def dt64_to_num(dt64):
    '''
    Converts Unix timestamp into numpy.datetime64
    '''
    tnum=(dt64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    return tnum

def num_to_dt64(tnum):
    '''
    Converts numpy.datetime64 into Unix timestamp
    '''
    dt64= np.datetime64('1970-01-01T00:00:00Z')+np.timedelta64(int(tnum*10**9), 'ns')
    return dt64

#%% Math
def mid(x):
    '''
    Mid points of vector
    '''
    return (x[1:]+x[:-1])/2


def nancorrcoef(x,y=None):
    '''
    Correlation matrix excluding nans.
    '''
    if y is None:
        reals=np.sum(np.isnan(x),axis=0)==0
        return np.corrcoef(x[:,reals])
    else:
        reals=~np.isnan(x+y)
        return np.corrcoef(x[reals],y[reals])
    
def round(x,resolution):
    return np.round(x/resolution)*resolution

def floor(x,resolution):
    return np.floor(x/resolution)*resolution

def ceil(x,resolution):
    return np.ceil(x/resolution)*resolution

def fill_nan(x):
    #fills isolated nan with the mean of the neighbors
    #2023/06/21: created, finalized
    x_int=x.copy()
    for i in range(1,len(x)-1):
        if ~np.isnan(x[i-1]) and ~np.isnan(x[i+1]) and np.isnan(x[i]):
            x_int[i]=(x[i-1]+x[i+1])/2
    return x_int


#%% Algebra
def hstack(a,b):
    if len(np.shape(b))==1:
        b=np.reshape(b,(len(b),1))
    if len(a)>0:
        ab=np.hstack((a,b))
    else:
        ab=b
    return ab


def vstack(a,b):
    '''
    Stack vertically vectors
    '''
    if len(a)>0:
        ab=np.vstack((a,b))
    else:
        ab=b
    return ab   

def len2(x):
    if 'int' in str(type(x)) or 'float' in str(type(x)):
        return 1
    elif 'list' in str(type(x)) or 'array' in str(type(x))  or 'str' in str(type(x)) or 'series' in str(type(x)):
        return len(x)
    else:
        raise ValueError
        
def general_inv(A):
    try:
        A_inv=np.linalg.inv(A)
    except:
        A_inv=np.matmul(np.linalg.inv(np.matmul(np.transpose(A),A)),np.transpose(A))
        
    return A_inv

def match_arrays(x1,x2,max_diff):
    d=np.abs(vec_diff(x1,x2))
    sel=d<max_diff
    matches=np.array([np.where(sel)[0],np.where(sel)[1]]).T
    
    single1=np.setdiff1d(np.arange(len(x1)),matches[:,0])
    single2=np.setdiff1d(np.arange(len(x2)),matches[:,1])
    matches=vstack(matches,np.array([single1,-999*(single1**0)]).T)
    matches=vstack(matches,np.array([-999*(single2**0),single2,]).T).astype(int)
    
    return matches

def vec_diff(x1,x2):
    a=np.tile(x1,(len(x2),1))
    b=np.tile(x2,(len(x1),1))
    return np.transpose(a)-b

#%% Signal processing
def lag(x,y,max_lag=10):
    lag=np.arange(-max_lag,max_lag+1)
    corr=np.zeros(len(lag))
    for l in lag:
        
        if l<0:
            x_lag=x[int(-l):]
            y_lag=y[:int(l)]
        elif l==0:
            x_lag=x.copy()
            y_lag=y.copy()
        else:
            x_lag=x[:-int(l)]
            y_lag=y[int(l):]
        corr[np.where(lag==l)[0][0]] = nancorrcoef(x_lag, y_lag)[0][1]
    
    return lag, corr

#%% Trigonometry
def cosd(x):
    return np.cos(x/180*np.pi)

def sind(x):
    return np.sin(x/180*np.pi)
    
def tand(x):
    return np.tan(x/180*np.pi)

def arctand(x):
    return np.arctan(x)*180/np.pi

def arccosd(x):
    return np.arccos(x)*180/np.pi

def arcsind(x):
    return np.arcsin(x)*180/np.pi

def ang_diff(angle1,angle2=None,unit='deg',mode='scalar'):
    try:
        if 'list' in str(type(angle1)) or 'int' in str(type(angle1)):
            angle1=np.array(angle1).astype(float)
        if 'list' in str(type(angle1)) or 'int' in str(type(angle2)):
            angle2=np.array(angle2).astype(float)
    except: pass
    if unit=='rad':
        angle1=angle1*180/np.pi
        angle2=angle2*180/np.pi
    angle1=angle1 % 360
    
    if mode=='scalar':
        angle2=angle2 % 360
    
   
    if mode=='scalar':
        dx=angle1-angle2
    elif mode=='vector':
        dx=np.diff(angle1)

    if len2(dx)>1:
        dx[dx>180]= dx[dx>180]-360
        dx[dx<-180]= dx[dx<-180]+360
    else:
        if dx>180:
            dx-=360
        if dx<-180:
            dx+=360
        
    if unit=='rad':
        dx=dx/180*np.pi
        
    return dx

def cart2pol(x,y):
    r=(x**2+y**2)**0.5
    th=arctand(y/x)
    if len2(th)>1:
        th[x<0]-=180
    else:
        if x<0:
            th-=180
    return r,th%360

#%% Graphics
def axis_equal():
    '''
    Makes axis of plot equal
    '''
    from mpl_toolkits.mplot3d import Axes3D
    ax=plt.gca()
    is_3d = isinstance(ax, Axes3D)
    if is_3d:
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        zlim=ax.get_zlim()
        ax.set_box_aspect((np.diff(xlim)[0],np.diff(ylim)[0],np.diff(zlim)[0]))
    else:
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        ax.set_box_aspect(np.diff(ylim)/np.diff(xlim))
        
        
def remove_labels(fig):
    '''
    Removes duplicated labels from multiplot
    '''
    axs=fig.axes

    for ax in axs:
        loc=ax.get_subplotspec()
        try:
            
            if loc.is_last_row()==False:
                ax.set_xticks(ax.get_xticks(),[])
                ax.set_xlabel('')
            if loc.is_first_col()==False:
                ax.set_yticks(ax.get_yticks(),[])
                ax.set_ylabel('')
        except:
            pass

def save_all_fig(name,cd,newfolder=False,resolution=300):
    '''
    Saves all current figures
    '''
    mkdir(os.path.join(cd,'figures'))
    if newfolder:
        mkdir(os.path.join(cd,'figures',name))
    figs = [plt.figure(n) for n in plt.get_fignums()]
    inc=0
    for fig in figs:
        if newfolder:
            fig.savefig(os.path.join(cd,'figures',name,'{i:02d}'.format(i=inc)+'.png'),dpi=resolution, bbox_inches='tight')
        else:
            fig.savefig(os.path.join(cd,'figures',name+'{i:02d}'.format(i=inc)+'.png'),dpi=resolution, bbox_inches='tight')
        inc+=1
        
def draw_turbine(x,y,D,wd):
    import matplotlib.image as mpimg
    from matplotlib import transforms
    from matplotlib import pyplot as plt
    img = mpimg.imread('C:/Users/SLETIZIA/OneDrive - NREL/Desktop/PostDoc/utils/Turbine5.png')
    ax=plt.gca()
    xlim=ax.get_xlim()
    ylim=ax.get_ylim()
    tr = transforms.Affine2D().scale(D/700,D/700).translate(-100*D/700,-370*D/700).rotate_deg(90-wd).translate(x,y)
    ax.imshow(img, transform=tr + ax.transData)
    plt.xlim(xlim)
    plt.ylim(ylim)
    
def add_gap(x,y,F,gap_x,gap_y,ext_x,ext_y):
    
    diff_x=np.append(np.diff(x)[0],np.diff(x))
    i_gap_x=np.concatenate([[0],np.where(diff_x>gap_x)[0]])
    
    if isinstance(x[0], np.datetime64):
        x2=np.array([],dtype= 'datetime64')
    else:    
        x2=[]
    F2=[]
    
    if ext_x!=0:
        x2=np.append(x2,x[0]-ext_x)
        F2=hstack(F2,np.reshape(F[:,0],(-1,1)))
    
    for i1,i2 in zip(i_gap_x[:-1],i_gap_x[1:]):
        x2=np.append(x2,x[i1:i2])
        F2=hstack(F2,F[:,i1:i2])
        
        if ext_x!=0:
            x2=np.append(x2,x[i2-1]+ext_x)
            F2=hstack(F2,np.reshape(F[:,i2-1],(-1,1)))
        
        x2=np.append(x2,x[i2-1]+diff_x[i2]/2)
        F2=hstack(F2,np.zeros((len(y),1))+np.nan)
        
        if ext_x!=0:
            x2=np.append(x2,x[i2]-ext_x)
            F2=hstack(F2,np.reshape(F[:,i2],(-1,1)))
            
    x2=np.append(x2,x[i_gap_x[-1]:])
    F2=hstack(F2,F[:,i_gap_x[-1]:])
    
    diff_y=np.append(np.diff(y)[0],np.diff(y))
    i_gap_y=np.concatenate([[0],np.where(diff_y>gap_y)[0]])
    
    if isinstance(y[0], np.datetime64):
        y2=np.array([],dtype= 'datetime64')
    else:    
        y2=[]
    F3=[]
    
    if ext_y!=0:
        y2=np.append(y2,y[0]-ext_y)
        F3=vstack(F3,np.reshape(F2[0,:],(1,-1)))
    for i1,i2 in zip(i_gap_y[:-1],i_gap_y[1:]):
        y2=np.append(y2,y[i1:i2])
        F3=vstack(F3,F2[i1:i2,:])
        
        if ext_y!=0:
            y2=np.append(y2,y[i2-1]+ext_y)
            F3=vstack(F3,np.reshape(F2[i2-1,:],(1,-1)))
        
        y2=np.append(y2,y[i2-1]+diff_y[i2]/2)
        F3=vstack(F3,np.zeros((1,len(x2)))+np.nan)
        
        if ext_y!=0:
            y2=np.append(y2,y[i2]-ext_y)
            F3=vstack(F3,np.reshape(F2[i2,:],(1,-1)))
            
    y2=np.append(y2,y[i_gap_y[-1]:])
    F3=vstack(F3,F2[i_gap_y[-1]:,:])
    
    return x2,y2,F3

def make_video(folder,output,fps=1):
    import cv2
    import os
    
    # Get the list of PNG images in the directory
    image_files = [file for file in os.listdir(folder) if file.endswith('.png')]
    
    print(str(len(image_files))+' images found')
    # Sort the image files to ensure proper ordering
    image_files.sort()
    
    # Get the first image to determine the video dimensions
    first_image = cv2.imread(os.path.join(folder, image_files[0]))
    height, width, _ = first_image.shape
    
    # Define the video codec and create a VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Use appropriate codec for your desired video format
    video_writer = cv2.VideoWriter(output, fourcc, fps, (width, height))  # Adjust the frame rate (here set to 25)
    
    # Iterate over the image files and write each frame to the video
    for image_file in image_files:
        image_path = os.path.join(folder, image_file)
        frame = cv2.imread(image_path)
        
        # Resize frame if its dimensions are different from the first image
        if frame.shape[0] != height or frame.shape[1] != width:
            frame = cv2.resize(frame, (width, height))
        
        video_writer.write(frame)
    
    # Release the VideoWriter and close the video file
    video_writer.release()
    
    print('Video saved as '+output)
    

#%% Machine learning
def RF_feature_selector(X,y,test_size=0.8,n_search=30,n_repeats=10,limits={}):
    '''
    Feature importance selector based on random forest. The optimal set of hyperparameters is optimized through a random search.
    Importance is evaluated throuhg the permutation method, which gives higher scores to fatures whose error metrics drops more after reshuffling.
    '''
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.inspection import permutation_importance
    from scipy.stats import randint
    from sklearn.model_selection import train_test_split
    from sklearn.model_selection import RandomizedSearchCV
    from sklearn.metrics import mean_absolute_error

    #build train/test datasets
    data = np.hstack((X, y.reshape(-1, 1)))

    data = data[~np.isnan(data).any(axis=1)]
    train_set, test_set = train_test_split(data, random_state=42, test_size=test_size)

    X_train = train_set[:,0:-1]
    y_train = train_set[:,-1]

    X_test = test_set[:,0:-1]
    y_test = test_set[:,-1]
    
    #default grid of hyperparamters (Bodini and Optis, 2020)
    if limits=={}:
        p_grid = {'n_estimators': randint(low=10, high=100), # number of trees
                  'max_features': randint(low=1,high= 6), # number of features to consider when looking for the best split
                  'min_samples_split' : randint(low=2, high=11),
                  'max_depth' : randint(low=4, high=10),
                  'min_samples_leaf' : randint(low=1, high=15)
            }
        
    # Choose cross-validation techniques for the inner and outer loops,
    # independently of the dataset.
    forest_reg = RandomForestRegressor()
    rnd_search = RandomizedSearchCV(forest_reg, param_distributions = p_grid, n_jobs = -1,
                                    n_iter=n_search, cv=5, scoring='neg_mean_squared_error')
    rnd_search.fit(X_train, y_train)
    print('Best set of hyperparameters found:')
    print(rnd_search.best_estimator_)

    predicted_test = rnd_search.best_estimator_.predict(X_test)
    test_mae = mean_absolute_error(y_test, predicted_test)
    print("Average testing MAE:", test_mae)

    predicted_train = rnd_search.best_estimator_.predict(X_train)
    train_mae = mean_absolute_error(y_train, predicted_train)
    print("Average training MAE:", train_mae)

    best_params=rnd_search.best_estimator_.get_params()    
        
    #random forest prediction with optimized hyperparameters
    reals=np.sum(np.isnan(np.hstack((X, y.reshape(-1, 1)))),axis=1)==0
    rnd_search.best_estimator_.fit(X[reals,:], y[reals])
        
    y_pred=y+np.nan
    y_pred[reals] = rnd_search.best_estimator_.predict(X[reals])
       
    reals=~np.isnan(y+y_pred)
    result = permutation_importance(rnd_search.best_estimator_, X[reals], y[reals], n_repeats=n_repeats, random_state=42, n_jobs=2)

    importance=result.importances_mean
    importance_std=result.importances_std
    
    return importance,importance_std,y_pred,test_mae,train_mae,best_params
