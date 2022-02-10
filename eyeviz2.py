import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
import ipywidgets as widgets
from ipywidgets import interact
import pandas as pd
from IPython.display import display, clear_output
from matplotlib.colors import LinearSegmentedColormap
cmap1 = LinearSegmentedColormap.from_list('rg',["r", "w", "g"], N=256)
from mpl_toolkits.axes_grid1 import make_axes_locatable

patientdf = pd.read_csv('TD_vf_vals.csv')
displacements = pd.read_csv('VF_Locations.csv')

displacements['xpos']=displacements['x_morphed']
displacements['ypos']=displacements['y_morphed']

#x = %system ls 'data/'
from os import listdir
from os.path import isfile, join
x = [f for f in listdir('data') if isfile(join('data',f))]
x.sort()

patientnums = [i.split('-')[0] for i in x]
patientthicks = [i.split('-')[1].replace('.hdf5','') for i in x]
h5files = ['data/'+i for i in x]

numsunique = []
for i in patientnums:
    if i not in numsunique:
        numsunique.append(i)
thickunique = []
for i in patientthicks:
    if i not in thickunique:
        thickunique.append(i)
		
h5files.sort()
h5f = h5py.File(h5files[0],'r')
for i in h5f['070/L'].keys():
    yunique = [int(j) for j in h5f['070/L/'+str(i)+''].keys()]
    yunique.sort(reverse=True)
    xunique = yunique.copy()
    break
h5f.close()

angles = pd.read_csv('angles.csv')

def closest(lst, K):
    x = []
    for i in range(len(lst)):
        x.append(abs(lst.iloc[i]-K))
    m = x.index(min(x))
    return m

patientdf.loc[:,'VISIT_DATE'] = pd.to_datetime(patientdf.loc[:,'VISIT_DATE'],format='%m/%d/%y')

def checkeyeside(patient_num,eyeside):
    clear_output(wait=True)
    if eyeside == 'R':
        if os.path.exists('infrared/'+patient_num+'_OD_SPEC.JPG'):
            im = plt.imread('infrared/'+patient_num+'_OD_SPEC.JPG')
        elif os.path.exists('infrared/'+patient_num+'_OD_SPEC.jpg'):
            im = plt.imread('infrared/'+patient_num+'_OD_SPEC.jpg')
        elif os.path.exists('infrared/'+patient_num+'_OD_SPEC.jpg'):
            im = plt.imread('infrared/'+patient_num+'_OD_SPEC.jpeg')
        else:
            clear_output(wait=True)
            #plt.figure().clf()
            print('No '+patient_num+' R infrared image exists')
            return('No '+patient_num+' R infrared image exists')
    else:
        if os.path.exists('infrared/'+patient_num+'_OS_SPEC.JPG'):
            im = plt.imread('infrared/'+patient_num+'_OS_SPEC.JPG')
        elif os.path.exists('infrared/'+patient_num+'_OS_SPEC.jpg'):
            im = plt.imread('infrared/'+patient_num+'_OS_SPEC.jpg')
        elif os.path.exists('infrared/'+patient_num+'_OS_SPEC.jpeg'):
            im = plt.imread('infrared/'+patient_num+'_OS_SPEC.jpeg')
        else:
            clear_output(wait=True)
            #plt.figure().clf()
            print('No '+patient_num+' L infrared image exists')
            return('No '+patient_num+' L infrared image exists')
    return 'ok'


def ploteye2(h5f,patient='070',eyeside='L',t=0,gridpts=640,tness='GCIPL',
            opacity1=1.0,opacity2=1.0,cutval=5.0):
    
    patient_num = patient
    localangle = np.abs(angles.loc[angles['AGPS']=='AGPS'+patient_num,'angle'].iloc[0]) *np.pi/180
    totalheight = 30*np.sin(localangle) + 25*np.cos(localangle)
    totalwidth = 30*np.cos(localangle)

    yfac=gridpts
    xfac=gridpts
    centerrangemin = int(0.2*gridpts)
    centerrangemax = int(0.8*gridpts)

    tstr = str(t)
    if len(h5f.keys()) == 0:
        return('Patient '+patient+
               ',all NaNs (or unable to be interpolated)')
    elif eyeside not in h5f[patient].keys():
        return('Patient '+patient+' eyeside '+eyeside+
               ',all NaNs (or unable to be interpolated)')
    elif tstr not in h5f[patient+'/'+eyeside].keys():
        return('Patient '+patient+' eyeside '+eyeside+' filenum '+tstr+
               ', all NaNs (or unable to be interpolated)')
    elif str(yfac) not in h5f[patient+'/'+eyeside+'/'+tstr].keys():
        return('Patient '+patient+' eyeside '+eyeside+' filenum '+tstr+' resolution '+str(yfac)+
               ', all NaNs (or unable to be interpolated)')
    else:

        if eyeside == 'R':
            if os.path.exists('infrared/'+patient_num+'_OD_SPEC.JPG'):
                im = plt.imread('infrared/'+patient_num+'_OD_SPEC.JPG')
            elif os.path.exists('infrared/'+patient_num+'_OD_SPEC.jpg'):
                im = plt.imread('infrared/'+patient_num+'_OD_SPEC.jpg')
            elif os.path.exists('infrared/'+patient_num+'_OD_SPEC.jpg'):
                im = plt.imread('infrared/'+patient_num+'_OD_SPEC.jpeg')
            else:
                clear_output(wait=True)
                #plt.figure().clf()
                print('No '+patient_num+' R infrared image exists')
                return('No '+patient_num+' R infrared image exists')
        else:
            if os.path.exists('infrared/'+patient_num+'_OS_SPEC.JPG'):
                im = plt.imread('infrared/'+patient_num+'_OS_SPEC.JPG')
            elif os.path.exists('infrared/'+patient_num+'_OS_SPEC.jpg'):
                im = plt.imread('infrared/'+patient_num+'_OS_SPEC.jpg')
            elif os.path.exists('infrared/'+patient_num+'_OS_SPEC.jpeg'):
                im = plt.imread('infrared/'+patient_num+'_OS_SPEC.jpeg')
            else:
                clear_output(wait=True)
                #plt.figure().clf()
                print('No '+patient_num+' L infrared image exists')
                return('No '+patient_num+' L infrared image exists')

        fig,ax = plt.subplots(1,4)
        fig.set_size_inches(4*5,5*totalheight/totalwidth)
        ax[0].set_axis_off()
        ax[0].imshow(im, aspect='auto',extent=[0, totalwidth, 0, totalheight])

        tstr_base = str(0)
        
        tmpdata_base = h5f[patient+'/'+eyeside+'/'+tstr_base+'/'+str(yfac)+'/'+str(xfac)+'/dataset'][:]
        d = pd.to_datetime(h5f[patient+'/'+eyeside+'/'+tstr_base+'/'+str(yfac)+'/'+str(xfac)+'/dataset'].attrs['date'],format='%Y-%m-%d')
        patient119Rv0 = patientdf.loc[patientdf['EID']=='AGPS'+patient_num+'_'+eyeside,:]
        m = closest(patient119Rv0['VISIT_DATE'],d)

        divider = make_axes_locatable(ax[0])
        cax = divider.append_axes("right", size="5%", pad=0.15)
        cax2 = divider.append_axes("right", size="5%", pad=0.6)

        p2 = ax[0].imshow(tmpdata_base,
                  aspect=1,
                  origin='lower', cmap='inferno', #extent=[0, 24, 0, 24], 
                       extent=[(totalwidth-24)/2,
                                          totalwidth-(totalwidth-24)/2,
                                          (totalheight-24)/2,
                                          totalheight-(totalheight-24)/2], 
                  vmin=np.amin(tmpdata_base[centerrangemin:centerrangemax,centerrangemin:centerrangemax]),
                       vmax=np.amax(tmpdata_base[centerrangemin:centerrangemax,centerrangemin:centerrangemax]),
                   alpha=opacity1)
        plt.colorbar(p2, cax=cax)   
        p3 = ax[0].scatter(x=totalwidth/2.0-0.5 + displacements['xpos'],
                    y=totalheight/2.0 - displacements['ypos'],
                    c=patient119Rv0.iloc[m][-68:] - patient119Rv0.iloc[0][-68:],
                    cmap = 'bwr',
                    s=25,
                    alpha=opacity2,
                    edgecolors='w',
                    vmax=15,vmin=-15);
        plt.colorbar(p3,cax=cax2)
        ax[0].text(36,0.45*totalheight,'OCT',rotation='vertical',fontsize=12)
        ax[0].text(47,0.1*totalheight,'VF difference from baseline',rotation='vertical',fontsize=12)
        ax[0].set_title('Baseline\n' + # = '+h5f[patient+'/'+eyeside+'/'+tstr+'/'+str(yfac)+'/'+str(xfac)+'/dataset'].attrs['date']+'\n'+
                     'Date = '+patient119Rv0['VISIT_DATE'].iloc[0].strftime('%Y-%m-%d'))

        ax[1].set_axis_off()
        ax[1].imshow(im, aspect='auto',extent=[0, totalwidth, 0, totalheight])

        tmpdata = h5f[patient+'/'+eyeside+'/'+tstr+'/'+str(yfac)+'/'+str(xfac)+'/dataset'][:]
        d = pd.to_datetime(h5f[patient+'/'+eyeside+'/'+tstr+'/'+str(yfac)+'/'+str(xfac)+'/dataset'].attrs['date'],format='%Y-%m-%d')
        patient119Rv0 = patientdf.loc[patientdf['EID']=='AGPS'+patient_num+'_'+eyeside,:]
        m = closest(patient119Rv0['VISIT_DATE'],d)

        divider = make_axes_locatable(ax[1])
        cax = divider.append_axes("right", size="5%", pad=0.15)
        cax2 = divider.append_axes("right", size="5%", pad=0.6)

        p2 = ax[1].imshow(tmpdata,
                  aspect=1,
                  origin='lower', cmap='inferno', #extent=[0, 24, 0, 24], 
                       extent=[(totalwidth-24)/2,
                                          totalwidth-(totalwidth-24)/2,
                                          (totalheight-24)/2,
                                          totalheight-(totalheight-24)/2], 
                  vmin=np.amin(tmpdata[centerrangemin:centerrangemax,centerrangemin:centerrangemax]),
                       vmax=np.amax(tmpdata[centerrangemin:centerrangemax,centerrangemin:centerrangemax]),
                   alpha=opacity1)
        plt.colorbar(p2, cax=cax)   
        p3 = ax[1].scatter(x=totalwidth/2.0-0.5 + displacements['xpos'],
                    y=totalheight/2.0 - displacements['ypos'],
                    c=patient119Rv0.iloc[m][-68:] - patient119Rv0.iloc[0][-68:],
                    cmap = 'bwr',
                    s=25,
                    alpha=opacity2,
                    edgecolors='w',
                    vmax=15,vmin=-15);
        plt.colorbar(p3,cax=cax2)
        ax[1].text(36,0.45*totalheight,'OCT',rotation='vertical',fontsize=12)
        ax[1].text(47,0.1*totalheight,'VF difference from baseline',rotation='vertical',fontsize=12)
        ax[1].set_title('OCT Date = '+h5f[patient+'/'+eyeside+'/'+tstr+'/'+str(yfac)+'/'+str(xfac)+'/dataset'].attrs['date']+'\n'+
                     'VF Date = '+patient119Rv0['VISIT_DATE'].iloc[m].strftime('%Y-%m-%d')+'\n'+
                    'VF Base Date = '+patient119Rv0['VISIT_DATE'].iloc[0].strftime('%Y-%m-%d'))

        ax[2].set_axis_off()
        ax[2].imshow(im, aspect='auto',extent=[0, totalwidth, 0, totalheight])

        masterdiff = (tmpdata[:] - tmpdata_base[:])
        masterdiff[np.isnan(masterdiff)] = 0
        abovecut = masterdiff > (-1*cutval)
        masterdiff[abovecut] = 0

        d = pd.to_datetime(h5f[patient+'/'+eyeside+'/'+tstr+'/'+str(yfac)+'/'+str(xfac)+'/dataset'].attrs['date'],format='%Y-%m-%d')
        patient119Rv0 = patientdf.loc[patientdf['EID']=='AGPS'+patient_num+'_'+eyeside,:]
        m = closest(patient119Rv0['VISIT_DATE'],d)

        divider = make_axes_locatable(ax[2])
        cax = divider.append_axes("right", size="5%", pad=0.15)
        cax2 = divider.append_axes("right", size="5%", pad=0.6)

        p2 = ax[2].imshow(abs(masterdiff),
                  aspect=1,
                  origin='lower', cmap='inferno', #extent=[0, 24, 0, 24], 
                       extent=[(totalwidth-24)/2,
                                          totalwidth-(totalwidth-24)/2,
                                          (totalheight-24)/2,
                                          totalheight-(totalheight-24)/2], 
                  vmin=0,
                       vmax=20,
                   alpha=opacity1)
        plt.colorbar(p2, cax=cax, ticks=[0,4,8,12,16,20])   
        p3 = ax[2].scatter(x=totalwidth/2.0-0.5 + displacements['xpos'],
                    y=totalheight/2.0 - displacements['ypos'],
                    c=patient119Rv0.iloc[m][-68:] - patient119Rv0.iloc[0][-68:],
                    cmap = 'bwr',
                    s=25,
                    alpha=opacity2,
                    edgecolors='w',
                    vmax=15,vmin=-15);
        plt.colorbar(p3,cax=cax2)
        ax[2].text(36,0.45*totalheight,'OCT',rotation='vertical',fontsize=12)
        ax[2].text(47,0.1*totalheight,'VF difference from baseline',rotation='vertical',fontsize=12)
        ax[2].set_title('Lower Vals')

        
        ax[3].set_axis_off()
        ax[3].imshow(im, aspect='auto',extent=[0, totalwidth, 0, totalheight])

        masterdiff = (tmpdata[:] - tmpdata_base[:])
        masterdiff[np.isnan(masterdiff)] = 0
        abovecut = masterdiff < (cutval)
        masterdiff[abovecut] = 0

        d = pd.to_datetime(h5f[patient+'/'+eyeside+'/'+tstr+'/'+str(yfac)+'/'+str(xfac)+'/dataset'].attrs['date'],format='%Y-%m-%d')
        patient119Rv0 = patientdf.loc[patientdf['EID']=='AGPS'+patient_num+'_'+eyeside,:]
        m = closest(patient119Rv0['VISIT_DATE'],d)

        divider = make_axes_locatable(ax[3])
        cax = divider.append_axes("right", size="5%", pad=0.15)
        cax2 = divider.append_axes("right", size="5%", pad=0.6)

        p2 = ax[3].imshow(abs(masterdiff),
                  aspect=1,
                  origin='lower', cmap='inferno', #extent=[0, 24, 0, 24], 
                       extent=[(totalwidth-24)/2,
                                          totalwidth-(totalwidth-24)/2,
                                          (totalheight-24)/2,
                                          totalheight-(totalheight-24)/2], 
                  vmin=0,
                       vmax=20,
                   alpha=opacity1)
        plt.colorbar(p2, cax=cax, ticks=[0,4,8,12,16,20])   
        p3 = ax[3].scatter(x=totalwidth/2.0-0.5 + displacements['xpos'],
                    y=totalheight/2.0 - displacements['ypos'],
                    c=patient119Rv0.iloc[m][-68:] - patient119Rv0.iloc[0][-68:],
                    cmap = 'bwr',
                    s=25,
                    alpha=opacity2,
                    edgecolors='w',
                    vmax=15,vmin=-15);
        plt.colorbar(p3,cax=cax2)
        ax[3].text(36,0.45*totalheight,'OCT',rotation='vertical',fontsize=12)
        ax[3].text(47,0.1*totalheight,'VF difference from baseline',rotation='vertical',fontsize=12)
        ax[3].set_title('Higher Vals')

        
        
        clear_output(wait=True)
        return fig
        return 'done'

