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

import subprocess

list_files = subprocess.run(["wget","-O","-","https://github.com/benjum/eye-to-share-d/archive/main.tar.gz","|","tar xz"])
list_files = subprocess.run(["mv","eye-to-share-d-main/00layers","./data"])
list_files = subprocess.run(["rm","-r","eye-to-share-d-main"])

patientdf = pd.read_csv('TD_vf_vals.csv')
displacements = pd.read_csv('VF_Locations.csv')

displacements['xpos']=displacements['x_morphed']
displacements['ypos']=displacements['y_morphed']

#x = %system ls 'data/'
from os import listdir
from os.path import isfile, join
x = [f for f in listdir('data') if isfile(join('data',f))]

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
		
h5f = h5py.File(h5files[0],'r')
for i in h5f['001/L'].keys():
    yunique = [int(j) for j in h5f['001/L/'+str(i)+''].keys()]
    yunique.sort(reverse=True)
    xunique = yunique.copy()
    break
h5f.close()

angles = pd.read_csv('angles_Complete.csv')

def closest(lst, K):
    x = []
    for i in range(len(lst)):
        x.append(abs(lst.iloc[i]-K))
    m = x.index(min(x))
    return m

patientdf.loc[:,'VISIT_DATE'] = pd.to_datetime(patientdf.loc[:,'VISIT_DATE'],format='%m/%d/%y')


# h5f = h5py.File(h5files[0],'r')
def ploteye(h5f,patient='001',eyeside='L',t=0,gridpts=640,tness='GCIPL',
            opacity1=1.0,opacity2=1.0):

    
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
            else:
                return('No '+patient_num+' R infrared image exists')
        else:
            if os.path.exists('infrared/'+patient_num+'_OS_SPEC.JPG'):
                im = plt.imread('infrared/'+patient_num+'_OS_SPEC.JPG')
            elif os.path.exists('infrared/'+patient_num+'_OS_SPEC.jpg'):
                im = plt.imread('infrared/'+patient_num+'_OS_SPEC.jpg')
            else:
                return('No '+patient_num+' L infrared image exists')

        fig = plt.figure(frameon=False)
        fig.set_size_inches(5,5*totalheight/totalwidth)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(im, aspect='auto',extent=[0, totalwidth, 0, totalheight])


        tmpdata = h5f[patient+'/'+eyeside+'/'+tstr+'/'+str(yfac)+'/'+str(xfac)+'/dataset'][:]
        d = pd.to_datetime(h5f[patient+'/'+eyeside+'/'+tstr+'/'+str(yfac)+'/'+str(xfac)+'/dataset'].attrs['date'],format='%Y-%m-%d')
        patient119Rv0 = patientdf.loc[patientdf['EID']=='AGPS'+patient_num+'_'+eyeside,:]
        m = closest(patient119Rv0['VISIT_DATE'],d)

        ax.imshow(tmpdata,
                  aspect=1,
                  origin='lower', cmap='inferno', #extent=[0, 24, 0, 24], 
                       extent=[(totalwidth-24)/2,
                                          totalwidth-(totalwidth-24)/2,
                                          (totalheight-24)/2,
                                          totalheight-(totalheight-24)/2], 
                  vmin=np.amin(tmpdata[centerrangemin:centerrangemax,centerrangemin:centerrangemax]),
                       vmax=np.amax(tmpdata[centerrangemin:centerrangemax,centerrangemin:centerrangemax]),
                   alpha=opacity1)

        p3 = ax.scatter(x=totalwidth/2.0-0.5 + displacements['xpos'],
                    y=totalheight/2.0 - displacements['ypos'],
                    c=patient119Rv0.iloc[m][-68:] - patient119Rv0.iloc[0][-68:],
                    cmap = cmap1,
                    s=100,
                    alpha=opacity2,
                    edgecolors='w',
                    vmax=15,vmin=-15);
        plt.colorbar(p3)
        ax.set_title('OCT Date = '+h5f[patient+'/'+eyeside+'/'+tstr+'/'+str(yfac)+'/'+str(xfac)+'/dataset'].attrs['date']+'\n'+
                     'VF Date = '+patient119Rv0['VISIT_DATE'].iloc[m].strftime('%Y-%m-%d')+'\n'+
                    'VF Raw Value = '+patient119Rv0['VISIT_DATE'].iloc[0].strftime('%Y-%m-%d'))
        display(fig)
        clear_output(wait=True)
        
# d1 = widgets.Dropdown(options=numsunique,description='Patient')
# d2 = widgets.Dropdown(options=['L','R'],description='Eyeside')
# d3 = widgets.Dropdown(options=thickunique,description='Layer',value='GCIPL')
# d4 = widgets.Dropdown(options=yunique,description='Grid Pts')
# s1 = widgets.IntSlider(min=0,
#                        max=len(h5f['001/L'].keys())-1)

# def valued1_changed(change):
#     global h5f
#     h5f.close()
#     h5f = h5py.File('data/'+change.new+'-'+d3.value+'.hdf5','r')
#     if len(h5f.keys()) == 0:
#         s1.max = 0
#     elif d2.value not in h5f[change.new].keys():
#         s1.max = 0
#     else:
#         s1.max = max(0,len(h5f[change.new+'/'+d2.value].keys())-1)
#     s1.value = 0
# d1.observe(valued1_changed, 'value')

# def valued2_changed(change):
#     if len(h5f.keys()) == 0:
#         s1.max = 0
#     elif change.new not in h5f[d1.value].keys():
#         s1.max = 0
#     else:
#         s1.max = max(0,len(h5f[d1.value+'/'+change.new].keys())-1)
#     s1.value = 0
# d2.observe(valued2_changed, 'value')

# def valued3_changed(change):
#     global h5f
#     h5f.close()
#     h5f = h5py.File('data/'+d1.value+'-'+change.new+'.hdf5','r')
#     if len(h5f.keys()) == 0:
#         s1.max = 0
#     elif d2.value not in h5f[d1.value].keys():
#         s1.max = 0
#     else:
#         s1.max = max(0,len(h5f[d1.value+'/'+d2.value].keys())-1)
#     s1.value = 0
# d3.observe(valued3_changed, 'value')

# interact(ploteye, 
#          patient=d1, 
#          eyeside=d2,
#          tness=d3,
#          gridpts=d4,
#          t=s1,opacity1=(0.0,1.0,0.1),opacity2=(0.0,1.0,0.1));
