# -*- coding: utf-8 -*-
"""
Created on 7/8/2020

@author: Patrick Cody | pac94@pitt.edu | pcody1@gmail.com
"""

# %% IMPORT LIBRARIES

import os
import numpy as np
#%matplotlib widget  #requires: jupyter-notebook and https://github.com/matplotlib/ipympl (doesn't work in jupyter lab)
from matplotlib import pyplot as plt    
import pandas as pd
import flow2P
#for updating flow2P:
#import importlib
#importlib.reload(flow2P)


#%% INPUT DATA DIRECTORY
tifDir = r'D:\Data\Patrick\processed2P\ZnT3_5-52kHz_G2\ZN0001'
# tifDir = r'D:\Data\Patrick\processed2P\CaMKII_5-25kHz\AA0104'
suite2PoutputDir = os.path.join(tifDir,'suite2p','plane0')


#%% JUST TRF
# animal = flow2P.getAnimal(tifDir)
# iscell = np.load(os.path.join(suite2PoutputDir,'iscell.npy'))
# dfTif,dfMap,_ = flow2P.getTifLegends(tifDir,animal)

# TRF = flow2P.animalTRF(dfMap,dfTif,iscell,suite2PoutputDir) #defaults to pre, TRF is TRFpre
# TRFpre = flow2P.animalTRF(dfMap,dfTif,iscell,suite2PoutputDir,treatment='pre')
# TRFpost = flow2P.animalTRF(dfMap,dfTif,iscell,suite2PoutputDir,treatment='post')


#%% analyze animal with stim output
animalOutput = flow2P.animal2P(tifDir,suite2PoutputDir)

#%% output plots
#plot average traces across ROI | PUT THIS IN flow2P as a plotting function?
groups = ['roi','treatment','PTonset_s','PTfreq_Hz','contrast']

#average across ROI
uROI = animalOutput.dfStim.groupby(groups)['spks'].apply(np.mean).\
    groupby(np.asarray(groups)[~np.isin(groups,'roi')].tolist()).apply(np.mean)

# for group_name, df_group in uROI.groupby(['treatment','PTonset_s']):
for group_name, df_group in uROI.groupby(np.asarray(groups)[~np.isin(groups,['roi','contrast'])].tolist()):
    plt.figure()
    print(group_name)
    plt.plot(uROI[(*group_name, 'low')])
    plt.plot(uROI[(*group_name, 'high')])
    plt.title(group_name)
    plt.legend(('low','high'))
    plt.ylabel('Firing Rate (Hz)')
    plt.xticks(np.arange(0,75,10),(np.arange(0,75,10)/5)-4)
    plt.xlabel('Time (s)')

#%%    
animalOutput.TRFcalc()

#%% get common BF:
animalOutput.TRF.freq[np.argmax([sum(np.isin(animalOutput.TRF.BF,fq))\
                                 for fq in animalOutput.TRF.freq])]
#or
from scipy import stats
stats.mode(np.asarray(animalOutput.TRF.BF))

#%% plot traces from cells with BF within 1 octave from pure-tone stimulus
octLim = 1

uROIwithin = animalOutput.dfStim[animalOutput.dfStim['BFoctDiff'].abs()<octLim]\
    .groupby(groups)['spks'].apply(np.mean)\
        .groupby(np.asarray(groups)[~np.isin(groups,'roi')].tolist()).apply(np.mean)

for group_name, df_group in uROIwithin\
    .groupby(np.asarray(groups)[~np.isin(groups,['roi','contrast'])].tolist()):
    
    plt.figure()
    print(group_name)
    plt.plot(uROIwithin[(*group_name, 'low')])
    plt.plot(uROIwithin[(*group_name, 'high')])
    plt.title(group_name)
    plt.legend(('low','high'))
    plt.ylabel('Firing Rate (Hz)')
    plt.xticks(np.arange(0,75,10),(np.arange(0,75,10)/5)-4)
    plt.xlabel('Time (s)')
# %%
