#%% IMPORTS
import os
import numpy as np
#%matplotlib widget  #requires: jupyter-notebook and https://github.com/matplotlib/ipympl (doesn't work in jupyter lab)
from matplotlib import pyplot as plt    
import pandas as pd
import flow2P
import re

#%% MULTIPLE ANIMAL SCRIPT
cohortDir = r'D:\Data\Patrick\processed2P\ZnT3_5-52kHz_G2'
animals = re.findall(r'[A-Z]{2}\d{4}',
    ' '.join(os.listdir(cohortDir)))

outs = [flow2P.animal2P(os.path.join(cohortDir,A),
    os.path.join(cohortDir,A,'suite2p','plane0')) for A in animals]
[out.TRFcalc() for out in outs]

dfStimCohort = pd.concat([out.dfStim for out in outs])
dfStimCohort.index.names
# %%

uROI = dfStimCohort[dfStimCohort['BFoctDiff'].abs()<1].groupby(\
    ['animal','roi','treatment','PTonset_s','contrast'])['spks']\
        .apply(np.mean).groupby(['treatment','PTonset_s','contrast'])\
            .apply(np.mean)

for group_name, df_group in uROI.groupby(['treatment','PTonset_s']):
    plt.figure()
    print(group_name)
    plt.plot(uROI[(*group_name, 'low')])
    plt.plot(uROI[(*group_name, 'high')])
    plt.title(group_name)
    plt.legend(('low','high'))
    plt.ylabel('Firing Rate (Hz)')
    plt.xticks(np.arange(0,75,10),(np.arange(0,75,10)/5)-4)
    plt.xlabel('Time (s)')