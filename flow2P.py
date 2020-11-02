# -*- coding: utf-8 -*-
"""
Created on Created on 7/8/2020

@author: Patrick Cody | pac94@pitt.edu | pcody1@gmail.com
"""

import os
import regex
from scipy.io import loadmat
import glob
import numpy as np
from itertools import compress    
from matplotlib import pyplot as plt
import pandas as pd




_stimParamRegex = regex.compile(r'(?|((?P<PTfreq_Hz>\d{4,5})Hz)|(stim_(?P<PTfreq_Hz>\d{1,2})kHz))'
                                   r'_(?P<PTampl_dB>\d{2})dB_(?P<PTpulseWidth_ms>\d{2,3})ms'
                                   r'(?|(_at_(?P<PTonset_s>\d)s)|(_at_(?P<PTonset_s>\dpt\d)s))')

_mapPulseParamRegex = regex.compile(r'(?P<freq_Hz>\d{4,5})Hz_(?P<ampl_dB>\d{2})dB_TestTone'
                                       r'_(?P<pulseWidth_ms>\d{2,3})msPulse_at_'
                                       r'((?P<onset_s>\d{1}\.\d{1,2})s|(?P<onset_s>\d{1,2})s)'
                                       r'_(?P<ramp_ms>\d{1,3})msRamp_(?<stimLen_ms>\d{2,5})'
                                       r'msTotal_Fs(?<sampleRate_kHz>\d{2,4})kHz')

_pp = ['pre','post']




def getAnimal(tifDir):
    """Infer animal name from tif directory."""
    
    if regex.search(r'[A-Z]{2}\d{4}',os.path.split(tifDir)[-1]) is not None:
        animal = regex.search(r'[A-Z]{2}\d{4}',os.path.split(tifDir)[-1])[0]
    else:
        animal = input('Enter animal name:\n')
    return animal



def calcOctDiff(freqA,freqB):
    """Calculates octaves from A to B"""
    
    return np.log2(np.asarray(freqA)/np.asarray(freqB))




def extractStimParams(pulseStruct,triggerParamsStruct):
    """Takes pulse and param struct from a tif .mat stim stamp 
    and outputs stimulus parameters as a dict"""
    
    stimParams = {'rawPulse': str(pulseStruct['pulsename']),
                 'Pulse': str(pulseStruct['pulsename'])[:-2],\
                 'delayToEphusTrig_s': float(triggerParamsStruct['stimDelay'])}

    stimParams.update(regex.search(_stimParamRegex,stimParams['Pulse']).groupdict())
        
    stimParams['PTfreq_Hz'] = int(stimParams['PTfreq_Hz'])*1000 if \
                              len(stimParams['PTfreq_Hz'])<3 else int(stimParams['PTfreq_Hz'])
    stimParams['PTonset_s'] = float(stimParams['PTonset_s'].replace('pt','.'))
         
    if 'DRC' in stimParams['Pulse']:           
            stimParams.update(regex.search(r'(?P<DRCamplRange_dB>\d{2}-\d{2})dB'
                                           r'_(?P<ephusPulseLen_s>\d{1,2})s_',
                                           stimParams['Pulse']).groupdict()) 
            stimParams['DRCamplDelta_dB'] = abs(eval(stimParams['DRCamplRange_dB']))
            stimParams['DRClen_ms'] = int(regex.search(r'(?P<DRClen>\d*)msDRC',
                                          stimParams['Pulse']).group('DRClen'))
            try:
                stimParams['DRCbandwidth_kHz'] = regex.search(r'(?P<kHzRng>\d-\d{2})kHz',
                                                              stimParams['Pulse']).group('kHzRng').replace('-',',')
        
            except AttributeError: 
                if '45-55' in stimParams['DRCamplRange_dB'] or '35-65' in stimParams['DRCamplRange_dB']:
                        stimParams['DRCbandwidth_kHz'] = '5,25'
                else:
                    raise ValueError('Could not determine DRC stimulus bandwidth')                          
    return stimParams




def extractMapParams(pulseStruct,triggerParamsStruct):
        """Takes pulse and param struct from tif .mat stim stamp 
        and outputs map stim params dict for multiple pulses.
        Assumes multiple pulses in 'pulse' struct/dict."""
        
        d = dict()
        for pNo in range(0,pulseStruct['pulsename'].shape[0]):
            for i,a in regex.search(_mapPulseParamRegex,
                                    pulseStruct['pulsename'][pNo]).groupdict().items():
                d.setdefault(i,[]).append(a)
        d['delayToEphusTrig_s'] = float(triggerParamsStruct['stimDelay'])
        d['ISI'] = float(triggerParamsStruct['ISI'])
        d['totalPulses'] = int(triggerParamsStruct['totalPulses'])
        return d

    
    
    
def getTifStimParams(tifFullFilePath):
    """Gets _Pulses.mat file assicaited with .tif provided, loads .mat, 
    discerns b/w map and stim pulse, and extracts respective pulses and parameters."""
    #have to loadmat here b/c need to discern map from stim params

    tempMat = loadmat(tifFullFilePath.replace('.tif','_Pulses.mat'),squeeze_me=True)
    if tempMat['pulse'].size>1: #map
        stim = tempMat['pulse']['pulseset'][0]
        stimParamDict = extractMapParams(tempMat['pulse'],tempMat['params'])
    else:
        stimParamDict = extractStimParams(tempMat['pulse'],tempMat['params'])
        stim = stimParamDict['rawPulse']
    return (stimParamDict,stim)




def getTifFRnFrames(tifDir):
        """Gets frame rate and number of frames from all .tifs in tifDir."""
        
        from ScanImageTiffReader import ScanImageTiffReader
        tifFiles = glob.glob(os.path.join(tifDir,'*.tif'))
        tifFrameRate = []
        nFrames = []
        with ScanImageTiffReader(tifFiles[0]) as reader:
            if not reader.metadata():
                SIv=3
                frRegex = regex.compile(r'frameRate=(?P<frameRate>\d.\d+)') #for ScanImage v3
            else:
                SIv=5
                frRegex = regex.compile(r'hRoiManager.scanFrameRate = (?P<frameRate>\d+)') #for ScanImage v5+
                
        for tifFile in tifFiles:           
            with ScanImageTiffReader(tifFile) as reader:
                if SIv==3:
                    tifFrameRate.append(regex.search(frRegex,reader.description(0)).group('frameRate'))
                elif SIv==5:
                    tifFrameRate.append(regex.search(frRegex,reader.metadata()).group('frameRate'))
                    
                nFrames.append(reader.shape()[0])

        return list(zip(tifFiles,np.array(tifFrameRate,dtype=float).round(1),nFrames))        
    


    
def getTreatmentDrug(tifDir):
        """Tries to infer drug treatment from tif directory."""
        
        treatments = ['ZX1','ACSF']
        try:  #look for presence of treatment name in directory
            drug = [i for (i,v) in zip(treatments,
                                       [any([treatment in item for item in os.listdir(tifDir)])\
                                        for treatment in treatments]) if v][0]
        
        except (IndexError, AttributeError):
            drug = input('Enter treatment name:\n')
            
        return drug
    
    
    
    
def getTifTreatment(tifDir,maxTimeBWtif=18,afterXtifs=10):
        """Infers pre/post treatment assicaited with each .tif file
        based upon time between successive tif files."""
        
        drug = getTreatmentDrug(tifDir)
        
        #if it has been more than 18 minutes between tif files 
        # and it's at least x tif files from first consider next post drug
        mTime = []
        treatmentKey = []
        ppID = 0
        tifNo = 0
        try:
            for entry in os.scandir(tifDir):
                if entry.is_file() and '.tif' in entry.name:
                    tifNo += 1
                    mTime.append(entry.stat().st_mtime)
        
                    if (len(mTime)>1 and (mTime[-1]-mTime[-2])/60)<maxTimeBWtif or tifNo<afterXtifs:
                        treatmentKey.append(_pp[ppID] + drug)
                    else:
                        ppID += 1
                        treatmentKey.append(_pp[ppID] + drug)
        
        except (IndexError):
            ppID = 0
            treatmentKey = []

            if len(glob.glob(os.path.join(tifDir,'lastPreTreatmentTif_*'))) != 0:
                lastPreTif = os.path.basename(glob.glob(os.path.join(tifDir,'lastPreTreatmentTif_*'))[0].replace('lastPreTreatmentTif_',''))
            else:
                lastPreTif = input('Enter last pre-treatment .tif file name:\n').replace('.tif','')
                from pathlib import Path
                Path(os.path.join(tifDir,'lastPreTreatmentTif_' + lastPreTif)).touch()
            
            for entry in os.scandir(tifDir):
                if entry.is_file() and '.tif' in entry.name:
                            
                    if entry.name == lastPreTif + '.tif':
                        treatmentKey.append(_pp[ppID] + drug)
                        ppID += 1
                    else:
                        treatmentKey.append(_pp[ppID] + drug)

        return treatmentKey



    
def getTifLegends(tifDir,animal):
        """Generates dataFrames of relevant params for tifs using getTifTreatment and getTifFRnFrames.
        Splits .tifs into dataFrames respective to mapping and experimental stimuli"""
        
        tifFRnFrames = getTifFRnFrames(tifDir)
        tifTreatments = getTifTreatment(tifDir)
        
        tifLegend = dict()
        tifStimLegend = dict()
        tifMapLegend = dict()
        mapIDX = []
        stimIDX = []
        for i,(tif,FR,nFrames) in enumerate(tifFRnFrames):
            tifLegend.setdefault('animal',[]).append(animal)
            tifLegend.setdefault('treatment',[]).append(tifTreatments[i])
            tifLegend.setdefault('fileName',[]).append(tif)
            tifLegend.setdefault('frameRate',[]).append(FR)
            tifLegend.setdefault('nFrames',[]).append(nFrames)
            for params,stim in [getTifStimParams(tif)]:
                tifLegend.setdefault('stim',[]).append(stim)
                if 'map' in stim.lower():
                    mapIDX.append(i)
                    tifMapLegend.setdefault('pulseSet',[]).append(stim)
                    for key,value in params.items():
                        tifMapLegend.setdefault(key,[]).append(value)
                else:
                    stimIDX.append(i)
                    for key,value in params.items():
                        tifStimLegend.setdefault(key,[]).append(value)
                        
        for key,value in tifLegend.items():
            tifMapLegend[key] = np.array(tifLegend[key])[mapIDX].tolist()
            tifStimLegend[key] = np.array(tifLegend[key])[stimIDX].tolist()           
            
        return pd.DataFrame.from_dict(tifLegend),\
               pd.DataFrame.from_dict(tifMapLegend),\
               pd.DataFrame.from_dict(tifStimLegend)




def dBdelta2lowHigh(dfStim,insertColLoc = 9):
    """Adds contrast column to dfStim in place respective to DRCamplDelta_dB"""
    
    delta_to_lohi = {np.min(dfStim['DRCamplDelta_dB']): 'low',
                 np.max(dfStim['DRCamplDelta_dB']): 'high'}
    dfStim.insert(insertColLoc,'contrast',dfStim['DRCamplDelta_dB'].map(delta_to_lohi))


    
    
def ROIasMultiIndex(dfStim,iscell):
        """Split out ROI as a separate dimension (as a 'MultiIndex') to make slicing/indexing easier. 
        Do this AFTER adding suite2p data to table (if adding suite2P data)."""
        
        #via: https://stackoverflow.com/questions/45846765/efficient-way-to-unnest-explode-multiple-list-columns-in-a-pandas-dataframe 
        #set index as those columns exluded from explode
        dfNew = dfStim.set_index(dfStim.columns[~np.isin(dfStim.columns,['F','Fneu','spks'])].tolist()).apply(pd.Series.explode).reset_index()
        dfNew['roi'] = np.tile(np.arange(0,int(np.sum(iscell[:,0]))),int(len(dfNew)/np.sum(iscell[:,0])))
        
        return dfNew.set_index(['fileName','roi'])     

    
    
    
#decide whether to use iscell or load from file every time
def suite2PvarByTif(dfTif,iscell,suite2PoutputDir,suite2PoutputVar = 'spks'):
    """Splits suite2P output vector from suite2PoutputVar
    into lists of frames respective to tif files."""
    
    frameSpan = zip(dfTif['nFrames'].cumsum()-dfTif['nFrames'],dfTif['nFrames'].cumsum())

    suite2PoutputVar = np.load(os.path.join(suite2PoutputDir,suite2PoutputVar + '.npy'))[iscell[:,0].nonzero()[0],:]
        
    suite2PoutputByTif = []
    for a,b in frameSpan:
        #a[start:stop]  # items start through stop-1 see: https://stackoverflow.com/questions/509211/understanding-slice-notation
        suite2PoutputByTif.append(suite2PoutputVar[:,a:b])
    
    return suite2PoutputByTif




def suite2Poutput2dfStim(dfStim,dfTif,iscell,suite2PoutputDir,suite2PoutputVars=['spks']):
    """Adds suite2PoutputVars to dfStim organized by tifs."""
    
    for varI in suite2PoutputVars:
        dfStim[varI] = list(compress(suite2PvarByTif(dfTif,iscell,suite2PoutputDir,varI),\
                                     ~dfTif['stim'].str.contains('map',case=0)))
        
        
        
        
class animal2P:
    """Container for experiment output from a single animal."""
    def __init__(self,tifDir,suite2PoutputDir=None):
        self.tifDir = tifDir
        self.animal = getAnimal(self.tifDir)
        self.drug = getTreatmentDrug(self.tifDir)
        self.suite2PoutputDir = suite2PoutputDir
        [self.dfTif,
         self.dfMap,
         self.dfStim] = getTifLegends(self.tifDir,self.animal)
        dBdelta2lowHigh(self.dfStim)
        if self.suite2PoutputDir is not None:
            self.iscell = np.load(os.path.join(suite2PoutputDir,'iscell.npy'))
            suite2Poutput2dfStim(self.dfStim,self.dfTif,self.iscell,self.suite2PoutputDir)
        self.dfStim = ROIasMultiIndex(self.dfStim,self.iscell)
        
    def TRFcalc(self,fLookPostOnset = 3):
        if np.logical_and(np.any(self.dfMap['treatment'].str.contains('pre')),\
                          np.any(self.dfMap['treatment'].str.contains('post'))):
            self.TRF = animalTRF(self.dfMap,self.dfTif,self.iscell,self.suite2PoutputDir,fLookPostOnset,'pre')
            self.TRFpost = animalTRF(self.dfMap,self.dfTif,self.iscell,self.suite2PoutputDir,fLookPostOnset,'post')
            self.dfStim = pd.DataFrame.join(pd.DataFrame(np.transpose([self.dfStim.index.get_level_values(1).unique(),self.TRF.BF,self.TRFpost.BF]),\
                columns=['roi','BF','BFpost']).set_index('roi'),self.dfStim)  
            self.dfStim['BFoctDiff'] = calcOctDiff(self.dfStim['BF'],self.dfStim['PTfreq_Hz'])
            self.dfStim['BFoctDiffPost'] = calcOctDiff(self.dfStim['BFpost'],self.dfStim['PTfreq_Hz'])
            
        else:
            self.TRF = animalTRF(self.dfMap,self.dfTif,self.iscell,self.suite2PoutputDir,fLookPostOnset)
            self.dfStim = pd.DataFrame.join(pd.DataFrame(np.transpose([self.dfStim.index.get_level_values(1).unique(),self.TRF.BF]),\
                columns=['roi','BF']).set_index('roi'),self.dfStim)
            self.dfStim['BFoctDiff'] = calcOctDiff(self.dfStim['BF'],self.dfStim['PTfreq_Hz'])
    

    
#maybe have iscell come directly from suite2PoutputDir
class animalTRF:
    """Container for a TRF from one mapping session."""
    
    def __init__(self,dfMap,dfTif,iscell,suite2PoutputDir,fLookPostOnset = 3,treatment='pre'):
        self.treatment = treatment
        dfMapInput = dfMap[dfMap['treatment'].str.contains(treatment)]
        self.fs = np.unique(dfMapInput['frameRate'])[0]
        self.trigDelay = np.unique(dfMapInput['delayToEphusTrig_s'])[0]
        self.ISI = np.unique(dfMapInput['ISI'])[0]
        self.nPulsesPerTif = np.unique(dfMapInput['totalPulses'])[0]
        self.framesPerPulse = self.ISI*self.fs
        self.fLookPostOnset = fLookPostOnset
        
        #get list of freq/ampl in order of presentation and unique req/ampl vals
        self.freqOrder = np.hstack(dfMapInput['freq_Hz']).astype(int)
        self.amplOrder = np.hstack(dfMapInput['ampl_dB']).astype(int)
        self.freq = np.unique(self.freqOrder)
        self.ampl = np.unique(self.amplOrder)
        
        self.calcTRF(dfMapInput,dfTif,iscell,suite2PoutputDir,treatment)
    
    
    
    
    def calcTRF(self,dfMapInput,dfTif,iscell,suite2PoutputDir,treatment):
        """Isolates relevant frames from mapping .tifs and organizes 
        by ROI, amplitude, and frequency in TRF and for each ROI
        averages max response across traces w/ same ampl/freq (uTRF). Also
        calculates BF for each ROI."""

        mapSpksByTif = list(compress(suite2PvarByTif(dfTif,iscell,suite2PoutputDir),\
                                     np.logical_and(dfTif['stim'].str.contains('map',case=0),\
                                                    dfTif['treatment'].str.contains(treatment,case=0))))
        
        spkStack = np.hstack(list(map(lambda x: \
                                    x[:,int(self.fs*self.trigDelay):\
                                        int(self.fs*self.ISI*self.nPulsesPerTif+(self.fs*self.trigDelay))],\
                                        mapSpksByTif)))
        
        #create indices for pulse traces in spkStack
        maxFramePreOnset = min(np.hstack(dfMapInput['onset_s'])).astype(float)*self.fs
        maxFramePostOnset = self.framesPerPulse - max(np.hstack(dfMapInput['onset_s'])).astype(float)*self.fs
        pTraceStart = np.hstack(dfMapInput['onset_s']).astype(float)*self.fs-maxFramePreOnset
        pTraceEnd = np.hstack(dfMapInput['onset_s']).astype(float)*self.fs+maxFramePostOnset

        #calc indices along length of spkStack
        self.pTraceStarts = ((np.cumsum(np.ones(len(np.hstack(dfMapInput['onset_s'])))*self.framesPerPulse)-self.framesPerPulse)\
                             + pTraceStart).astype(int)
        self.pTraceEnds = ((np.cumsum(np.ones(len(np.hstack(dfMapInput['onset_s'])))*self.framesPerPulse)-self.framesPerPulse)\
                           + pTraceEnd).astype(int)
      
        #sort spkStack by pulses: spkTraces is (num pulses, num ROI, num frames in pulse)
        spkTraces = np.zeros((len(self.freqOrder),int(sum(iscell[:,0])),int(maxFramePreOnset + maxFramePostOnset)))
        for count,(start,stop) in enumerate(zip(self.pTraceStarts,self.pTraceEnds)):
            spkTraces[count,:,:] = spkStack[:,np.arange(start,stop)]
            
        #get BF and roi oct diff from BF
        #for each roi for each freq/ampl combination gets associated pulse traces and finds mean of max within a window
        #for each roi gets BF from 
        self.BF = []
        self.TRF = []
        self.uTRF = []
        
        for roi in np.arange(0,int(sum(iscell[:,0]))):
            
            self.uTRF.append(np.empty((len(self.freq),len(self.ampl)))) #average TRF for a given ROI
            self.TRF.append(np.empty((len(self.freq),len(self.ampl)),dtype=object)) #container for traces at a given freq/ampl
            
            for fq in self.freq:
                for amp in self.ampl:
                    
                    self.uTRF[roi][self.freq==fq,self.ampl==amp] = \
                        np.mean(np.max(spkTraces[np.logical_and(self.freqOrder==fq,self.amplOrder==amp),
                        roi,int(maxFramePreOnset):int(maxFramePreOnset)+self.fLookPostOnset],axis=1))

                    self.TRF[roi][(self.freq==fq).nonzero()[0][0],\
                                  (self.ampl==amp).nonzero()[0][0]] = \
                        spkTraces[np.logical_and(self.freqOrder==fq,self.amplOrder==amp),roi,:]
            
            self.BF.append(self.freq[self.uTRF[roi].mean(axis=1).argmax()])

    
    def plotTRFmeanAcrossCells(self):
        plt.figure()
        plt.imshow(np.flipud(np.asarray(self.uTRF).mean(axis=0).T)) 
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude (dB SPL)')
        plt.xticks(ticks=np.arange(0,len(self.freq),3),labels=self.freq[np.arange(0,len(self.freq),3)])
        plt.yticks(ticks=np.arange(0,len(self.ampl)),labels=self.ampl[::-1])
        
    def plotRoiTRF(self,ROIid):
        plt.figure()
        plt.imshow(np.flipud(self.uTRF[ROIid].T))
        #consider making this a function
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude (dB SPL)')
        plt.xticks(ticks=np.arange(0,len(self.freq),3),labels=self.freq[np.arange(0,len(self.freq),3)])
        plt.yticks(ticks=np.arange(0,len(self.ampl)),labels=self.ampl[::-1])