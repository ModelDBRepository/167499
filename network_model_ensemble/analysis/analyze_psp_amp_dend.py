import sys
import os.path
import glob
import numpy as np

def compute_avg_psp_amplitudes(baseName, identifier, NREALIZATIONS):
    
    fnames = []
    scan_directory(baseName, fnames, identifier)
    realizationPairs = []
    for i in range(NREALIZATIONS):
        controlName, L1inactName = None, None
        for fname in fnames:
            rIndex = fname.find('realization') + 11
            NR = int(fname[rIndex:rIndex+2])
            if NR == i:
                if 'L1inact' in fname:
                    L1inactName = fname
                elif 'control' in fname:
                    controlName = fname
        realizationPairs.append((controlName, L1inactName))
    
    summaryData = {}
    for pair in realizationPairs:
        print 'Analyzing realization pair:'
        print 'control: %s' % pair[0]
        print 'L1 inactivated: %s' % pair[1]
        
        dataControl = np.loadtxt(pair[0], skiprows=2, unpack=True)
        dataL1inact = np.loadtxt(pair[1], skiprows=2, unpack=True)
        goodTracesControl = []
        goodTracesL1inact = []
        for i in range(1, len(dataControl)):
            goodTracesControl.append(dataControl[i])
            goodTracesL1inact.append(dataL1inact[i])
                
        rIndex = pair[0].find('realization') + 11
        rStr = pair[0][rIndex:rIndex+2]
        if not summaryData.has_key(rStr):
            summaryData[rStr] = {}
        
        tOffset = 100.0
        dt = 0.025
        tStim = 200.0
        tBegin = 10.0
        windows = [40.0]
        for window in windows:
            beginBin = int((tStim+tBegin-tOffset)/dt+0.5)
            endBin = int((tStim+tBegin+window-tOffset)/dt+0.5)
            if not summaryData[rStr].has_key(window):
                summaryData[rStr][window] = {}
            for i in range(len(goodTracesControl)):
                controlTrace = goodTracesControl[i]
                L1inactTrace = goodTracesL1inact[i]
                controlPSP = np.max(controlTrace[beginBin:endBin])
                L1inactPSP = np.max(L1inactTrace[beginBin:endBin])
                try:
                    summaryData[rStr][window]['L1inact'].append(L1inactPSP)
                    summaryData[rStr][window]['control'].append(controlPSP)
                except KeyError:
                    summaryData[rStr][window]['L1inact'] = [L1inactPSP]
                    summaryData[rStr][window]['control'] = [controlPSP]
        
    summaryName = os.path.join(baseName, identifier)
    
    with open(summaryName+'_psp_amplitudes_summary.csv', 'w') as summaryFile:
        header = 'realization\twindow\tcontrol avg\tcontrol std\tL1 inactivated avg\tL1 inactivated std\n'
        summaryFile.write(header)
        rStrings = summaryData.keys()
        rStrings.sort()
        for rStr in rStrings:
            windows = summaryData[rStr].keys()
            windows.sort()
            for window in windows:
                controlMean = np.mean(summaryData[rStr][window]['control'])
                controlStd = np.std(summaryData[rStr][window]['control'])
                L1inactMean = np.mean(summaryData[rStr][window]['L1inact'])
                L1inactStd = np.std(summaryData[rStr][window]['L1inact'])
                line = '%s\t%.1f\t%.5f\t%.5f\t%.5f\t%.5f\n' % (rStr, window, controlMean, controlStd, L1inactMean, L1inactStd)
                summaryFile.write(line)

def scan_directory(path, fnames, identifier):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory(fname, fnames, identifier)
        elif identifier in fname and fname.endswith('vm_dend_traces.csv'):
            fnames.append(fname)
        else:
            continue

if __name__=='__main__':
    if len(sys.argv) == 4:
        pathName = sys.argv[1]
        identifier = sys.argv[2]
        NREALIZATIONS = int(sys.argv[3])
        compute_avg_psp_amplitudes(pathName, identifier, NREALIZATIONS)