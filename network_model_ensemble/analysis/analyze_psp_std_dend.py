import sys
import os.path
import glob
import numpy as np

def compute_psp_std_windows(baseName, identifier, NREALIZATIONS):
    
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
        
        vmStdControl = np.std(goodTracesControl, axis=0)
        vmStdL1inact = np.std(goodTracesL1inact, axis=0)
        
        rIndex = pair[0].find('realization') + 11
        rStr = pair[0][rIndex:rIndex+2]
        if not summaryData.has_key(rStr):
            summaryData[rStr] = {}
        
        tOffset = 100.0
        dt = 0.025
        tStim = 200.0
        tBegin = 15.0
        windows = [35.0]
        vmStdsControl = []
        vmStdsL1inact = []
        for window in windows:
            beginBin = int((tStim+tBegin-tOffset)/dt+0.5)
            endBin = int((tStim+tBegin+window-tOffset)/dt+0.5)
            avgVmStdControl = np.mean(vmStdControl[beginBin:endBin])
            vmStdsControl.append(avgVmStdControl)
            avgVmStdL1inact = np.mean(vmStdL1inact[beginBin:endBin])
            vmStdsL1inact.append(avgVmStdL1inact)
            
            if not summaryData[rStr].has_key(window):
                summaryData[rStr][window] = {}
            summaryData[rStr][window]['L1inact'] = avgVmStdL1inact
            summaryData[rStr][window]['control'] = avgVmStdControl
    
    summaryName = os.path.join(baseName, identifier)
    
    with open(summaryName+'_vm_std_summary.csv', 'w') as summaryFile:
        header = 'realization\twindow\tcontrol\tL1 inactivated\n'
        summaryFile.write(header)
        rStrings = summaryData.keys()
        rStrings.sort()
        for rStr in rStrings:
            windows = summaryData[rStr].keys()
            windows.sort()
            for window in windows:
                line = rStr
                line += '\t'
                line += str(window)
                line += '\t'
                line += str(summaryData[rStr][window]['control'])
                line += '\t'
                line += str(summaryData[rStr][window]['L1inact'])
                line += '\n'
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
        compute_psp_std_windows(pathName, identifier, NREALIZATIONS)