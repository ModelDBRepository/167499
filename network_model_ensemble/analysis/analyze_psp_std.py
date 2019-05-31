import sys
import os.path
import glob
import numpy as np

def compute_psp_std_windows(baseName, summaryName, NREALIZATIONS, case):
    eps = 1e-6
    TRACESPERSTATE = 1000
    spikeThreshold = -38.0 + eps
    
    suffix = 'all_traces.csv'
    fnames = []
    scan_directory(baseName, fnames, suffix)
    realizationPairs = []
    for i in range(NREALIZATIONS):
        controlName, L1inactName = None, None
        for fname in fnames:
            rIndex = fname.find('realization') + 11
            NR = int(fname[rIndex:rIndex+2])
            if NR == i:
                if case == 1:
                    if 'L1inact' in fname:
                        L1inactName = fname
                    elif 'control' in fname:
                        controlName = fname
                elif case == 2:
                    if 'control1' in fname:
                        controlName = fname
                    elif 'control2' in fname:
                        L1inactName = fname
        realizationPairs.append((controlName, L1inactName))
    
    summaryData = {}
    for pair in realizationPairs:
        print 'Analyzing realization pair:'
        if case == 1:
            print 'control: %s' % pair[0]
            print 'L1 inactivated: %s' % pair[1]
        elif case == 2:
            print 'control 1: %s' % pair[0]
            print 'control 2: %s' % pair[1]
        
        dataControl = np.loadtxt(pair[0], skiprows=2, unpack=True)
        dataL1inact = np.loadtxt(pair[1], skiprows=2, unpack=True)
        goodTracesControl = []
        goodTracesL1inact = []
        for i in range(1, len(dataControl)):
            maxV1 = np.max(dataControl[i])
            maxV2 = np.max(dataL1inact[i])
            if maxV1 > spikeThreshold or maxV2 > spikeThreshold:
                print '\tabove threshold!'
                print '\tmaxV1 = %.2f' % maxV1
                print '\tmaxV2 = %.2f' % maxV2
                continue
            goodTracesControl.append(dataControl[i])
            goodTracesL1inact.append(dataL1inact[i])
        print 'Discarded %d traces above spike threshold' % (len(dataControl)-1-len(goodTracesControl))
        
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
            if case == 1:
                summaryData[rStr][window]['L1inact'] = avgVmStdL1inact
                summaryData[rStr][window]['control'] = avgVmStdControl
            elif case == 2:
                summaryData[rStr][window]['control2'] = avgVmStdL1inact
                summaryData[rStr][window]['control1'] = avgVmStdControl
    
    with open(summaryName, 'w') as summaryFile:
        if case == 1:
            header = 'realization\twindow\tcontrol\tL1 inactivated\n'
        elif case == 2:
            header = 'realization\twindow\tcontrol1\tcontrol2\n'
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
                if case == 1:
                    line += '\t'
                    line += str(summaryData[rStr][window]['control'])
                    line += '\t'
                    line += str(summaryData[rStr][window]['L1inact'])
                elif case == 2:
                    line += '\t'
                    line += str(summaryData[rStr][window]['control1'])
                    line += '\t'
                    line += str(summaryData[rStr][window]['control2'])
                line += '\n'
                summaryFile.write(line)

def compute_psp_std_windows_gGABA(baseName, summaryName):
    eps = 1e-6
    GGABAList = ['0.0','0.4','0.8','1.6','3.0','6.0','10.0','15.0']
    spikeThreshold = -38.0 + eps
    
    suffix  = 'all_traces.csv'
    fnames = []
    scan_directory(baseName, fnames, suffix)
    realizationPairs = []
    for i in range(len(GGABAList)):
        controlName, L1inactName = None, None
        for fname in fnames:
            rIndex = fname.find('gGABA') + 6
            endIndex = fname.find('/', rIndex)
            rStr = fname[rIndex:endIndex]
            if rStr == GGABAList[i]:
                if 'L1inact' in fname:
                    L1inactName = fname
                elif 'control' in fname:
                    controlName = fname
        realizationPairs.append((controlName, L1inactName))
    
    summaryData = {}
    for pair in realizationPairs:
        print 'Analyzing realization:'
        print pair[0]
        
        dataControl = np.loadtxt(pair[0], skiprows=2, unpack=True)
        goodTracesControl = []
        for i in range(1, len(dataControl)):
            maxV = np.max(dataControl[i])
            if maxV > spikeThreshold:
                print '\tabove threshold!'
                print '\tmaxV = %.2f' % maxV
                continue
            goodTracesControl.append(dataControl[i])
        print 'Discarded %d traces above spike threshold' % (len(dataControl)-1-len(goodTracesControl))
        
        vmStdControl = np.std(goodTracesControl, axis=0)
        rIndex = pair[0].find('gGABA') + 6
        endIndex = pair[0].find('/', rIndex)
        rStr = pair[0][rIndex:endIndex]
        if not summaryData.has_key(rStr):
            summaryData[rStr] = {}
        
        tOffset = 100.0
        dt = 0.025
        tStim = 200.0
        tBegin = 15.0
        windows = [35.0]
        vmStdsControl = []
        for window in windows:
            beginBin = int((tStim+tBegin-tOffset)/dt+0.5)
            endBin = int((tStim+tBegin+window-tOffset)/dt+0.5)
            avgVmStdControl = np.mean(vmStdControl[beginBin:endBin])
            vmStdsControl.append(avgVmStdControl)
            
            if not summaryData[rStr].has_key(window):
                summaryData[rStr][window] = {}
            summaryData[rStr][window]['control'] = avgVmStdControl
    
    with open(summaryName, 'w') as summaryFile:
        header = 'gGABA\twindow\tVm STD\n'
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
                line += '\n'
                summaryFile.write(line)

def scan_directory(path, fnames, suffix):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory(fname, fnames, suffix)
        elif fname.endswith(suffix):
            fnames.append(fname)
        else:
            continue

if __name__=='__main__':
    if len(sys.argv) == 3:
        pathName = sys.argv[1]
        summaryName = sys.argv[2]
        compute_psp_std_windows_gGABA(pathName, summaryName)
    if len(sys.argv) == 5:
        pathName = sys.argv[1]
        summaryName = sys.argv[2]
        NREALIZATIONS = int(sys.argv[3])
        case = int(sys.argv[4])
        compute_psp_std_windows(pathName, summaryName, NREALIZATIONS, case)