import sys
import os.path
import glob
import numpy as np

def compute_avg_psp_amplitudes(baseName, summaryName, NREALIZATIONS, case):
    eps = 1e-6
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
        
        rIndex = fname.find('realization') + 11
        rStr = fname[rIndex:rIndex+2]
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
            PSPsControl = []
            PSPsL1inact = []
            for i in range(len(goodTracesControl)):
                PSPControl = np.max(goodTracesControl[i][beginBin:endBin])
                PSPsControl.append(PSPControl)
                PSPL1inact = np.max(goodTracesL1inact[i][beginBin:endBin])
                PSPsL1inact.append(L1inact)
            avgPSPControl = np.mean(PSPsControl)
            avgPSPL1inact = np.mean(PSPsL1inact)
            
            if not summaryData[rStr].has_key(window):
                summaryData[rStr][window] = {}
            if case == 1:
                summaryData[rStr][window]['L1inact'] = avgPSPL1inact
                summaryData[rStr][window]['control'] = avgPSPControl
            elif case == 2:
                summaryData[rStr][window]['control2'] = avgPSPL1inact
                summaryData[rStr][window]['control1'] = avgPSPControl
    
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

def compute_avg_psp_amplitudes_unpaired(baseName, summaryName):
    eps = 1e-6
    spikeThreshold = -38.0 + eps
    
    suffix = 'all_traces.csv'
    fnames = []
    scan_directory(baseName, fnames, suffix)
    
    summaryData = {}
    for fname in fnames:
        print 'Analyzing realization:'
        print fname
        
        data = np.loadtxt(fname, skiprows=2, unpack=True)
        goodTraces = []
        for i in range(1, len(data)):
            maxV = np.max(data[i])
            if maxV > spikeThreshold:
                print '\tabove threshold!'
                print '\tmaxV = %.2f' % maxV
                continue
            goodTraces.append(data[i])
        print 'Discarded %d traces above spike threshold' % (len(data)-1-len(goodTraces))
        
        rIndex = fname.find('realization') + 11
        rStr = fname[rIndex:rIndex+2]
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
            PSPs = []
            for i in range(len(goodTraces)):
                PSP = np.max(goodTraces[i][beginBin:endBin])
                PSPs.append(PSP)
            avgPSP = np.mean(PSPs)
            summaryData[rStr][window] = avgPSP
    
    with open(summaryName, 'w') as summaryFile:
        header = 'realization\twindow\tPSP\n'
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
                line += str(summaryData[rStr][window])
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
        compute_avg_psp_amplitudes_unpaired(pathName, summaryName)
    if len(sys.argv) == 5:
        pathName = sys.argv[1]
        summaryName = sys.argv[2]
        NREALIZATIONS = int(sys.argv[3])
        case = int(sys.argv[4])
        compute_avg_psp_amplitudes(pathName, summaryName, NREALIZATIONS, case)