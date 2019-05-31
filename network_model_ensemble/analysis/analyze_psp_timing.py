import sys
import os.path
import glob
import numpy as np

def compute_avg_psp_times_unpaired(baseName, summaryName):
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
        t = data[0]
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
            PSPTimes = []
            for i in range(len(goodTraces)):
                PSPBin = np.argmax(goodTraces[i][beginBin:endBin]) + beginBin
                PSPTimes.append(t[PSPBin]-tStim)
            avgPSPTime = np.mean(PSPTimes)
            stdPSPTime = np.std(PSPTimes)
            
            if not summaryData[rStr].has_key(window):
                summaryData[rStr][window] = {}
            summaryData[rStr][window]['avg'] = avgPSPTime
            summaryData[rStr][window]['std'] = stdPSPTime
    
    with open(summaryName, 'w') as summaryFile:
        header = 'realization\twindow\tavg peak latency\tstd peak latency\n'
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
                line += str(summaryData[rStr][window]['avg'])
                line += '\t'
                line += str(summaryData[rStr][window]['std'])
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
        compute_avg_psp_times_unpaired(pathName, summaryName)