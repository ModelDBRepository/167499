import sys
import os.path
import glob
import numpy as np

def compute_ongoing_vm_parameters(baseName, summaryName):
    suffix = 'all_traces.csv'
    fnames = []
    scan_directory(baseName, fnames, suffix)
    
    summaryData = {}
    for fname in fnames:
        data = np.loadtxt(fname, skiprows=2, unpack=True)
        t = data[0]
        downUpBorder = len(data)//2 + 1
        
        rIndex = fname.find('realization') + 11
        rStr = fname[rIndex:rIndex+2]
        if not summaryData.has_key(rStr):
            summaryData[rStr] = {}
        
        tOffset = 100.0
        dt = 0.025
        tStim = 200.0
        tBegin = -50.0
        windows = [50.0]
        for window in windows:
            beginBin = int((tStim+tBegin-tOffset)/dt+0.5)
            endBin = int((tStim+tBegin+window-tOffset)/dt+0.5)
            vmUp = []
            stdUp = []
            vmDown = []
            stdDown = []
            for i in range(1, downUpBorder):
                avgVm = np.mean(data[i][beginBin:endBin])
                stdVm = np.std(data[i][beginBin:endBin])
                vmDown.append(avgVm)
                stdDown.append(stdVm)
            for i in range(downUpBorder, len(data)):
                avgVm = np.mean(data[i][beginBin:endBin])
                stdVm = np.std(data[i][beginBin:endBin])
                vmUp.append(avgVm)
                stdUp.append(stdVm)
            avgVmUp = np.mean(vmUp)
            avgVmDown = np.mean(vmDown)
            avgStdUp = np.mean(stdUp)
            avgStdDown = np.mean(stdDown)
            
            if not summaryData[rStr].has_key(window):
                summaryData[rStr][window] = {}
                summaryData[rStr][window]['up'] = {}
                summaryData[rStr][window]['down'] = {}
            summaryData[rStr][window]['up']['avg'] = avgVmUp
            summaryData[rStr][window]['up']['std'] = avgStdUp
            summaryData[rStr][window]['down']['avg'] = avgVmDown
            summaryData[rStr][window]['down']['std'] = avgStdDown
    
    with open(summaryName, 'w') as summaryFile:
        header = 'realization\twindow\tstate\tVm avg\tVm STD\n'
        summaryFile.write(header)
        rStrings = summaryData.keys()
        rStrings.sort()
        for rStr in rStrings:
            windows = summaryData[rStr].keys()
            windows.sort()
            for window in windows:
                states = summaryData[rStr][window].keys()
                states.sort()
                for state in states:
                    line = rStr
                    line += '\t'
                    line += str(window)
                    line += '\t'
                    line += state
                    line += '\t'
                    line += str(summaryData[rStr][window][state]['avg'])
                    line += '\t'
                    line += str(summaryData[rStr][window][state]['std'])
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
        compute_ongoing_vm_parameters(pathName, summaryName)
