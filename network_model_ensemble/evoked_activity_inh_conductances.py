'''
L2 neuron model
evoked activity with/without L1 PW neurons
Analysis of sensitivity to inhibitory conductance strength
No variation in functional connectivity
of excitatory synapses between conditions
(i.e., activation times of excitatory synapses
are replayed between conditions)

written by Robert Egger (robert.egger@tuebingen.mpg.de)
(c) 2013-2015 Max Planck Society
'''

import sys
import time
import os, os.path
import glob
import neuron
import single_cell_parser as scp
import single_cell_analyzer as sca
import numpy as np
hasMatplotlib = True
try:
    import matplotlib.pyplot as plt
except ImportError:
    hasMatplotlib = False
h = neuron.h

def evoked_activity_replay_synapses_new(simName, cellName, gGABA, ongoingName, evokedParamName, replayFolder):
    '''
    pre-stimulus ongoing activity
    and evoked activity replayed from
    control trial of all variable case
    '''
    neuronParameters = scp.build_parameters(cellName)
    ongoingParameters = scp.build_parameters(ongoingName)
    evokedNWParameters = scp.build_parameters(evokedParamName)
    scp.load_NMODL_parameters(neuronParameters)
    scp.load_NMODL_parameters(ongoingParameters)
    scp.load_NMODL_parameters(evokedNWParameters)
    cellParam = neuronParameters.neuron
    paramOngoing = ongoingParameters.network
    paramEvoked = evokedNWParameters.network
    
    cell = scp.create_cell(cellParam, scaleFunc=dendriteScalingUniform)
            
    uniqueID = str(os.getpid())
    dirName = simName
    dirName += '_gGABA_'
    dirName += '%.1f' % gGABA
    if not simName.endswith('/'):
        dirName += '/'
    dirName += time.strftime('%Y%m%d-%H%M')
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    
    vTraces = []
    tTraces = []
    
    nSweeps = 2000
    #nSweeps = 2
    tOffset = 100.0 # avoid numerical transients
    tStim = 200.0
    tStop = 250.0
    neuronParameters.sim.tStop = tStop
    dt = neuronParameters.sim.dt
    offsetBin = int(tOffset/dt + 0.5)
    
    spikeThresh = -38.0 # Petersen, AS
    tStim = 200.0
    
    #===========================================================================
    # Load synapse activation times for all trials of corresponding
    # network realization from 'replayFolder'.
    # Replay all synapses and only vary L1D1 evoked gGABA
    #===========================================================================
    synInfoNames = []
    scan_directory(replayFolder, synInfoNames, '_synapses.csv')
    
    nRun = 0
    while nRun < nSweeps:
        synParameters = paramOngoing
        synParametersEvoked = paramEvoked
        
        synInfoName = ''
        nRunStr = 'run%04d' % nRun
        for name in synInfoNames:
            if nRunStr in name:
                synInfoName = name
                break
        
        synParameters.update(synParametersEvoked)
        for synType in synParameters.keys():
            synParameters[synType].synapses.releaseProb = 1.0
        synParameters['L1D1'].synapses.receptors.gaba_syn.weight = gGABA
        
        print 'Replaying network activity from file %s' % synInfoName
        replayNW = scp.NetworkMapper(cell, synParameters)
        replayNW.reconnect_saved_synapses(synInfoName)
        if gGABA == 0.0:
            if cell.synapses.has_key('L1D1'):
                for syn in cell.synapses['L1D1']:
                    syn.disconnect_hoc_synapse()
                cell.remove_synapses('L1D1')
            if cell.synapses.has_key('L1D1_ongoing'):
                for syn in cell.synapses['L1D1_ongoing']:
                    syn.disconnect_hoc_synapse()
                cell.remove_synapses('L1D1_ongoing')
        
        if cell.synapses.has_key('L1D1'):
            for i in range(len(cell.synapses['L1D1'])):
                syn = cell.synapses['L1D1'][i]
                if syn.is_active():
                    L1strength = syn.netcons[0].weight[0]
                    print 'L1 synapse weight = %.1fnS' % L1strength
        
        print 'Testing evoked response properties run %d of %d' % (nRun+1, nSweeps)
        tVec = h.Vector()
        tVec.record(h._ref_t)
        startTime = time.time()
        scp.init_neuron_run(neuronParameters.sim)
        stopTime = time.time()
        simdt = stopTime - startTime
        print 'NEURON runtime: %.2f s' % simdt
        
        vmSoma = np.array(cell.soma.recVList[0])
        t = np.array(tVec)
        begin = int((tStim+15.0)/dt+0.5)
        end = int((tStim+50.0)/dt+0.5)
        maxV = np.max(vmSoma[begin:end])
        vTraces.append(np.array(vmSoma[offsetBin:])), tTraces.append(np.array(t[offsetBin:]))
        nRun += 1
        #===================================================================
        # no max V check here!!! we want the trace to be computed in any case!
        #===================================================================
        
        cell.re_init_cell()
        cell.remove_synapses('all')
        replayNW.re_init_network()

        print '-------------------------------'
    
    vTraces = np.array(vTraces)
    print 'computing Vm STD and histogram'
    vStd = np.std(vTraces, axis=0)
    peakWindow, avgPeak = sca.compute_mean_psp_amplitude(vTraces, tStim=200.0-tOffset, dt=neuronParameters.sim.dt)
    windows, avgVmStd = sca.compute_vm_std_windows(vStd, tStim=200.0-tOffset, dt=neuronParameters.sim.dt)
    hist, bins = sca.compute_vm_histogram(vTraces)
    scp.write_all_traces(dirName+'/'+uniqueID+'_vm_all_traces.csv', t[offsetBin:], vTraces)
    scp.write_sim_results(dirName+'/'+uniqueID+'_vm_std.csv', t[offsetBin:], vStd)
    scp.write_sim_results(dirName+'/'+uniqueID+'_vm_avg_psp.csv', peakWindow, avgPeak)
    scp.write_sim_results(dirName+'/'+uniqueID+'_vm_std_windows.csv', windows, avgVmStd)
    scp.write_sim_results(dirName+'/'+uniqueID+'_vm_hist.csv', hist, bins[:-1])
    
    print 'writing simulation parameter files'
    neuronParameters.save(dirName+'/'+uniqueID+'_neuron_model.param')
    ongoingParameters.save(dirName+'/'+uniqueID+'_network_model_upstate.param')
    evokedNWParameters.save(dirName+'/'+uniqueID+'_network_model_evoked_upstate.param')
    
    if hasMatplotlib:
        ax = []
        plt.figure()
        for i in range(nSweeps):
            ax.append(plt.plot(tTraces[i], vTraces[i], 'k'))
        plt.xlabel('t [ms]')
        plt.ylabel('Vm [mV]')
        plt.savefig(dirName+'/'+uniqueID+'_all_traces.pdf')
        plt.figure()
        plt.plot(tTraces[0], vStd, 'k')
        plt.xlabel('t [ms]')
        plt.ylabel('Vm STD [mV]')
        plt.savefig(dirName+'/'+uniqueID+'_vm_std.pdf')

def scan_directory(path, fnames, suffix):
    for fname in glob.glob(os.path.join(path, '*')):
        if os.path.isdir(fname):
            scan_directory(fname, fnames, suffix)
        elif fname.endswith(suffix):
            fnames.append(fname)
        else:
            continue

def dendriteScalingUniform(cell):
    dendScale = 1/1.2
    for sec in cell.sections:
        if sec.label == 'Dendrite' or sec.label == 'ApicalDendrite':
            dummy = h.pt3dclear(sec=sec)
            for i in range(sec.nrOfPts):
                x, y, z = sec.pts[i]
                sec.diamList[i] = sec.diamList[i]*dendScale
                d = sec.diamList[i]
                dummy = h.pt3dadd(x, y, z, d, sec=sec)

if __name__ == '__main__':
    if len(sys.argv) == 7:
        name = sys.argv[1]
        cellName = sys.argv[2]
        gGABA = float(sys.argv[3])
        ongoingName = sys.argv[4]
        evokedName = sys.argv[5]
        replayFolderName = sys.argv[6]
        evoked_activity_replay_synapses_new(name, cellName, gGABA, ongoingName, evokedName, replayFolderName)
    else:
        print 'Error! Number of arguments is %d; should be 5!' % (len(sys.argv)-1)
    
    
    
    