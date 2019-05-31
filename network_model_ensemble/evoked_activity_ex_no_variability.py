'''
L2 neuron model
evoked activity with/without L1 PW neurons
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

def evoked_activity_replay_synapses(simName, cellName, ongoingUpName, ongoingDownName, evokedUpParamName, evokedDownParamName):
    '''
    pre-stimulus ongoing activity
    and evoked activity
    L1 inactivation simulated by replaying exact
    synapse activation times, only without L1D1 evoked
    '''
    neuronParameters = scp.build_parameters(cellName)
    upParameters = scp.build_parameters(ongoingUpName)
    downParameters = scp.build_parameters(ongoingDownName)
    evokedUpNWParameters = scp.build_parameters(evokedUpParamName)
    evokedDownNWParameters = scp.build_parameters(evokedDownParamName)
    scp.load_NMODL_parameters(neuronParameters)
    scp.load_NMODL_parameters(upParameters)
    scp.load_NMODL_parameters(downParameters)
    scp.load_NMODL_parameters(evokedUpNWParameters)
    cellParam = neuronParameters.neuron
    paramUp = upParameters.network
    paramDown = downParameters.network
    paramEvokedUp = evokedUpNWParameters.network
    paramEvokedDown = evokedDownNWParameters.network
    
    cell = scp.create_cell(cellParam, scaleFunc=dendriteScalingUniform)
            
    uniqueID = str(os.getpid())
    dirName = simName
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
    # control runs: create trials with random activity;
    # save exact synapse activation times for replay;
    # save L1D1 activation times in separate file for analysis
    #===========================================================================
    if 'control' in simName:
        nRun = 0
        while nRun < nSweeps:
            # 50% down states, 50% up states
            if nRun < 0.5*nSweeps:
                synParameters = paramDown
                synParametersEvoked = paramEvokedDown
            else:
                synParameters = paramUp
                synParametersEvoked = paramEvokedUp
            
            ongoingNW = scp.NetworkMapper(cell, synParameters, neuronParameters.sim)
            ongoingNW.create_network()
            evokedNW = scp.NetworkMapper(cell, synParametersEvoked, neuronParameters.sim)
            evokedNW.create_saved_network()
            
            synTypes , synTypeL1D1 = [], []
            for synType in cell.synapses.keys():
                # replay all synapses except for L1D1 evoked/ongoing
                if 'L1D1' in synType:
                    synTypeL1D1.append(synType)
                else:
                    synTypes.append(synType)
            synTypes.sort()
            
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
            #===================================================================
            # discard trials above threshold!
            #===================================================================
            begin = int((tStim+15.0)/dt+0.5)
            end = int((tStim+50.0)/dt+0.5)
            maxV = np.max(vmSoma[begin:end])
            if maxV < spikeThresh:
                vTraces.append(np.array(vmSoma[offsetBin:])), tTraces.append(np.array(t[offsetBin:]))
                
                print 'writing simulation results'
                fname = 'simulation_'
                fname += uniqueID
                fname += '_run%04d' % nRun
                if nRun < 0.5*nSweeps:
                    fname += '_down_state'
                else:
                    fname += '_up_state'
                synName = dirName + '/' + fname + '_synapses.csv'
                synNameL1D1 = dirName + '/' + fname + '_synapsesL1D1.csv'
                print 'computing active synapse properties'
                sca.compute_synapse_distances_times(synName, cell, t, synTypes)
                sca.compute_synapse_distances_times(synNameL1D1, cell, t, synTypeL1D1)
                
                nRun += 1
            else:
                print 'Trial above threshold; running new trial'
            
            cell.re_init_cell()
            cell.remove_synapses('all')
            ongoingNW.re_init_network()
            evokedNW.re_init_network()
    
            print '-------------------------------'
    
    #===========================================================================
    # Load synapse activation times for all trials of corresponding
    # network realization. They should NOT have the L1D1 activation times
    #===========================================================================
    elif 'L1inact' in simName:
        tmpStr = simName
        controlBasePath = tmpStr.replace('L1inact', 'control')
        synInfoNames = []
        scan_directory(controlBasePath, synInfoNames, '_synapses.csv')
        
        nRun = 0
        while nRun < nSweeps:
            # 50% down states, 50% up states
            if nRun < 0.5*nSweeps:
                synParameters = paramDown
                synParametersEvoked = paramEvokedDown
            else:
                synParameters = paramUp
                synParametersEvoked = paramEvokedUp
            
            synInfoName = ''
            nRunStr = 'run%04d' % nRun
            for name in synInfoNames:
                if nRunStr in name:
                    synInfoName = name
                    break
            
            synParameters.update(synParametersEvoked)
            for synType in synParameters.keys():
                synParameters[synType].synapses.releaseProb = 1.0
            
            print 'Replaying network activity from file %s' % synInfoName
            replayNW = scp.NetworkMapper(cell, synParameters)
            replayNW.reconnect_saved_synapses(synInfoName)
            
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
            #===================================================================
            # no max V check here!!! we want the trace to be computed in any case!
            #===================================================================
            vTraces.append(np.array(vmSoma[offsetBin:])), tTraces.append(np.array(t[offsetBin:]))
            nRun += 1
            
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
    upParameters.save(dirName+'/'+uniqueID+'_network_model_upstate.param')
    downParameters.save(dirName+'/'+uniqueID+'_network_model_downstate.param')
    evokedUpNWParameters.save(dirName+'/'+uniqueID+'_network_model_evoked_upstate.param')
    evokedDownNWParameters.save(dirName+'/'+uniqueID+'_network_model_evoked_downstate.param')
    
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
        ongoingUpName = sys.argv[3]
        ongoingDownName = sys.argv[4]
        evokedUpName = sys.argv[5]
        evokedDownName = sys.argv[6]
        evoked_activity_replay_synapses(name, cellName, ongoingUpName, ongoingDownName, evokedUpName, evokedDownName)
    else:
        print 'Error! Number of arguments is %d; should be 6!' % (len(sys.argv)-1)
    
    
    
    