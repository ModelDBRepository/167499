Model of a cortical layer (L) 2 pyramidal neuron embedded in an anatomically realistic network of
two barrel columns in rat vibrissal cortex. This model is used to investigate the effects of
spatially and temporally specific inhibition from L1 inhibitory interneurons on the sensory-evoked
subthreshold responses of the L2 pyramidal neuron, and can be used to create simulation results
underlying Figures 3D, 4B, 4C and 4E from
Egger, Schmitt et al.:  Robustness of sensory-evoked excitation is increased by inhibitory inputs
                        to distal apical tuft dendrites (PNAS 2015).

Contact: Robert Egger (robert.egger@tuebingen.mpg.de)


REQUIREMENTS:
    Required and tested in a Linux environment with:
    NEURON 7.2 (http://www.neuron.yale.edu/ftp/neuron/versions/v7.2/)
        NEURON needs to be installed such that it can be imported as a python module;
        see http://www.davison.webfactional.com/notes/installation-neuron-python/
    python 2.7 (http://www.python.org/getit/)
    numpy 1.6 (http://www.scipy.org/install.html)
    sumatra 0.3 (http://neuralensemble.org/sumatra/)
    
    Optional:
    matplotlib 1.1 (http://www.scipy.org/install.html)
    
    Older versions of numpy/matplotlib are very likely to work, but have not been tested.
    NEURON < 7.2 does NOT work, due to differences in the python/hoc interface. 7.3/7.4 are likely to work,
    but has not been tested.

SETUP:
    Make sure NEURON is installed with full python support, i.e., it can be imported as a python module.
    If you are not sure whether this is the case in your installation, start an interactive python
    session, type
    import neuron
    and hit enter. If installed correctly, information about the NEURON version installed on your system
    will be printed.
    If you get an error message, you may need to re-compile neuron from the source code. A detailed
    installation guide can be found at
    http://www.davison.webfactional.com/notes/installation-neuron-python/
    
    Compile the NEURON model mechanisms: change into the "mechanisms" folder and run nrnivmodl
    
    Adjust installation paths:
    1. Adjust the installation paths in the cell model .param files:
        - open "cell_model/waters_l23_ar_dend.param", find the "neuron" "filename" parameter value and
        change the path to the absolute path of the neuron reconstruction .hoc file on your system.
        - open "cell_model/waters_l23_ar_dend_recording_dend_dense.param", find the "neuron" "filename"
        parameter value and change the path to the absolute path of the neuron reconstruction .hoc
        file on your system.
        - open "cell_model/waters_l23_ar_dend_recording_dend_dense.param", find the "recordingSites"
        parameter value and change the path to the absolute path of the file
        "cell_model/apical_recording_sites.landmarkAscii"
    2. Adjust the PYTHONPATH variable:
        - open the file "network_model_ensemble/pythonpath.sh" and change the path containing
        "model_publication/lib" to the absolute installation path of this model on your system.
        E.g., if you unpacked the model in the folder "/home/user/model", the line should read:
        export PYTHONPATH=/home/user/model/lib:$PYTHONPATH

RUN:
    To produce results underlying the analyses in Figures 3D, 4B, 4C and 4E, open the network_model_ensemble
    folder. Next, run the respective bash script (e.g., run_Figure3D.sh).
    Upon execution of the script, a folder called "results" is automatically created in the appropriate
    figure directory (e.g., Figure3D/results), and the analysis is automatically started and saved into
    .csv files in the "results" folder.
    The individual membrane potential traces and information about location and activation times of
    all synapses for all simulations can be found in the subfolders of the "results" folder.
    The analysis scripts can be found in the "analysis" directory.
    
    WARNING: running the scripts for Figure3D, Figure4B, Figure4C_1, Figure4C_2 and Figure4C_3 can take
    a long time, since between 50,000 and 100,000 simulation trials are started in these cases. If you
    have access to a multicore system, we recommend breaking up the simulations into smaller batches that
    can be run in parallel, depending on the number of cores available. In this case special care needs
    to be taken of the simulations Figure4C_1, Figure4C_2 and Figure4C_3. Here, the second half of the
    simulations (separated by a blank line in the shell script) can only be run after the first half has
    completed, because functional connectivity configurations are loaded from the results of the first half.
    For clarity, we have provided two example scripts showing a configuration that would allow splitting up
    the simulations for Figure4C_2 on a system with two cores (run_Figure4C_2_parallel1.sh and
    run_Figure4C_2_parallel2.sh).

Neuron morphology and synapses:
    The L2 pyramidal neuron morphology .hoc file, the synapse distribution files, as well as screenshots
    of the subcellular distribution of synapses from all presynaptic cell types on the L2 pyramidal
    neuron can be found in the "NeuroNet" directory and its subfolders.

Acknowledgements:
    We thank Dr. Christiaan de Kock (VU University Amsterdam) for providing the L2 pyramidal
    neuron reconstruction.
