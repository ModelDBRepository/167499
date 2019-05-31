#!/bin/bash

modeldir=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)
source $modeldir"/pythonpath.sh"
cd $modeldir"/Figure4C_4"

cellModelName=$modeldir"/../cell_model/waters_l23_ar_dend.param"

ongoingUp=$modeldir"/Figure4C_4/up_state_parameters_new_boutons.param"
replayDir=$modeldir"/Figure4C_4/replay"

python $modeldir'/evoked_activity_inh_conductances.py' $modeldir'/Figure4C_4/results/control' $cellModelName 0.0 $ongoingUp $modeldir'/Figure4C_4/control_up_state_functional_map_20130326-1528_23746.param' $replayDir
python $modeldir'/evoked_activity_inh_conductances.py' $modeldir'/Figure4C_4/results/control' $cellModelName 0.4 $ongoingUp $modeldir'/Figure4C_4/control_up_state_functional_map_20130326-1528_23746.param' $replayDir
python $modeldir'/evoked_activity_inh_conductances.py' $modeldir'/Figure4C_4/results/control' $cellModelName 0.8 $ongoingUp $modeldir'/Figure4C_4/control_up_state_functional_map_20130326-1528_23746.param' $replayDir
python $modeldir'/evoked_activity_inh_conductances.py' $modeldir'/Figure4C_4/results/control' $cellModelName 1.6 $ongoingUp $modeldir'/Figure4C_4/control_up_state_functional_map_20130326-1528_23746.param' $replayDir
python $modeldir'/evoked_activity_inh_conductances.py' $modeldir'/Figure4C_4/results/control' $cellModelName 3.0 $ongoingUp $modeldir'/Figure4C_4/control_up_state_functional_map_20130326-1528_23746.param' $replayDir
python $modeldir'/evoked_activity_inh_conductances.py' $modeldir'/Figure4C_4/results/control' $cellModelName 6.0 $ongoingUp $modeldir'/Figure4C_4/control_up_state_functional_map_20130326-1528_23746.param' $replayDir
python $modeldir'/evoked_activity_inh_conductances.py' $modeldir'/Figure4C_4/results/control' $cellModelName 10.0 $ongoingUp $modeldir'/Figure4C_4/control_up_state_functional_map_20130326-1528_23746.param' $replayDir
python $modeldir'/evoked_activity_inh_conductances.py' $modeldir'/Figure4C_4/results/control' $cellModelName 15.0 $ongoingUp $modeldir'/Figure4C_4/control_up_state_functional_map_20130326-1528_23746.param' $replayDir

echo "*********************"
echo "Analysis: PSP SD"
python $modeldir"/analysis/analyze_psp_std.py" $modeldir"/Figure4C_4/results" $modeldir"/Figure4C_4/results/summary_psp_std.csv"
