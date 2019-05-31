#!/bin/bash

modeldir=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)
source $modeldir"/pythonpath.sh"
cd $modeldir"/Figure4E_2"

cellModelName=$modeldir"/../cell_model/waters_l23_ar_dend_recording_dend_dense.param"

ongoingUp=$modeldir"/Figure4E_2/up_state_noNMDA_parameters_new_boutons.param"
ongoingDown=$modeldir"/Figure4E_2/down_state_noNMDA_parameters_new_boutons.param"
ongoingUpL1inact=$modeldir"/Figure4E_2/up_state_noNMDA_parameters_new_boutons_L1D1_inactivated.param"
ongoingDownL1inact=$modeldir"/Figure4E_2/down_state_noNMDA_parameters_new_boutons_L1D1_inactivated.param"

python $modeldir'/evoked_activity_replay_dend_rec.py' $modeldir'/Figure4E_2/results/realization00/control' $cellModelName $ongoingUp $ongoingDown 'evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23746.param' 'evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23746.param'

python $modeldir'/evoked_activity_replay_dend_rec.py' $modeldir'/Figure4E_2/results/realization00/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact 'evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23746.param' 'evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23746.param'

echo "*********************"
echo "Analysis: PSP amplitude"
for i in {0..34}
do
    realization=$(printf "%03d" $i)
    python $modeldir"/analysis/analyze_psp_amp_dend.py" $modeldir"/Figure4E_2/results" "apical_recording_sites_ID_"$realization 1
done
echo "Analysis: PSP SD"
for i in {0..34}
do
    realization=$(printf "%03d" $i)
    python $modeldir"/analysis/analyze_psp_std_dend.py" $modeldir"/Figure4E_2/results" "apical_recording_sites_ID_"$realization 1
done
