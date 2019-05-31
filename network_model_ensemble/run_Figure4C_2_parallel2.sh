#!/bin/bash

modeldir=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)
source $modeldir"/pythonpath.sh"
cd $modeldir"/Figure4C_2"

cellModelName=$modeldir"/../cell_model/waters_l23_ar_dend.param"

ongoingUp=$modeldir"/Figure4C_2/up_state_parameters_new_boutons.param"
ongoingDown=$modeldir"/Figure4C_2/down_state_parameters_new_boutons.param"
ongoingUpL1inact=$modeldir"/Figure4C_2/up_state_parameters_new_boutons_L1D1_inactivated.param"
ongoingDownL1inact=$modeldir"/Figure4C_2/down_state_parameters_new_boutons_L1D1_inactivated.param"

python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization25/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23650.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23650.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization26/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23654.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23654.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization27/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23658.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23658.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization28/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23662.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23662.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization29/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23666.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23666.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization30/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23670.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23670.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization31/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23674.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23674.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization32/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23678.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23678.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization33/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23682.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23682.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization34/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23686.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23686.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization35/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23690.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23690.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization36/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23694.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23694.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization37/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23698.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23698.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization38/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23702.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23702.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization39/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23706.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23706.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization40/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23710.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23710.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization41/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23714.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23714.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization42/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23718.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23718.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization43/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23722.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23722.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization44/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23726.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23726.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization45/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23730.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23730.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization46/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23734.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23734.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization47/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23738.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23738.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization48/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23742.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23742.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization49/control' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23746.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23746.param'

python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization25/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23650.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23650.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization26/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23654.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23654.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization27/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23658.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23658.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization28/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23662.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23662.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization29/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23666.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23666.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization30/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23670.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23670.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization31/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23674.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23674.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization32/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23678.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23678.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization33/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23682.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23682.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization34/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23686.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23686.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization35/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1527_23690.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1527_23690.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization36/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23694.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23694.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization37/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23698.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23698.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization38/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23702.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23702.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization39/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23706.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23706.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization40/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23710.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23710.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization41/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23714.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23714.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization42/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23718.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23718.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization43/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23722.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23722.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization44/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23726.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23726.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization45/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23730.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23730.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization46/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23734.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23734.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization47/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23738.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23738.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization48/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23742.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23742.param'
python $modeldir'/evoked_activity_ex_no_variability.py' $modeldir'/Figure4C_2/results/realization49/L1inact' $cellModelName $ongoingUpL1inact $ongoingDownL1inact $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_up_state_functional_map_20130326-1528_23746.param' $modeldir'/Figure4C_2/evoked_activity_new_boutons_L1inact_down_state_functional_map_20130326-1528_23746.param'

echo "*********************"
echo "Analysis: PSP SD"
python $modeldir"/analysis/analyze_psp_std.py" $modeldir"/Figure4C_2/results" $modeldir"/Figure4C_2/results/summary_psp_std.csv" 50 1
