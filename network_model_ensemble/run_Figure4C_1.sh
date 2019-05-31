#!/bin/bash

modeldir=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)
source $modeldir"/pythonpath.sh"
cd $modeldir"/Figure4C_1"

cellModelName=$modeldir"/../cell_model/waters_l23_ar_dend.param"

ongoingUp=$modeldir"/Figure4C_1/up_state_parameters_new_boutons.param"
ongoingDown=$modeldir"/Figure4C_1/down_state_parameters_new_boutons.param"

python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization00/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1526_23551.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1526_23551.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization01/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1526_23554.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1526_23554.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization02/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1526_23558.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1526_23558.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization03/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23562.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23562.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization04/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23566.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23566.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization05/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23570.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23570.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization06/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23574.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23574.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization07/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23578.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23578.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization08/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23582.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23582.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization09/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23586.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23586.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization10/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23590.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23590.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization11/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23594.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23594.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization12/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23598.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23598.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization13/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23602.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23602.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization14/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23606.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23606.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization15/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23610.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23610.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization16/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23614.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23614.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization17/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23618.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23618.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization18/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23622.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23622.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization19/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23626.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23626.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization20/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23630.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23630.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization21/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23634.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23634.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization22/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23638.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23638.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization23/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23642.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23642.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization24/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23646.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23646.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization25/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23650.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23650.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization26/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23654.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23654.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization27/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23658.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23658.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization28/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23662.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23662.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization29/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23666.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23666.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization30/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23670.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23670.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization31/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23674.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23674.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization32/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23678.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23678.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization33/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23682.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23682.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization34/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23686.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23686.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization35/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23690.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23690.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization36/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23694.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23694.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization37/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23698.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23698.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization38/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23702.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23702.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization39/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23706.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23706.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization40/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23710.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23710.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization41/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23714.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23714.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization42/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23718.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23718.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization43/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23722.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23722.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization44/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23726.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23726.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization45/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23730.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23730.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization46/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23734.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23734.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization47/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23738.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23738.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization48/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23742.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23742.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization49/control1' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23746.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23746.param'

python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization00/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1526_23551.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1526_23551.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization01/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1526_23554.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1526_23554.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization02/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1526_23558.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1526_23558.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization03/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23562.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23562.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization04/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23566.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23566.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization05/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23570.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23570.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization06/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23574.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23574.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization07/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23578.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23578.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization08/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23582.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23582.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization09/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23586.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23586.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization10/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23590.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23590.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization11/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23594.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23594.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization12/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23598.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23598.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization13/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23602.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23602.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization14/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23606.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23606.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization15/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23610.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23610.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization16/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23614.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23614.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization17/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23618.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23618.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization18/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23622.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23622.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization19/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23626.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23626.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization20/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23630.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23630.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization21/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23634.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23634.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization22/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23638.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23638.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization23/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23642.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23642.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization24/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23646.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23646.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization25/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23650.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23650.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization26/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23654.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23654.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization27/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23658.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23658.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization28/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23662.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23662.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization29/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23666.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23666.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization30/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23670.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23670.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization31/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23674.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23674.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization32/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23678.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23678.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization33/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23682.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23682.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization34/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23686.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23686.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization35/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1527_23690.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1527_23690.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization36/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23694.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23694.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization37/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23698.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23698.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization38/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23702.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23702.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization39/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23706.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23706.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization40/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23710.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23710.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization41/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23714.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23714.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization42/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23718.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23718.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization43/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23722.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23722.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization44/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23726.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23726.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization45/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23730.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23730.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization46/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23734.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23734.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization47/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23738.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23738.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization48/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23742.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23742.param'
python $modeldir'/evoked_activity_inh_no_variability.py' $modeldir'/Figure4C_1/results/realization49/control2' $cellModelName $ongoingUp $ongoingDown $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_up_state_functional_map_20130326-1528_23746.param' $modeldir'/Figure4C_1/evoked_activity_new_boutons_control_down_state_functional_map_20130326-1528_23746.param'


echo "*********************"
echo "Analysis: PSP SD"
python $modeldir"/analysis/analyze_psp_std.py" $modeldir"/Figure4C_1/results" $modeldir"/Figure4C_1/results/summary_psp_std.csv" 50 2