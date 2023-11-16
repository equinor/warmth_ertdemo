#
# prepare
#

git clone git@github.com:equinor/warmth_ertdemo.git
cd warmth_ertdemo

source /prog/res/komodo/stable/enable
komodoenv _env
source _env/enable
pip install .

#
# step0 - demonstrate a warmth simulation
# 

cd step0
python mywell.py
cd ..

#
# Step 1 - using ERT just to run 5 realizations
#

cd step1
ls -al
# show mywell.ert
# show MYWELL_WVAL
chmod +x mywell_eval.py
ert gui mywell.ert
# run ensemble
tree mywell_out
ls -al mywell_out/re*/it*/mywell_temp.out
head -5 mywell_out/re*/it*/mywell_temp.out
rm -rf storage mywell_out logs .ert_runpath_list 
cd ..

#
# Step 2 - let ERT generate parameters
#

cd step2
ls -al
ert ensemble_experiment mywell.ert
# show mywell.ert
# show mywell_eval.py
# show params.tmpl
# show params_priors
tree mywell_out
head mywell_out/re*/it*/params.json
ert gui mywell.ert
rm -rf storage mywell_out logs .ert_runpath_list 
cd ..

#
# Step 3 - add observations and do history matching with ES-MDA
#

cd step3
ert es_mda mywell.ert
# show mywell_temp_data.txt
# show observations
# show mywell.ert
ert gui mywell.ert
rm -rf storage mywell_out logs .ert_runpath_list 
cd ..
