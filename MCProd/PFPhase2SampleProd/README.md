# Quick start

# Local test
Check
- CMSSW release
- cmsDriver settings
- output locations
- number of events
in shell script and jdl file (e.g. singlePion.sh & submitSinglePion_noPU_12_1.jdl)

Make sure the output directory is prepared.
~~~
bash singlePion.sh 1 101 10 /store/user/hatake/PF/CMSSW_12_1_0_pre3_HF/Run3_2021 >& run_1.log &
# 1st argument: job index
# 2nd argument: dummy job ID. meaningful only for condor jobs
# 3rd argument: # of events
# 4th argument: output location (e.g. under /cms/data/store/user/..., but start from /store/user). Check also prefix such as gsiftp://kodiak-se.baylor.edu//cms/data/ in shell script.
~~~

# Batch job submission

- Create a grid proxy file to be shipped
- Create condor log file storage directory
- Submit!
~~~
export X509_USER_PROXY=/tmp/x509up_u$(id -u); voms-proxy-init -voms cms -valid 192:00; cp $X509_USER_PROXY x509up
mkdir logs_singlePion # depends on the directory name specified in submission jdl script
condor_submit submitSinglePion_noPU_12_1.jdl
~~~
