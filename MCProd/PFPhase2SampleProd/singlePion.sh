#!/bin/bash
echo "Starting the job"
hostname
export X509_USER_PROXY=$PWD/x509up
export SCRAM_ARCH=slc7_amd64_gcc900
source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_12_1_0_pre3/src ] ; then
  echo release CMSSW_12_1_0_pre3 already exists
else
  scram p CMSSW CMSSW_12_1_0_pre3
fi
cd CMSSW_12_1_0_pre3/src
eval `scram runtime -sh`

filenum=$1
jobnum=$2
nevents=$3
# use 10000 events/lumi
declare -i lumi=$filenum*10000/$nevents+1
declare -i firstEvent=$filenum*$nevents
echo "Event" $firstEvent "lumi" $lumi
basepath=$4
# Download fragment from McM
#curl -s -k https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/JME-RunIIWinter19PFCalib16GS-00001 --retry 3 --create-dirs -o Configuration/GenProduction/python/JME-RunIIWinter19PFCalib16GS-00001-fragment.py
#[ -s Configuration/GenProduction/python/JME-RunIIWinter19PFCalib16GS-00001-fragment.py ] || exit $?;
#scram b
#cd ../..

gsfile=step_gs_${filenum}.root
digifile=step_digi_${filenum}.root
recofile=step_reco_SinglePi_Eta2p7_3p5_${jobnum}_${filenum}.root

# reco
cmsDriver.py SinglePiPt25Eta1p7_2p7_cfi --python_filename=step_gs.py -s GEN,SIM \
    --conditions auto:phase2_realistic_T21 --beamspot HLLHC --datatier GEN-SIM --eventcontent FEVTDEBUG \
    --geometry Extended2026D76 --era Phase2C11I13M9 \
    --nThreads 8 \
    --fileout file:$gsfile -n $nevents \
    --customise_commands "process.generator.PGunParameters.MinEta=2.7\nprocess.generator.PGunParameters.MaxEta=3.5\nprocess.generator.PGunParameters.MaxPt=200.01\nprocess.generator.PGunParameters.MinPt=199.99\nfrom IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper;randHelper=RandomNumberServiceHelper(process.RandomNumberGeneratorService);randHelper.populate()\nprocess.source.firstLuminosityBlock=cms.untracked.uint32(${lumi})\nprocess.source.firstEvent=cms.untracked.uint32(${firstEvent})"

# digi step
cmsDriver.py step_digi -s DIGI:pdigi_valid,L1TrackTrigger,L1,DIGI2RAW,HLT:@fake2 \
    --conditions auto:phase2_realistic_T21 --datatier GEN-SIM-DIGI-RAW --eventcontent FEVTDEBUGHLT \
    --geometry Extended2026D76 --era Phase2C11I13M9 \
    --nThreads 8 \
    --filein file:$gsfile --fileout file:$digifile -n -1

# reco step
cmsDriver.py step_reco -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT \
    --conditions auto:phase2_realistic_T21 --datatier GEN-SIM-RECO,MINIAODSIM \
    --eventcontent FEVTDEBUGHLT,MINIAODSIM --geometry Extended2026D76 --era Phase2C11I13M9 \
    --nThreads 8 \
    --filein file:$digifile --fileout file:$recofile -n -1 $customizeReco

#xrdcp $recofile root://eoscms.cern.ch/${basepath}/$recofile
eval `scram unsetenv -sh`; gfal-copy $recofile gsiftp://kodiak-se.baylor.edu//cms/data/${basepath}/$recofile

rm *.root
