
#
# Number of events
###########
#
EVENTS=1000

#
# PT range
###########
#
Dataset=0p01to10
DAS="dbs:/DoublePhoton_FlatPT-0p01to10_13p6TeV/Run3Winter24GS-NoMaterial_133X_mcRun3_2024_realistic_v8-v4/GEN-SIM"

#Dataset=10to500
#DAS="dbs:/DoublePhoton_FlatPT-10to500_13p6TeV/Run3Winter24GS-NoMaterial_133X_mcRun3_2024_realistic_v8-v4/GEN-SIM"

#
# PU scenarios
###########
#
#Pileup=NoPileUp
#PUName=NoPileUp
# also adjust customise_commands of DIGI step
#
Pileup=Flat_10_50_25ns
PUName=FlatPU0to80ZM
#
#customDIGI="process.RAWSIMoutput.outputCommands.extend(['keep *_*articleFlowCluster*_*_*'])"
#customDIGI="process.RAWSIMoutput.outputCommands.extend(['keep  *_*articleFlowCluster*_*_*']); process.mix.input.nbPileupEvents.probFunctionVariable = cms.vint32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80) \n process.mix.input.nbPileupEvents.probValue = cms.vdouble(0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457)" --pileup_input "dbs:/MinBias_TuneCP5_13p6TeV-pythia8/Run3Winter24GS-133X_mcRun3_2024_realistic_v7-v1/GEN-SIM"

#
# MCM entries
# https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=DoublePhoton_FlatPT-0p01to10_13p6TeV
# https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=DoublePhoton_FlatPT-10to500_13p6TeV
# https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=DoublePhoton_FlatPT-500to1000_13p6TeV&prepid=*Winter24*
# https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=DoublePhoton_FlatPT-1000to1500_13p6TeV&prepid=*Winter24*
# and make sure to choose the sample with GS --geometry DB:ZeroMaterial
#

#-----
# DIGI
#-----
#cmsDriver.py  --python_filename EGM-Run3Winter24Digi-00057_KH_cfg.py --eventcontent RAWSIM --pileup NoPileUp \
#	      --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW --fileout file:EGM-Run3Winter24Digi-00057.root \
#	      --conditions 133X_mcRun3_2024_realistic_v9 --step DIGI,L1,DIGI2RAW,HLT:2023v12 --geometry DB:ZeroMaterial \
#	      --filein "file:EGM-Run3Winter24GS-00002.root" \
#	      --era Run3_2023 --no_exec --mc -n -1

# Keeping *_*articleFlowCluster*_*_* in outputs

#         --pileup Flat_10_50_25ns \
#	 --customise_commands "process.mix.input.nbPileupEvents.probFunctionVariable = cms.vint32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80) \n process.mix.input.nbPileupEvents.probValue = cms.vdouble(0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457)"	 
#	 --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW --fileout file:EGM-Run3Winter24Digi-00037.root \
#	 --pileup_input "dbs:/MinBias_TuneCP5_13p6TeV-pythia8/Run3Winter24GS-133X_mcRun3_2024_realistic_v7-v1/GEN-SIM" \
#	 --conditions 133X_mcRun3_2024_realistic_v9
	 
cmsDriver.py  --python_filename EGM-Run3Winter24Digi-${Dataset}-${PUName}_KH_cfg.py --eventcontent RAWSIM --pileup ${Pileup} \
	      --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW --fileout file:EGM-Run3Winter24Digi-${Dataset}-${PUName}.root \
	      --conditions 133X_mcRun3_2024_realistic_v9 --step DIGI,L1,DIGI2RAW,HLT:2023v12 --geometry DB:ZeroMaterial \
	      --filein ${DAS} \
	      --era Run3_2023 --no_exec --mc -n ${EVENTS} --nThreads 8 \
	      --customise_commands "process.RAWSIMoutput.outputCommands.extend(['keep  *_*articleFlowCluster*_*_*']); process.mix.input.nbPileupEvents.probFunctionVariable = cms.vint32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80) \n process.mix.input.nbPileupEvents.probValue = cms.vdouble(0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457,0.0123457)" --pileup_input "dbs:/MinBias_TuneCP5_13p6TeV-pythia8/Run3Winter24GS-133X_mcRun3_2024_realistic_v7-v1/GEN-SIM"
	      #--customise_commands "process.RAWSIMoutput.outputCommands.extend(['keep *_*articleFlowCluster*_*_*'])"

cmsRun EGM-Run3Winter24Digi-${Dataset}-${PUName}_KH_cfg.py

echo ${customDIGI}

#-----
# RECO
#-----
# cmsDriver.py  --python_filename EGM-Run3Winter24Reco-00057_KH_cfg.py --eventcontent AODSIM \
# 	      --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM \
# 	      --fileout file:EGM-Run3Winter24Reco-00057.root --conditions 133X_mcRun3_2024_realistic_v9 \
# 	      --customise_commands "process.AODSIMoutput.outputCommands.extend(['keep  *_particleFlowCluster*_*_*', 'keep *_mix_MergedTrackTruth_*', 'keep *_ecalDigis*_*_*', 'keep *_genParticles*_*_*']); process.AODSIMoutput.outputCommands.append('keep *_particleFlowClusterECAL_*_*'); process.AODSIMoutput.outputCommands.append('keep *_mix_MergedTrackTruth_*'); process.AODSIMoutput.outputCommands.append('keep EBDigiCollection_ecalDigis*_*_*'); process.AODSIMoutput.outputCommands.append('keep EEDigiCollection_ecalDigis*_*_*'); process.AODSIMoutput.outputCommands.append('keep *SrFlagsSorted*_ecalDigis*_*_*'); process.AODSIMoutput.outputCommands.append('keep recoGenParticle*_genParticles_*_*')" --step RAW2DIGI,L1Reco,RECO,RECOSIM --geometry DB:ZeroMaterial \
# 	      --filein "file:EGM-Run3Winter24Digi-00057.root" \
# 	      --era Run3_2023 --no_exec --mc -n -1

# Keeping *_*articleFlowCluster*_*_* in outputs
cmsDriver.py  --python_filename EGM-Run3Winter24Reco-${Dataset}-${PUName}_KH_cfg.py --eventcontent AODSIM \
	      --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM \
	      --fileout file:EGM-Run3Winter24Reco-${Dataset}-${PUName}.root --conditions 133X_mcRun3_2024_realistic_v9 \
	      --customise_commands "process.AODSIMoutput.outputCommands.extend(['keep  *_particleFlowCluster*_*_*', 'keep *_mix_MergedTrackTruth_*', 'keep *_ecalDigis*_*_*', 'keep *_genParticles*_*_*']); process.AODSIMoutput.outputCommands.append('keep *_*articleFlowCluster*_*_*'); process.AODSIMoutput.outputCommands.append('keep *_mix_MergedTrackTruth_*'); process.AODSIMoutput.outputCommands.append('keep EBDigiCollection_ecalDigis*_*_*'); process.AODSIMoutput.outputCommands.append('keep EEDigiCollection_ecalDigis*_*_*'); process.AODSIMoutput.outputCommands.append('keep *SrFlagsSorted*_ecalDigis*_*_*'); process.AODSIMoutput.outputCommands.append('keep recoGenParticle*_genParticles_*_*')" \
	      --step RAW2DIGI,L1Reco,RECO,RECOSIM --geometry DB:ZeroMaterial \
	      --filein "file:EGM-Run3Winter24Digi-${Dataset}-${PUName}.root" \
	      --era Run3_2023 --no_exec --mc -n -1 --nThreads 8

cmsRun EGM-Run3Winter24Reco-${Dataset}-${PUName}_KH_cfg.py

#-----
# PAT
#-----

# Keeping *_*articleFlowCluster*_*_* in outputs
cmsDriver.py  --python_filename EGM-Run3Winter24MiniAOD-${Dataset}-${PUName}_KH_cfg.py --eventcontent MINIAODSIM \
	      --customise Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM \
	      --fileout file:EGM-Run3Winter24MiniAOD-${Dataset}-${PUName}.root --conditions 133X_mcRun3_2024_realistic_v9 \
	      --customise_commands "process.MINIAODSIMoutput.outputCommands.extend(['keep *_*articleFlowCluster*_*_*', 'keep *_mix_MergedTrackTruth_*', 'keep  *_ecalDigis*_*_*', 'keep *_genParticles*_*_*']); process.MINIAODSIMoutput.outputCommands.append('keep *_particleFlowClusterECAL_*_*'); process.MINIAODSIMoutput.outputCommands.append('keep *_mix_MergedTrackTruth_*'); process.MINIAODSIMoutput.outputCommands.append('keep EBDigiCollection_ecalDigis*_*_*'); process.MINIAODSIMoutput.outputCommands.append('keep EEDigiCollection_ecalDigis*_*_*'); process.MINIAODSIMoutput.outputCommands.append('keep *SrFlagsSorted*_ecalDigis*_*_*'); process.MINIAODSIMoutput.outputCommands.append('keep recoGenParticle*_genParticles_*_*')" \
	      --step PAT --geometry DB:ZeroMaterial \
	      --filein "file:EGM-Run3Winter24Reco-${Dataset}-${PUName}.root" \
	      --era Run3_2023 --no_exec --mc -n -1 --nThreads 8

cmsRun EGM-Run3Winter24MiniAOD-${Dataset}-${PUName}_KH_cfg.py
