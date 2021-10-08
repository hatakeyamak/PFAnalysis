#
# Usage example:
# python run_simtrack.py
#
# import ROOT in batch mode
import sys
import math
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from ctypes import c_uint8

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.FWLiteEnabler.enable() 

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

# 
# Import CMS python class definitions such as Process, Source, and EDProducer
import FWCore.ParameterSet.Config as cms
#
genparticles, genparLabel = Handle("std::vector<reco::GenParticle>"), "genParticles"
simtracks, simtrackLabel = Handle("std::vector<SimTrack>"), "g4SimHits"
caloparticles, caloparLabel = Handle("std::vector<CaloParticle>"), "mix:MergedCaloTruth"
vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlinePrimaryVertices"
pfcands, pfcandLabel = Handle("std::vector<reco::PFCandidate>"), "particleFlow"
pfclusterhf, pfclusterhfLabel = Handle("std::vector<reco::PFCluster>"), "particleFlowClusterHF"
pfrechithf, pfrechithfLabel = Handle("std::vector<reco::PFRecHit>"), "particleFlowRecHitHF"
pftracks, pftrackLabel = Handle("std::vector<reco::PFRecTrack>"), "pfTrack"

# Phase2 SinglePion PU200
#events = Events('root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRWinter20DIGI/SinglePion_PT0to200/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/10000/DD90F2CA-90F5-E449-96AC-6A868594B25E.root')
# Phase2 SinglePion NoPU
#events = Events('root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRWinter20DIGI/SinglePion_PT0to200/GEN-SIM-DIGI-RAW/NoPU_110X_mcRun4_realistic_v3-v2/50000/E3B81E3A-9B35-A54B-A6E0-2A4D41CD83A5.root')
# Run 3
events = Events('root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_0_0_patch1/RelValSinglePiPt60GeV/GEN-SIM-DIGI-RAW/110X_mcRun3_2021_realistic_v6-v1/10000/BA6D7F51-9A29-0C4D-8480-08AF3874C33F.root')

for iev,event in enumerate(events):
    if iev >= 100: break 
    event.getByLabel(genparLabel, genparticles)
    #event.getByLabel(caloparLabel, caloparticles)
    event.getByLabel(simtrackLabel, simtracks)
    #event.getByLabel(vertexLabel, vertices)
    #event.getByLabel(pfcandLabel, pfcands)
    #event.getByLabel(pfclusterhfLabel, pfclusterhf)
    #event.getByLabel(pfrechithfLabel, pfrechithf)
    #event.getByLabel(pftrackLabel, pftracks)

    # Vertices 
    # H_NPV.Fill(len(vertices.product()))
    # print "len(vertices.product())", len(vertices.product()) 
    # if len(vertices.product()) == 0 or vertices.product()[0].ndof() < 4:
    #     print "Event has no good primary vertex."
    #     #continue
    # else:
    #     PV = vertices.product()[0]
    #     print "PV at x,y,z = %+5.3f, %+5.3f, %+6.3f, ndof: %.1f " % (PV.x(), PV.y(), PV.z(), PV.ndof())

    print("== event: ",iev)
    
    # Gen particles
    for i,j in enumerate(genparticles.product()):  # loop over gen candidates
        print("GenParticle: %3d: pt %5.1f eta %5.2f phi %5.2f pdgid %6d " % ( i, j.pt(), j.eta(), j.phi(), j.pdgId() ))
        
    # Simtracks
    for i,j in enumerate(simtracks.product()):  # loop over gen candidates
        if j.genpartIndex() > 0.:
            print("SimTrack:    %3d: pt %5.1f eta %5.2f phi %5.2f pdgid %6d genparIdx %4d" \
                % ( i, j.momentum().pt(), j.momentum().eta(), j.momentum().phi(), j.type(), j.genpartIndex() ))
            print("  trackerSurface position R %5.2f Z %5.2f " \
                % ( j.trackerSurfacePosition().R(), j.trackerSurfacePosition().Z() ))
            print("  trackerSurface momentum pt %5.1f eta %5.2f phi %5.2f " \
                % ( j.trackerSurfaceMomentum().pt(), j.trackerSurfaceMomentum().eta(), j.trackerSurfaceMomentum().phi() ))

    # Calo particles
    # for i,j in enumerate(caloparticles.product()):  # loop over calo candidates
    #     if j.status()>0.:
    #         print "CaloParticle:  %3d: pt %5.1f eta %5.2f phi %5.2f status %4d " % ( i, j.pt(), j.eta(), j.phi(), j.status() )

    # PF HF clusters    
    # for i,j in enumerate(pfclusterhf.product()):  # loop over pf clusters
    #     pfc = j
    #     # layer 11:, 12:
    #     print "HFCluster: iev %3d pfclus %3d: pt %5.1f eta %5.2f phi %5.2f layer %3d" % ( iev, i, pfc.pt(), pfc.eta(), pfc.phi(), pfc.layer() )
    #     fracs = pfc.recHitFractions()
    #     hfracs = pfc.hitsAndFractions()
    #     for e in range(0, fracs.size()):  # loop over rechits
    #         id = hfracs[e].first.rawId()
    #         print "  ", id, pfc.recHitFractions()[e].recHitRef().detId()

    # PF HF rechits    
    # for i,j in enumerate(pfrechithf.product()):  # loop over pf rechits
    #     print "PFRecHitHF: iev %3d pfrechits %3d: energy %5.1f detid %15d layer %3d" % ( iev, i, j.energy(), j.detId(), j.layer() )      
    # PF tracks
    # for i,j in enumerate(pftracks.product()):  # loop over pf rechits
    #     print "PFTrack: iev %3d pftrack %3d: pt %5.1f eta %5.2f phi %5.2f charge %3d" % ( iev, i, j.trackRef().pt(), j.trackRef().eta(), j.trackRef().phi(), j.charge() )
        

