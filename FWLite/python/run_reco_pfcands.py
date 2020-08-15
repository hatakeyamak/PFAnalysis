# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from ctypes import c_uint8

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.AutoLibraryLoader.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

#genparticles, genparLabel = Handle("std::vector<reco::GenParticle>"), "genParticles"
#caloparticles, caloparLabel = Handle("std::vector<CaloParticle>"), "mix:MergedCaloTruth"
#vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "inclusiveSecondaryVertices"
#pfcands, pfcandLabel = Handle("std::vector<reco::PFCandidate>"), "particleFlowEGamma"
#vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlinePrimaryVertices"
pfcands, pfcandLabel = Handle("std::vector<reco::PFCandidate>"), "particleFlow"
#pfclusterhf, pfclusterhfLabel = Handle("std::vector<reco::PFCluster>"), "particleFlowClusterHF"
#pfrechithf, pfrechithfLabel = Handle("std::vector<reco::PFRecHit>"), "particleFlowRecHitHF"
#pftracks, pftrackLabel = Handle("std::vector<reco::PFRecTrack>"), "pfTrack"

gedPhotons, gedPhotonLabel = Handle("std::vector<reco::Photon>"), "gedPhotons"
gedGsfElectrons, gedGsfElectronLabel = Handle("std::vector<reco::GsfElectron>"), "gedGsfElectrons"
#svcands, svcandLabel = Handle("std::vector<reco::VertexCompositePtrCandidate>"), "inclusiveCandidateSecondaryVertices"

#pfcandPtScore = Handle("edm::ValueMap<float>")
#verticesScore = Handle("edm::ValueMap<float>")

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
#events = Events('file:../inputFiles/F6F192B6-19E8-1245-913A-545FF4D712A1.root')
events = Events('file:step3_20evt.root')
#events = Events('file:step3_20evt_org.root')
#events = Events('file:step3_rerereco_all.root')
#events = Events('file:step3.root')
#events = Events('file:step3_withparticleFlowTmpBarrel.root')

for iev,event in enumerate(events):
    if iev >= 20: break 
    #event.getByLabel(vertexLabel, vertices)
    #event.getByLabel(vertexLabel, verticesScore)
    event.getByLabel(pfcandLabel, pfcands)
    event.getByLabel(gedPhotonLabel, gedPhotons)
    event.getByLabel(gedGsfElectronLabel, gedGsfElectrons)
    #event.getByLabel(svcandLabel, svcands)
    #event.getByLabel(pfcandHcalDepthLabel,hcalDepthScore)
    #event.getByLabel(pfcandPtLabel,pfcandPtScore)

    # pfcands
    npf=0
    for i,j in enumerate(pfcands.product()):
        npf=npf+1
        print "pfcands: pt %17.13f eta %18.14f pdgId %5d " % ( j.pt(), j.eta(), j.pdgId())
    print "npfcand: ",npf

    # photons
    npf=0
    for i,j in enumerate(gedPhotons.product()):
        npf=npf+1
        print "gedPhotons: pt %17.13f eta %18.14f pdgId %5d " % ( j.pt(), j.eta(), j.pdgId())
    print "nphoton: ",npf

    # electrons
    npf=0
    for i,j in enumerate(gedGsfElectrons.product()):
        npf=npf+1
        print "gedGsfElectrons: pt %17.13f eta %18.14f pdgId %5d " % ( j.pt(), j.eta(), j.pdgId())
    print "nelectron: ",npf

    # taus

    # jets

    # met
    
    
