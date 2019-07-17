#
# Usage example:
# python plot_reco_pfHN.py
#
# import ROOT in batch mode
import sys
import math
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
from ROOT import TF1, TF2, TH1, TH2, TH2F, TProfile, TAxis, TMath, TEllipse, TStyle, TFile, TColor, TSpectrum, TCanvas, TPad, TVirtualFitter, gStyle
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from ctypes import c_uint8

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.FWLiteEnabler.enable() 

# Create histograms, etc.
ROOT.gROOT.SetStyle('Plain') # white background
H_NPV = ROOT.TH1F ("NPV","NPV",101,-0.5,100.5)
H_GenPi_Pt = ROOT.TH1F("H_GenPi_Pt","H_GenPi_Pt",50,0.,50.)
H_GenPi_Eta = ROOT.TH1F("H_GenPi_Eta","H_GenPi_Eta",100,-5.0,5.0)

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

# 
# Import CMS python class definitions such as Process, Source, and EDProducer
import FWCore.ParameterSet.Config as cms
ROOT.gSystem.Load("libDataFormatsHcalDetId.so")
#from DataFormats.HcalDetId import DetId

#muons, muonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"
#electrons, electronLabel = Handle("std::vector<pat::Electron>"), "slimmedElectrons"
#photons, photonLabel = Handle("std::vector<pat::Photon>"), "slimmedPhotons"
#taus, tauLabel = Handle("std::vector<pat::Tau>"), "slimmedTaus"
#jets, jetLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
#fatjets, fatjetLabel = Handle("std::vector<pat::Jet>"), "slimmedJetsAK8"
#mets, metLabel = Handle("std::vector<pat::MET>"), "slimmedMETs"
#vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlineSlimmedPrimaryVertices"
#
genparticles, genparLabel = Handle("std::vector<reco::GenParticle>"), "genParticles"
caloparticles, caloparLabel = Handle("std::vector<CaloParticle>"), "mix:MergedCaloTruth"
vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlinePrimaryVertices"
pfcands, pfcandLabel = Handle("std::vector<reco::PFCandidate>"), "particleFlow"
pfclusterhf, pfclusterhfLabel = Handle("std::vector<reco::PFCluster>"), "particleFlowClusterHF"
pfrechithf, pfrechithfLabel = Handle("std::vector<reco::PFRecHit>"), "particleFlowRecHitHF"
pftracks, pftrackLabel = Handle("std::vector<reco::PFRecTrack>"), "pfTrack"

pfcandPtScore = Handle("edm::ValueMap<float>")
verticesScore = Handle("edm::ValueMap<float>")

if len(sys.argv)>1:
    output=sys.argv[1]
else:
    output="PF_HF_SinglePi_Pt20"
        
# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
if output == "PF_HF_SinglePi_Pt20":
    events = Events('../29088.0_SinglePiPt25Eta1p7_2p7+SinglePiPt25Eta1p7_2p7_2023D41_GenSimHLBeamSpotFull+DigiFullTrigger_2023D41+RecoFullGlobal_2023D41+HARVESTFullGlobal_2023D41/step3.root')
# elif output == "NH_CaloEnergy_SingleMuon_2016C":
#     events = Events('root://cmsxrootd.fnal.gov//store/data/Run2016C/SingleMuon/AOD/07Aug17-v1/110001/44B1F7D7-3F80-E711-AE45-001E67E6F4A9.root')
# elif output == "NH_CaloEnergy_SingleMuon_2017C":
#     events = Events('root://cmsxrootd.fnal.gov//store/data/Run2017C/SingleMuon/AOD/17Nov2017-v1/70000/A073FB9C-6CDA-E711-9551-02163E01A465.root')
# elif output == "NH_CaloEnergy_SingleMuon_2018C":
#     events = Events('root://cmsxrootd.fnal.gov//store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00001/56DF33B4-0A6E-2747-99BC-3FD669B097FF.root')
    
for iev,event in enumerate(events):
    if iev >= 100: break 
    #event.getByLabel(muonLabel, muons)
    #event.getByLabel(electronLabel, electrons)
    #event.getByLabel(photonLabel, photons)
    #event.getByLabel(tauLabel, taus)
    #event.getByLabel(jetLabel, jets)
    #event.getByLabel(fatjetLabel, fatjets)
    #event.getByLabel(metLabel, mets)
    event.getByLabel(genparLabel, genparticles)
    event.getByLabel(caloparLabel, caloparticles)
    event.getByLabel(vertexLabel, vertices)
    event.getByLabel(pfcandLabel, pfcands)
    event.getByLabel(pfclusterhfLabel, pfclusterhf)
    event.getByLabel(pfrechithfLabel, pfrechithf)
    event.getByLabel(pftrackLabel, pftracks)
    #print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())

    # Vertices 
    H_NPV.Fill(len(vertices.product()))
    print "len(vertices.product())", len(vertices.product()) 
    if len(vertices.product()) == 0 or vertices.product()[0].ndof() < 4:
        print "Event has no good primary vertex."
        #continue
    else:
        PV = vertices.product()[0]
        print "PV at x,y,z = %+5.3f, %+5.3f, %+6.3f, ndof: %.1f " % (PV.x(), PV.y(), PV.z(), PV.ndof())

    # Gen particles
    for i,j in enumerate(genparticles.product()):  # loop over gen candidates
        print "GenParticle: iev %3d genpars %3d: pt %5.1f eta %5.2f phi %5.2f " % ( iev, i, j.pt(), j.eta(), j.phi() )

    # Calo particles
    for i,j in enumerate(caloparticles.product()):  # loop over calo candidates
        print "CaloParticle: iev %3d calopars %3d: pt %5.1f eta %5.2f phi %5.2f depth=" % ( iev, i, j.pt(), j.eta(), j.phi() )
             
    # PF HF clusters    
    for i,j in enumerate(pfclusterhf.product()):  # loop over pf clusters
        pfc = j
        # layer 11:, 12:
        print "HFCluster: iev %3d pfclus %3d: pt %5.1f eta %5.2f phi %5.2f layer %3d" % ( iev, i, pfc.pt(), pfc.eta(), pfc.phi(), pfc.layer() )
        fracs = pfc.recHitFractions()
        hfracs = pfc.hitsAndFractions()
        for e in range(0, fracs.size()):  # loop over rechits
            #print "  ", fracs[e], hfracs[e], hfracs[e].first, hfracs[e].second
            id = hfracs[e].first.rawId()
            print "  ", id, pfc.recHitFractions()[e].recHitRef().detId()
            #print HcalDetId(id).ieta(),HcalDetId(id).iphi(),HcalDetId(id).depth()
            # extracting ieta,iphi from detid: not working

    # PF HF rechits    
    for i,j in enumerate(pfrechithf.product()):  # loop over pf rechits
        # print "PFRecHitHF: iev %3d pfrechits %3d: energy %5.1f pt %5.1f detid %15d eta %5.2f phi %5.2f layer %3d" % ( iev, i, j.energy(), math.sqrt(j.pt2()), j.detId(), j.positionREP().eta(), j.positionREP().eta(), j.layer() )
        # access to HcalDetId and its geometrical position: not working.
        print "PFRecHitHF: iev %3d pfrechits %3d: energy %5.1f detid %15d layer %3d" % ( iev, i, j.energy(), j.detId(), j.layer() )        
        
    # PF tracks
    for i,j in enumerate(pftracks.product()):  # loop over pf rechits
        print "PFTrack: iev %3d pftrack %3d: pt %5.1f eta %5.2f phi %5.2f charge %3d" % ( iev, i, j.trackRef().pt(), j.trackRef().eta(), j.trackRef().phi(), j.charge() )
        j.calculatePositionREP()
        #print j.extrapolatedPoint(9) # get information on extrapolated points: not working
        
# Set up canvas : 
w = 1400 
h =  700
can  = ROOT.TCanvas("can", "histograms   ", w, h)

#####
ROOT.gPad.SetLogy()
ROOT.gPad.SetLogy(0)
H_NPV.Draw()
can.SaveAs(output+"_NPV.pdf")
can.SaveAs(output+"_NPV.png")
can.SaveAs(output+"_NPV.root")

