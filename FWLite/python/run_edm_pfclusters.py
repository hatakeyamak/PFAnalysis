import math

# Define a few useful functions
#===============================
def deltaR2( e1, p1, e2=None, p2=None):
    """Take either 4 arguments (eta,phi, eta,phi) or two objects that have 'eta', 'phi' methods)"""
    if (e2 == None and p2 == None):
        return deltaR2(e1.eta(),e1.phi(), p1.eta(), p1.phi())
    de = e1 - e2
    dp = deltaPhi(p1, p2)
    return de*de + dp*dp

def deltaR( *args ):
    return math.sqrt( deltaR2(*args) )

def deltaPhi( p1, p2):
    '''Computes delta phi, handling periodic limit conditions.'''
    res = p1 - p2
    while res > math.pi:
        res -= 2*math.pi
    while res < -math.pi:
        res += 2*math.pi
    return res

# import ROOT in batch mode
#===============================
import sys
if (len(sys.argv)-1)!=3:
    print("Need two arguements: for example, python3 run_edm_pfclusters.py 0p01to10/10to500 NoPileUp/FlatPU0to80ZM Custom/Default")
PTrange = sys.argv[1]
PUscenario = sys.argv[2]
Option = sys.argv[3]
print(PTrange,PUscenario,Option)

debug=False

import ROOT
#........................ 20/2/2024
from ROOT import TF1, TF2, TH1, TH2, TH2F, TProfile, TAxis, TMath, TEllipse, TStyle, TFile, TColor, TSpectrum, TCanvas, TPad, TVirtualFitter, gStyle
#........................
ROOT.gROOT.SetBatch(True)

# load FWLite C++ libraries
#===============================
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.AutoLibraryLoader.enable()

# Create histograms, etc.
#.................................. 20/2/2024
ROOT.gROOT.SetStyle('Plain') # white background
#H_NPV = ROOT.TH1F ("NPV","NPV",101,-0.5,100.5)
H_GenGM_Pt = ROOT.TH1F("H_GenPi_Pt","H_GenPi_Pt",50,0.,200.)
H_GenGM_Eta = ROOT.TH1F("H_GenPi_Eta","H_GenPi_Eta",100,-5.0,5.0)
H_GenGM_Pt_vs_Resp = ROOT.TProfile("H_GenPi_Pt_vs_Resp","H_GenPi_Pt_vs_Resp",50,0.,200.,0.8,1.2)
H_GenGM_Pt_vs_RespE = ROOT.TProfile("H_GenPi_Pt_vs_RespE","H_GenPi_Pt_vs_RespE",50,0.,200.,0.8,1.2)
H_GenGM_Pt_vs_RespCE = ROOT.TProfile("H_GenPi_Pt_vs_RespCE","H_GenPi_Pt_vs_RespCE",50,0.,200.,0.8,1.2)
H_GenGM_Pt_vs_RespCEVsGen = ROOT.TProfile("H_GenPi_Pt_vs_RespCEVsGen","H_GenPi_Pt_vs_RespCEVsGen",50,0.,200.,0.8,1.2)
H_GenGM_Eta_vs_Resp = ROOT.TProfile("H_GenPi_Eta_vs_Resp","H_GenPi_Eta_vs_Resp",40,-4,4.,0.8,1.2)
H_GenGM_Eta_vs_RespE = ROOT.TProfile("H_GenPi_Eta_vs_RespE","H_GenPi_Eta_vs_RespE",40,-4.,4.,0.8,1.2)
H_GenGM_Eta_vs_RespCE = ROOT.TProfile("H_GenPi_Eta_vs_RespCE","H_GenPi_Eta_vs_RespCE",40,-4,4.,0.8,1.2)
H_GenGM_Eta_vs_RespCEVsGen = ROOT.TProfile("H_GenPi_Eta_vs_RespCEVsGen","H_GenPi_Eta_vs_RespCEVsGen",40,-4.,4.,0.8,1.2)

# load FWlite python libraries
#===============================
from DataFormats.FWLite import Handle, Events

trackingparticles, trackingparLabel = Handle("std::vector<TrackingParticle>"), "mix:MergedTrackTruth"
genparticles, genparLabel = Handle("std::vector<reco::GenParticle>"), "genParticles"
pfcluster, pfclusterLabel = Handle("std::vector<reco::PFCluster>"), "particleFlowClusterECAL"
pfclusterU, pfclusterULabel = Handle("std::vector<reco::PFCluster>"), "particleFlowClusterECALUncorrected"
pfclusterHLT, pfclusterHLTLabel = Handle("std::vector<reco::PFCluster>"), "hltParticleFlowClusterECALUnseeded"
pfclusterL1HLT, pfclusterL1HLTLabel = Handle("std::vector<reco::PFCluster>"), "hltParticleFlowClusterECALL1Seeded"

#caloparticles, caloparLabel = Handle("std::vector<CaloParticle>"), "mix:MergedCaloTruth"
#vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "inclusiveSecondaryVertices"
#pfcands, pfcandLabel = Handle("std::vector<reco::PFCandidate>"), "particleFlowEGamma"
#vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlinePrimaryVertices"
#pfcands, pfcandLabel = Handle("std::vector<reco::PFCandidate>"), "particleFlow"
#pfcluster, pfclusterLabel = Handle("std::vector<reco::PFCluster>"), "particleFlowClusterECAL"
#pfclusterU, pfclusterULabel = Handle("std::vector<reco::PFCluster>"), "particleFlowClusterECALUncorrected"
#pfclusterHLT, pfclusterHLTLabel = Handle("std::vector<reco::PFCluster>"), "hltParticleFlowClusterECALUnseeded"
#pfclusterL1HLT, pfclusterL1HLTLabel = Handle("std::vector<reco::PFCluster>"), "hltParticleFlowClusterECALL1Seeded"
#pfrechithf, pfrechithfLabel = Handle("std::vector<reco::PFRecHit>"), "particleFlowRecHitHF"
#pftracks, pftrackLabel = Handle("std::vector<reco::PFRecTrack>"), "pfTrack"

#gedPhotons, gedPhotonLabel = Handle("std::vector<reco::Photon>"), "gedPhotons"
#gedGsfElectrons, gedGsfElectronLabel = Handle("std::vector<reco::GsfElectron>"), "gedGsfElectrons"
#svcands, svcandLabel = Handle("std::vector<reco::VertexCompositePtrCandidate>"), "inclusiveCandidateSecondaryVertices"

#pfcandPtScore = Handle("edm::ValueMap<float>")
#verticesScore = Handle("edm::ValueMap<float>")

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
#=================================================
#================================================

#events = Events('file:EGM-Run3Winter24Reco-0p01to10-FlatPU0to80ZM.root')
events = Events('file:EGM-Run3Winter24Reco-%s-%s.root' %(PTrange, PUscenario))

# read list of default miniaod files in das
list_files = open('files_%s_%s.txt' %(PTrange, PUscenario), "r")
tmp = list_files.read() 
list = tmp.split("\n")
list.pop() #
fulllist = ['root://cmsxrootd.fnal.gov//'+x for x in list]
#fulllist = ['file:/cms/data/'+x for x in list]
print(fulllist)
res = fulllist[:20] # only 2 files
print(res)
if Option=="Default":
    events = Events(res)

# event loop starts    
for iev,event in enumerate(events):
    #if iev >= 10: break
    if iev%1000==1: print(iev)
    event.getByLabel(trackingparLabel, trackingparticles)
    event.getByLabel(genparLabel, genparticles)
    event.getByLabel(pfclusterLabel, pfcluster)
    event.getByLabel(pfclusterULabel, pfclusterU)
    if Option=="Custom":
        event.getByLabel(pfclusterHLTLabel, pfclusterHLT)
        event.getByLabel(pfclusterL1HLTLabel, pfclusterL1HLT)
    #event.getByLabel(puinfoLabel, puinfo)
    #event.getByLabel(vertexLabel, vertices)
    #event.getByLabel(vertexLabel, verticesScore)
    #event.getByLabel(pfcandLabel, pfcands)
    #event.getByLabel(gedPhotonLabel, gedPhotons)
    #event.getByLabel(gedGsfElectronLabel, gedGsfElectrons)
    #event.getByLabel(svcandLabel, svcands)
    #event.getByLabel(pfcandHcalDepthLabel,hcalDepthScore)
    #event.getByLabel(pfcandPtLabel,pfcandPtScore)
    #event.getByLabel(pfclusterLabel, pfcluster)
    #event.getByLabel(pfclusterULabel, pfclusterU)
    #event.getByLabel(pfclusterHLTLH_GenGM_Ptabel, pfclusterHLT)
    #event.getByLabel(pfclusterL1HLTLabel, pfclusterL1HLT)

    # PU
    #puinfo.product()
    #for i,pu in enumerate():  # loop over gen candidates
    #    print(pu.getTrueNumInteractions(),pu.getPU_NumInteractions())

    # Tracking particles
    for i,j in enumerate(trackingparticles.product()):  # loop over tracking particles

        goodPhoton=False
        promptPhoton=False
        if j.status()==1 and j.pdgId()==22: # and j.decayVertices().size()>0 and j.genParticles().size()>0:
            promptPhoton=True # this is a photon matached to a gen-level photon (prompt)
            if debug:
                print("TrackingParticle:  %3d: pt %5.1f eta %5.2f phi %5.2f status %4d " % ( i, j.pt(), j.eta(), j.phi(), j.status()))
                print(j.status(),j.pdgId(),j.decayVertices().at(0).position().rho(),j.decayVertices().at(0).position().z(),j.genParticles().at(0).pt())
            # check that this photon didn't decay before entering ECAL (EB case and EE case)
            if abs(j.decayVertices().at(0).position().z())<304.5 and j.decayVertices().at(0).position().rho()>123.8:
                goodPhoton=True # decay vertex is within EB
            if abs(j.decayVertices().at(0).position().z())>317.0:
                goodPhoton=True # decay vertex is within EE
            if goodPhoton:
                igen=j.genParticles().at(0)

        if goodPhoton==False:
            if promptPhoton:
                findEle=False
                #print(j.decayVertices().at(0).daughterTracks().size())
                for x in range (0, len(j.decayVertices().at(0).daughterTracks())):
                    if debug: print(j.decayVertices().at(0).daughterTracks().at(x).pdgId())
                    if abs(j.decayVertices().at(0).daughterTracks().at(x).pdgId())==11:
                        findEle=True
                    if abs(j.decayVertices().at(0).daughterTracks().at(x).pdgId())!=11:
                        if debug: print("Non e+e- decay of photon????",j.decayVertices().at(0).daughterTracks().at(x).pdgId()) 
                if findEle==False: print("prompt but converts before entering ECAL")
            continue

    # Gen particles ........................................
    #for i,igen in enumerate(genparticles.product()):  # loop over gen candidates

        if debug: print("GenParticle: iev %3d genpars %3d: pt %5.1f eta %5.2f phi %5.2f " % ( iev, i, igen.pt(), igen.eta(), igen.phi() ))

        H_GenGM_Pt.Fill(igen.pt())
        H_GenGM_Eta.Fill(igen.eta())

        nMatchOff=0
        nMatchHLT=0

        minR=1.0
        for i,pfc in enumerate(pfcluster.product()):  # loop over pf clusters
            if pfc.pt()>1.0 and deltaR(igen,pfc) <=0.1:
                if debug: print("Cluster   : iev %3d pfclus %3d: pt %5.1f eta %5.2f phi %5.2f layer %3d E %5.2f %5.2f %5.2f" % ( iev, i, pfc.pt(), pfc.eta(), pfc.phi(), pfc.layer(), pfc.correctedEnergy(), pfc.energy(), pfc.correctedEnergy()/pfc.energy() ))
                if deltaR(igen,pfc)<minR:
                    minR=deltaR(igen,pfc)
                    pfc_off=pfc
                    nMatchOff=nMatchOff+1

        # for i,pfc in enumerate(pfclusterU.product()):  # loop over pf clusters
        #     #pfc = j
        #     if pfc.pt()>1.0 and deltaR(igen,pfc) <=0.1:
        #         print("ClusterU  : iev %3d pfclus %3d: pt %5.1f eta %5.2f phi %5.2f layer %3d E %5.2f %5.2f %5.2f" % ( iev, i, pfc.pt(), pfc.eta(), pfc.phi(), pfc.layer(), pfc.correctedEnergy(), pfc.energy(), pfc.correctedEnergy()/pfc.energy() ))

        minR=1.0
        if Option=="Custom": # only when customized miniaod is used
            for i,pfc in enumerate(pfclusterHLT.product()):  # loop over pf clusters
                if pfc.pt()>1.0 and deltaR(igen,pfc) <=0.1:
                    if debug: print("ClusterHLT: iev %3d pfclus %3d: pt %5.1f eta %5.2f phi %5.2f layer %3d E %5.2f %5.2f %5.2f" % ( iev, i, pfc.pt(), pfc.eta(), pfc.phi(), pfc.layer(), pfc.correctedEnergy(), pfc.energy(), pfc.correctedEnergy()/pfc.energy() ))
                    if deltaR(igen,pfc)<minR:
                        minR=deltaR(igen,pfc)
                        pfc_hlt=pfc
                        nMatchHLT=nMatchHLT+1
                        
        # ===================================================

        if nMatchOff>0 and Option=="Default":
            H_GenGM_Pt_vs_RespCEVsGen.Fill(igen.pt(),pfc_off.correctedEnergy()/igen.energy())
            H_GenGM_Eta_vs_RespCEVsGen.Fill(igen.eta(),pfc_off.correctedEnergy()/igen.energy())
            
        if nMatchOff>0 and nMatchHLT>0 and Option=="Custom":
            H_GenGM_Pt_vs_Resp.Fill(igen.pt(),pfc_off.pt()/pfc_hlt.pt())
            H_GenGM_Pt_vs_RespE.Fill(igen.pt(),pfc_off.energy()/pfc_hlt.energy())
            H_GenGM_Pt_vs_RespCE.Fill(igen.pt(),pfc_off.correctedEnergy()/pfc_hlt.correctedEnergy())
            H_GenGM_Eta_vs_Resp.Fill(igen.eta(),pfc_off.pt()/pfc_hlt.pt())
            H_GenGM_Eta_vs_RespE.Fill(igen.eta(),pfc_off.energy()/pfc_hlt.energy())
            H_GenGM_Eta_vs_RespCE.Fill(igen.eta(),pfc_off.correctedEnergy()/pfc_hlt.correctedEnergy())

# Saving histograms in output root file
f = ROOT.TFile.Open('hist_%s_%s_%s.root' %(PTrange, PUscenario, Option),"RECREATE")
H_GenGM_Pt.Write()
H_GenGM_Eta.Write()
H_GenGM_Pt_vs_Resp.Write()
H_GenGM_Pt_vs_RespE.Write()
H_GenGM_Pt_vs_RespCE.Write()
H_GenGM_Pt_vs_RespCEVsGen.Write()
H_GenGM_Eta_vs_Resp.Write()
H_GenGM_Eta_vs_RespE.Write()
H_GenGM_Eta_vs_RespCE.Write()
H_GenGM_Eta_vs_RespCEVsGen.Write()
f.Write()
f.Close()
