import math

#
# Define a few useful functions
#

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

#
# import ROOT in batch mode
#
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
from ROOT import TF1, TF2, TH1, TH2, TH2F, TProfile, TAxis, TMath, TEllipse, TStyle, TFile, TColor, TSpectrum, TCanvas, TPad, TVirtualFitter, gStyle, TLorentzVector
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from ctypes import c_uint8

#
# load FWLite C++ libraries
#
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.AutoLibraryLoader.enable()

#
# Create histograms, etc.
#
'''
ROOT.gROOT.SetStyle('Plain') # white background
H_GenPi_Pt = ROOT.TH1F("H_GenPi_Pt","H_GenPi_Pt",50,0.,50.)
H_GenPi_Eta = ROOT.TH1F("H_GenPi_Eta","H_GenPi_Eta",100,-5.0,5.0)

H_GenPi_Eta_Sml_pT1 = ROOT.TH1F("Response_Sml_pT1" , "0.0 < pT < 2.0",40,0.0,2.0)
H_GenPi_Eta_Sml_pT2 = ROOT.TH1F("Response_Sml_pT2" , "2.0 < pT < 4.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_pT3 = ROOT.TH1F("Response_Sml_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_pT4 = ROOT.TH1F("Response_Sml_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_pT5 = ROOT.TH1F("Response_Sml_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_pT1 = ROOT.TH1F("Response_Med_pT1" , "0.0 < pT < 2.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_pT2 = ROOT.TH1F("Response_Med_pT2" , "2.0 < pT < 4.0",40,0.0,2.0)
H_GenPi_Eta_Med_pT3 = ROOT.TH1F("Response_Med_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_pT4 = ROOT.TH1F("Response_Med_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_pT5 = ROOT.TH1F("Response_Med_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_pT1 = ROOT.TH1F("Response_Lrg_pT1" , "0.0 < pT < 2.0",40,0.0,2.0)
H_GenPi_Eta_Lrg_pT2 = ROOT.TH1F("Response_Lrg_pT2" , "2.0 < pT < 4.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_pT3 = ROOT.TH1F("Response_Lrg_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_pT4 = ROOT.TH1F("Response_Lrg_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_pT5 = ROOT.TH1F("Response_Lrg_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)

H_GenPi_Eta_Sml_tpT1 = ROOT.TH1F("Tracking_Response_Sml_pT1" , "0.0 < pT < 2.0",40,0.0,2.0)
H_GenPi_Eta_Sml_tpT2 = ROOT.TH1F("Tracking_Response_Sml_pT2" , "2.0 < pT < 4.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_tpT3 = ROOT.TH1F("Tracking_Response_Sml_pT3" , "4.0 < pT < 6.0", 40, 0.0, 2.0,)
H_GenPi_Eta_Sml_tpT4 = ROOT.TH1F("Tracking_Response_Sml_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_tpT5 = ROOT.TH1F("Tracking_Response_Sml_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_tpT1 = ROOT.TH1F("Tracking_Response_Med_pT1" , "0.0 < pT < 2.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_tpT2 = ROOT.TH1F("Tracking_Response_Med_pT2" , "2.0 < pT < 4.0",40,0.0,2.0)
H_GenPi_Eta_Med_tpT3 = ROOT.TH1F("Tracking_Response_Med_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_tpT4 = ROOT.TH1F("Tracking_Response_Med_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_tpT5 = ROOT.TH1F("Tracking_Response_Med_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_tpT1 = ROOT.TH1F("Tracking_Response_Lrg_pT1" , "0.0 < pT < 2.0",40,0.0,2.0)
H_GenPi_Eta_Lrg_tpT2 = ROOT.TH1F("Tracking_Response_Lrg_pT2" , "2.0 < pT < 4.0",40, 0.0, 2.0)
H_GenPi_Eta_Lrg_tpT3 = ROOT.TH1F("Tracking_Response_Lrg_pT3" , "4.0 < pT < 6.0",40, 0.0, 2.0)
H_GenPi_Eta_Lrg_tpT4 = ROOT.TH1F("Tracking_Response_Lrg_pT4" , "6.0 < pT < 8.0",40, 0.0, 2.0)
H_GenPi_Eta_Lrg_tpT5 = ROOT.TH1F("Tracking_Response_Lrg_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)

H_GenPi_Eta_Sml_LorentzT1 = ROOT.TH1F("Lorentz_Response_Sml_pT1" , "0.0 < pT < 2.0",40,0.0,2.0)
H_GenPi_Eta_Sml_LorentzT2 = ROOT.TH1F("Lorentz_Response_Sml_pT2" , "2.0 < pT < 4.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_LorentzT3 = ROOT.TH1F("Lorentz_Response_Sml_pT3" , "4.0 < pT < 6.0", 40, 0.0, 2.0,)
 vfhH_GenPi_Eta_Sml_LorentzT4 = ROOT.TH1F("Lorentz_Response_Sml_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_LorentzT5 = ROOT.TH1F("Lorentz_Response_Sml_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_LorentzT1 = ROOT.TH1F("Lorentz_Response_Med_pT1" , "0.0 < pT < 2.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_LorentzT2 = ROOT.TH1F("Lorentz_Response_Med_pT2" , "2.0 < pT < 4.0",40,0.0,2.0)
H_GenPi_Eta_Med_LorentzT3 = ROOT.TH1F("Lorentz_Response_Med_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_LorentzT4 = ROOT.TH1F("Lorentz_Response_Med_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_LorentzT5 = ROOT.TH1F("Lorentz_Response_Med_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_LorentzT1 = ROOT.TH1F("Lorentz_Response_Lrg_pT1" , "0.0 < pT < 2.0",40,0.0,2.0)
H_GenPi_Eta_Lrg_LorentzT2 = ROOT.TH1F("Lorentz_Response_Lrg_pT2" , "2.0 < pT < 4.0",40, 0.0, 2.0)
H_GenPi_Eta_Lrg_LorentzT3 = ROOT.TH1F("Lorentz_Response_Lrg_pT3" , "4.0 < pT < 6.0",40, 0.0, 2.0)
H_GenPi_Eta_Lrg_LorentzT4 = ROOT.TH1F("Lorentz_Response_Lrg_pT4" , "6.0 < pT < 8.0",40, 0.0, 2.0)
H_GenPi_Eta_Lrg_LorentzT5 = ROOT.TH1F("Lorentz_Response_Lrg_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)

H_GenPi_Eta_Sml_TrkMerpT1 = ROOT.TH1F("TrkMegResponse_Sml_pT1" , "0.0 < pT < 2.0",40,0.0,2.0)
H_GenPi_Eta_Sml_TrkMerpT2 = ROOT.TH1F("TrkMegResponse_Sml_pT2" , "2.0 < pT < 4.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_TrkMerpT3 = ROOT.TH1F("TrkMegResponse_Sml_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_TrkMerpT4 = ROOT.TH1F("TrkMegResponse_Sml_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_TrkMerpT5 = ROOT.TH1F("TrkMegResponse_Sml_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_TrkMegpT1 = ROOT.TH1F("TrkMegResponse_Med_pT1" , "0.0 < pT < 2.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_TrkMegpT2 = ROOT.TH1F("TrkMegResponse_Med_pT2" , "2.0 < pT < 4.0",40,0.0,2.0)
H_GenPi_Eta_Med_TrkMegpT3 = ROOT.TH1F("TrkMegResponse_Med_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_TrkMegpT4 = ROOT.TH1F("TrkMegResponse_Med_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_TrkMegpT5 = ROOT.TH1F("TrkMegResponse_Med_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_TrkMegpT1 = ROOT.TH1F("TrkMegResponse_Lrg_pT1" , "0.0 < pT < 2.0",40,0.0,2.0)
H_GenPi_Eta_Lrg_TrkMegpT2 = ROOT.TH1F("TrkMegResponse_Lrg_pT2" , "2.0 < pT < 4.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_TrkMegpT3 = ROOT.TH1F("TrkMegResponse_Lrg_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_TrkMegpT4 = ROOT.TH1F("TrkMegResponse_Lrg_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_TrkMegpT5 = ROOT.TH1F("TrkMegResponse_Lrg_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)

H_GenPi_Eta_Sml_TrkTrkpT1 = ROOT.TH1F("Response_Sml_pT1" , "0.0 < pT < 2.0",40,0.0,2.0)
H_GenPi_Eta_Sml_TrkTrkpT2 = ROOT.TH1F("Response_Sml_pT2" , "2.0 < pT < 4.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_TrkTrkpT3 = ROOT.TH1F("Response_Sml_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_TrkTrkpT4 = ROOT.TH1F("Response_Sml_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Sml_TrkTrkpT5 = ROOT.TH1F("Response_Sml_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_TrkTrkpT1 = ROOT.TH1F("Response_Med_pT1" , "0.0 < pT < 2.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_TrkTrkpT2 = ROOT.TH1F("Response_Med_pT2" , "2.0 < pT < 4.0",40,0.0,2.0)
H_GenPi_Eta_Med_TrkTrkpT3 = ROOT.TH1F("Response_Med_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_TrkTrkpT4 = ROOT.TH1F("Response_Med_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Med_TrkTrkpT5 = ROOT.TH1F("Response_Med_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_TrkTrkpT1 = ROOT.TH1F("Response_Lrg_pT1" , "0.0 < pT < 2.0",40,0.0,2.0)
H_GenPi_Eta_Lrg_TrkTrkpT2 = ROOT.TH1F("Response_Lrg_pT2" , "2.0 < pT < 4.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_TrkTrkpT3 = ROOT.TH1F("Response_Lrg_pT3" , "4.0 < pT < 6.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_TrkTrkpT4 = ROOT.TH1F("Response_Lrg_pT4" , "6.0 < pT < 8.0",40 , 0.0, 2.0)
H_GenPi_Eta_Lrg_TrkTrkpT5 = ROOT.TH1F("Response_Lrg_pT5" , "8.0 < pT < 10.0",40 , 0.0, 2.0)

H_GenPi_Resp = ROOT.TH1F("H_GenPi_Resp", "Response", 100, -3.0, 3.0)
'''
#
# load FWlite python libraries
#
from DataFormats.FWLite import Handle, Events

genparticles, genparLabel = Handle("std::vector<reco::GenParticle>"), "genParticles"
#caloparticles, caloparLabel = Handle("std::vector<CaloParticle>"), "mix:MergedCaloTruth"
#vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "inclusiveSecondaryVertices"
pfcands, pfcandLabel = Handle("std::vector<reco::PFCandidate>"), "particleFlowEGamma"
#vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlinePrimaryVertices"
pfcands, pfcandLabel = Handle("std::vector<reco::PFCandidate>"), "particleFlow"
#pfclusterhf, pfclusterhfLabel = Handle("std::vector<reco::PFCluster>"), "particleFlowClusterHF"
#pfrechithf, pfrechithfLabel = Handle("std::vector<reco::PFRecHit>"), "particleFlowRecHitHF"
pftracks, pftrackLabel = Handle("std::vector<reco::PFRecTrack>"), "pfTrack"

trackstersEM, trackstersEMLabel = Handle ("std::vector<ticl::Trackster>"), "ticlTrackstersEM"
trackstersHad, trackstersHadLabel = Handle ("std::vector<ticl::Trackster>"), "ticlTrackstersHAD"
trackstersTrk, trackstersTrkLabel = Handle ("std::vector<ticl::Trackster>"), "ticlTrackstersTrk"
trackstersTrkEm, trackstersTrkEmLabel = Handle ("std::vector<ticl::Trackster>"), "ticlTrackstersTrkEM"
trackstersMerge, trackstersMergeLabel = Handle ("std::vector<ticl::Trackster>"), "ticlTrackstersMerge"
ticlcands, ticlCandLabel = Handle ("std::vector<TICLCandidate>"), "ticlTrackstersMerge"

gedPhotons, gedPhotonLabel = Handle("std::vector<reco::Photon>"), "gedPhotons"
gedGsfElectrons, gedGsfElectronLabel = Handle("std::vector<reco::GsfElectron>"), "gedGsfElectrons"
#svcands, svcandLabel = Handle("std::vector<reco::VertexCompositePtrCandidate>"), "inclusiveCandidateSecondaryVertices"

#pfcandPtScore = Handle("edm::ValueMap<float>")
#verticesScore = Handle("edm::ValueMap<float>")

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
#events = Events('file:../inputFiles/F6F192B6-19E8-1245-913A-545FF4D712A1.root')
#events = Events('file:/afs/cern.ch/user/h/hatake/work/public/ForTim/RelValSinglePiFlatPt0p7To10_CMSSW_11_3_0_pre5_GEN-SIM-RECO.root')
events = Events('file:/afs/cern.ch/user/h/hatake/work/public/ForTim/RelValSinglePiFlatPt0p7To10_CMSSW_11_3_0_pre5_GEN-SIM-RECO.root')
#events = Events('file:step3_20evt_org.root')
#events = Events('file:step3_rerereco_all.root')
#events = Events('file:step3.root')
#events = Events('file:step3_withparticleFlowTmpBarrel.root')

#
# Event loop starts
#
for iev,event in enumerate(events):
    #if iev >= 10: break
    #event.getByLabel(vertexLabel, vertices)
    #event.getByLabel(vertexLabel, verticesScore)
    event.getByLabel(genparLabel, genparticles)
    event.getByLabel(pfcandLabel, pfcands)
    event.getByLabel(gedPhotonLabel, gedPhotons)
    event.getByLabel(gedGsfElectronLabel, gedGsfElectrons)
    event.getByLabel(trackstersMergeLabel, trackstersMerge)
    event.getByLabel(ticlCandLabel, ticlcands)
    #event.getByLabel(svcandLabel, svcands)
    #event.getByLabel(pfcandHcalDepthLabel,hcalDepthScore)
    #event.getByLabel(pfcandPtLabel,pfcandPtScore)

    print "******************************************************************************** ievent: ",iev

    ncount=0
    for k,kticl in enumerate(ticlcands.product()):
        ncount = ncount+1
        #print "k,kticl",k,kticl
    print ncount

    ncount=0
    for k,kticl in enumerate(trackstersMerge.product()):
        ncount = ncount+1
        #print "k,kticl",k,kticl
    print ncount

    #
    # Loop pver gen particles
    #
    for i,igen in enumerate(genparticles.product()):
        #print "i,igen",i,igen
        print("** gen pt, eta, phi %7.2f %7.2f %7.2f" % \
            (igen.pt(), igen.eta(), igen.phi()))

        #
        # Find best-matched pf candidate
        #
        mindR=100
        jpfMindR=-1
        for j,jpf in enumerate(pfcands.product()):
            #print "j,jpf",j,jpf            
            if (jpf.charge()>0 or jpf.charge()<0):
            #if (jpf.charge()>0 or jpf.charge()<0) \
            #   and (abs(jpf.pdgId())==211): # without pdgId check, get reference error.
                if deltaR(igen,jpf) <=0.1:
                    if jpf.pt()>0.1:
                        if deltaR(igen,jpf) < mindR:
                            mindR = deltaR(igen,jpf)
                            jpfMindR = j
        #
        # Analyze best matched PF candidate
        #
        if jpfMindR >= 0: # matched candidate is found
            jpf = pfcands.product()[jpfMindR]
            trkresp = 0.
            if jpf.gsfTrackRef().isNonnull():
                print "gsftrack. skip."
            else: 
                jtrk = jpf.trackRef()
                trkresp = jtrk.pt()/igen.pt()                
            print(" jpf,mindR,   pf   pt/eta/phi/pdgId, resp: %6d %7.4f, %7.2f %7.2f %7.2f %8d, %7.3f %7.3f" % \
                  (jpfMindR,mindR, \
                   jpf.pt(), jpf.eta(), jpf.phi(), jpf.pdgId(), \
                   jpf.pt()/igen.pt(), \
                   trkresp))

        #
        # Find best-matched ticlCandidate
        #
        mindR=100
        kticlMindR=-1
        for k,kticl in enumerate(ticlcands.product()):
            #print "k,kticl",k,kticl
            if kticl.charge()>0 or kticl.charge()<0:
                if deltaR(igen,kticl) <=0.1:
                    if kticl.pt()>0.1:
                        if deltaR(igen,kticl) < mindR:
                            mindR = deltaR(igen,kticl)
                            kticlMindR =k
        #
        # Analyze best matched ticlCandidate/trackstersMerge
        #
        if kticlMindR >= 0: # matched candidate is found
            #
            # first check ticl candidate
            kticl = ticlcands.product()[kticlMindR]
            ktrk = kticl.trackPtr()
            trkresp = ktrk.pt()/igen.pt()                
            print(" kticl,mindR, ticl pt/eta/phi/pdgId, resp: %6d %7.4f, %7.2f %7.2f %7.2f %8d, %7.3f %7.3f" % \
                  (kticlMindR,mindR, \
                   kticl.pt(), kticl.eta(), kticl.phi(), kticl.pdgId(), \
                   kticl.pt()/igen.pt(),
                   trkresp))
            #
            # now also check tracksters
            raw_energy=0.
            regressed_energy=0.
            for ktrkstr in kticl.tracksters():
                raw_energy += ktrkstr.raw_energy()
                regressed_energy += ktrkstr.regressed_energy()
            ktrkstr = trackstersMerge.product()[kticlMindR]
            mpion = 0.13957
            momentum = 0.
            momentum_regressed = 0.
            if raw_energy > 0:
                momentum = math.sqrt(raw_energy*raw_energy-mpion*mpion)
                momentum_regressed = momentum
            if regressed_energy > 0:
                momentum_regressed = math.sqrt(regressed_energy*regressed_energy-mpion*mpion)
            trkstr_p4 = TLorentzVector()
            trkstr_p4_regressed = TLorentzVector()
            trkstr_p4.SetPxPyPzE(momentum*ktrk.momentum().unit().x(),\
                                 momentum*ktrk.momentum().unit().y(),\
                                 momentum*ktrk.momentum().unit().z(),\
                                 raw_energy)
            trkstr_p4_regressed.SetPxPyPzE(momentum_regressed*ktrk.momentum().unit().x(),\
                                 momentum_regressed*ktrk.momentum().unit().y(),\
                                 momentum_regressed*ktrk.momentum().unit().z(),\
                                 regressed_energy)
            print(" ktrkstr raw_energy, regressed_energy, pt (raw, regressed): %7.2f %7.2f, %7.2f %7.2f" % \
                  (ktrkstr.raw_energy(), ktrkstr.regressed_energy(), \
                   trkstr_p4.Pt(), trkstr_p4_regressed.Pt()))

    # loop over gen particles first
    '''for i,igen in enumerate(genparticles.product()):
        mindR=100
        igenMindR=-1

        # loop over pf candidates, first loop to find mindR pair

        jsum = TLorentzVector()
        jpfTLV = TLorentzVector()
        jpfTLV.SetPtEtaPhiM(0, 0, 0, 0)
        for j,jpf in enumerate(pfcands.product()):
           # if jpf.charge()>0 or jpf.charge()<0:
                if deltaR(igen,jpf) <=0.2:
                    if jpf.pt()>0.1:
                        print 'pfcands: ', jpf.pt(), 'id ',jpf.pdgId()
                        jpfTLV.SetPtEtaPhiM(jpf.pt(),jpf.eta(),jpf.phi(),jpf.mass())
                        jsum+= jpfTLV
                        if deltaR(igen,jpf) < mindR:
                            mindR = deltaR(igen,jpf)
                            igenMindR =j
        print 'gen pt',igen.pt()
        print 'TLV pt', jsum.Pt()

        if abs(igen.eta()) > 1.3 and abs(igen.eta()) < 2.1:
            if igen.pt() < 2.0:
                H_GenPi_Eta_Sml_LorentzT1.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 2.0 and igen.pt() <4.0:
                H_GenPi_Eta_Sml_LorentzT2.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 4.0 and igen.pt() <6.0:
                H_GenPi_Eta_Sml_LorentzT3.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 6.0 and igen.pt() <8.0:
                H_GenPi_Eta_Sml_LorentzT4.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 8.0 and igen.pt() <10.0:
                H_GenPi_Eta_Sml_LorentzT5.Fill(jsum.Pt()/igen.pt())

        if abs(igen.eta()) > 2.1 and abs(igen.eta()) < 2.5:
            if igen.pt() < 2.0:
                H_GenPi_Eta_Med_LorentzT1.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 2.0 and igen.pt() <4.0:
                H_GenPi_Eta_Med_LorentzT2.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 4.0 and igen.pt() <6.0:
                H_GenPi_Eta_Med_LorentzT3.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 6.0 and igen.pt() <8.0:
                H_GenPi_Eta_Med_LorentzT4.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 8.0 and igen.pt() <10.0:
                H_GenPi_Eta_Med_LorentzT5.Fill(jsum.Pt()/igen.pt())

        if abs(igen.eta()) > 2.5 and abs(igen.eta()) < 3.0:
            if igen.pt() < 2.0:
                H_GenPi_Eta_Lrg_LorentzT1.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 2.0 and igen.pt() <4.0:
                H_GenPi_Eta_Lrg_LorentzT2.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 4.0 and igen.pt() <6.0:
                H_GenPi_Eta_Lrg_LorentzT3.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 6.0 and igen.pt() <8.0:
                H_GenPi_Eta_Lrg_LorentzT4.Fill(jsum.Pt()/igen.pt())
            if igen.pt() > 8.0 and igen.pt() <10.0:
                H_GenPi_Eta_Lrg_LorentzT5.Fill(jsum.Pt()/igen.pt())


    # loop over pf candidates, 2nd loop for further analyses

        for j,jpf in enumerate(pfcands.product()):
            if jpf.charge()>0 or jpf.charge()<0:
                if deltaR(igen,jpf) <=0.2:
                    if deltaR(igen,jpf) < mindR:
                        mindR = deltaR(igen,jpf)
                        igenMindR =j

                if j == igenMindR:
                    if abs(igen.pt()) > 0.1:
                        H_GenPi_Resp.Fill(jpf.pt()/igen.pt())
                        if abs(igen.eta()) > 1.3 and abs(igen.eta()) < 2.1:
                            if igen.pt() < 2.0:
                                H_GenPi_Eta_Sml_pT1.Fill(jpf.pt()/igen.pt())
                                H_GenPi_Eta_Sml_tpT1.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 2.0 and igen.pt() < 4.0:
                                H_GenPi_Eta_Sml_pT2.Fill(jpf.pt()/igen.pt())
                                H_GenPi_Eta_Sml_tpT2.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 4.0 and igen.pt() <6.0:
                                H_GenPi_Eta_Sml_pT3.Fill((jpf.pt()+(-0.8123670471028798*jpf.pt() +  4.676169889317416))/igen.pt())
                                H_GenPi_Eta_Sml_tpT3.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 6.0 and igen.pt() < 8.0:
                                H_GenPi_Eta_Sml_pT4.Fill((jpf.pt()+ (-0.847206245687777*jpf.pt() + 6.6082879097191105))/igen.pt())
                                H_GenPi_Eta_Sml_tpT4.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 8.0 and igen.pt() < 10.0:
                                H_GenPi_Eta_Sml_pT5.Fill((jpf.pt()+(-0.9387780712767415*jpf.pt() + 8.831262389664225))/igen.pt())
                                H_GenPi_Eta_Sml_tpT5.Fill(jpf.trackRef().pt()/igen.pt())

                        if abs(igen.eta()) > 2.1 and abs(igen.eta()) < 2.5:
                            if igen.pt() < 2.0:
                                H_GenPi_Eta_Med_pT1.Fill(jpf.pt()/igen.pt())
                                H_GenPi_Eta_Med_tpT1.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 2.0 and igen.pt() <4.0:
                                H_GenPi_Eta_Med_pT2.Fill(jpf.pt()/igen.pt())
                                H_GenPi_Eta_Med_tpT2.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 4.0 and igen.pt() <6.0:
                                H_GenPi_Eta_Med_pT3.Fill((jpf.pt()+(-0.8123670471028798*jpf.pt() +  4.676169889317416))/igen.pt())
                                H_GenPi_Eta_Med_tpT3.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 6.0 and igen.pt() <8.0:
                                H_GenPi_Eta_Med_pT4.Fill((jpf.pt()+ (-0.847206245687777*jpf.pt() + 6.6082879097191105))/igen.pt())
                                H_GenPi_Eta_Med_tpT4.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 8.0 and igen.pt() <10.0:
                                H_GenPi_Eta_Med_pT5.Fill((jpf.pt()+(-0.9387780712767415*jpf.pt() + 8.831262389664225))/igen.pt())
                                H_GenPi_Eta_Med_tpT5.Fill(jpf.trackRef().pt()/igen.pt())


                        if abs(igen.eta()) > 2.5 and abs(igen.eta()) < 3.0:
                            if igen.pt() < 2.0:
                                H_GenPi_Eta_Lrg_pT1.Fill(jpf.pt()/igen.pt())
                                H_GenPi_Eta_Lrg_tpT1.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 2.0 and igen.pt() <4.0:
                                H_GenPi_Eta_Lrg_pT2.Fill(jpf.pt()/igen.pt())
                                H_GenPi_Eta_Lrg_tpT2.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 4.0 and igen.pt() <6.0:
                                H_GenPi_Eta_Lrg_pT3.Fill((jpf.pt()+(-0.8123670471028798*jpf.pt() +  4.676169889317416))/igen.pt())
                                H_GenPi_Eta_Lrg_tpT3.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 6.0 and igen.pt() <8.0:
                                H_GenPi_Eta_Lrg_pT4.Fill((jpf.pt()+ (-0.847206245687777*jpf.pt() + 6.6082879097191105))/igen.pt())
                                H_GenPi_Eta_Lrg_tpT4.Fill(jpf.trackRef().pt()/igen.pt())
                            if igen.pt() > 8.0 and igen.pt() <10.0:
                                H_GenPi_Eta_Lrg_pT5.Fill((jpf.pt()+(-0.9387780712767415*jpf.pt() + 8.831262389664225))/igen.pt())
                                H_GenPi_Eta_Lrg_tpT5.Fill(jpf.trackRef().pt()/igen.pt())

    #print "match: ",igen.pt(),sum_pt

    #photons
    npf=0
    for i,j in enumerate(gedPhotons.product()):
       npf=npf+1
    #   print "gedPhotons: pt %17.13f eta %18.14f pdgId %5d " % ( j.pt(), j.eta(), j.pdgId())
    #   print "nphoton: ",npf

    # electrons
    npf=0
    for i,j in enumerate(gedGsfElectrons.product()):
       npf=npf+1
    #       print "gedGsfElectrons: pt %17.13f eta %18.14f pdgId %5d " % ( j.pt(), j.eta(), j.pdgId())
    #       print "nelectron: ",npf

    # taus
    
    # jets

    # met

#
# Event-loop is over
#

f = ROOT.TFile.Open("myfile.root","RECREATE")

H_GenPi_Pt.Write()
H_GenPi_Eta.Write()

H_GenPi_Eta_Sml_pT1.Write()
H_GenPi_Eta_Sml_pT2.Write()
H_GenPi_Eta_Sml_pT3.Write()
H_GenPi_Eta_Sml_pT4.Write()
H_GenPi_Eta_Sml_pT5.Write()
H_GenPi_Eta_Med_pT1.Write()
H_GenPi_Eta_Med_pT2.Write()
H_GenPi_Eta_Med_pT3.Write()
H_GenPi_Eta_Med_pT4.Write()
H_GenPi_Eta_Med_pT5.Write()
H_GenPi_Eta_Lrg_pT1.Write()
H_GenPi_Eta_Lrg_pT2.Write()
H_GenPi_Eta_Lrg_pT3.Write()
H_GenPi_Eta_Lrg_pT4.Write()
H_GenPi_Eta_Lrg_pT5.Write()

H_GenPi_Eta_Sml_tpT1.Write()
H_GenPi_Eta_Sml_tpT2.Write()
H_GenPi_Eta_Sml_tpT3.Write()
H_GenPi_Eta_Sml_tpT4.Write()
H_GenPi_Eta_Sml_tpT5.Write()
H_GenPi_Eta_Med_tpT1.Write()
H_GenPi_Eta_Med_tpT2.Write()
H_GenPi_Eta_Med_tpT3.Write()
H_GenPi_Eta_Med_tpT4.Write()
H_GenPi_Eta_Med_tpT5.Write()
H_GenPi_Eta_Lrg_tpT1.Write()
H_GenPi_Eta_Lrg_tpT2.Write()
H_GenPi_Eta_Lrg_tpT3.Write()
H_GenPi_Eta_Lrg_tpT4.Write()
H_GenPi_Eta_Lrg_tpT5.Write()

H_GenPi_Eta_Sml_LorentzT1.Write()
H_GenPi_Eta_Sml_LorentzT2.Write()
H_GenPi_Eta_Sml_LorentzT3.Write()
H_GenPi_Eta_Sml_LorentzT4.Write()
H_GenPi_Eta_Sml_LorentzT5.Write()
H_GenPi_Eta_Med_LorentzT1.Write()
H_GenPi_Eta_Med_LorentzT2.Write()
H_GenPi_Eta_Med_LorentzT3.Write()
H_GenPi_Eta_Med_LorentzT4.Write()
H_GenPi_Eta_Med_LorentzT5.Write()
H_GenPi_Eta_Lrg_LorentzT1.Write()
H_GenPi_Eta_Lrg_LorentzT2.Write()
H_GenPi_Eta_Lrg_LorentzT3.Write()
H_GenPi_Eta_Lrg_LorentzT4.Write()
H_GenPi_Eta_Lrg_LorentzT5.Write()

H_GenPi_Resp.Write()


#ROOT.TFile.Close(f)
f.Write()
f.Close()
'''
