#
# Usage example:
# python plot_reco_pfcands_ForBU.py PFH_CaloEnergy_FlatQCD_PU_2018 
# python plot_reco_pfcands_ForBU.py PFH_CaloEnergy_FlatQCD_NoPU_2018
# python plot_reco_pfcands_ForBU.py PFH_CaloEnergy_JetHT_2016H-10_6_0_MatrixTest >& PFH_CaloEnergy_JetHT_2016H-10_6_0_MatrixTest_v01.log &
# python plot_reco_pfcands_ForBU.py PFH_CaloEnergy_JetHT_2018D-10_6_0_MatrixTest >& PFH_CaloEnergy_JetHT_2018D-10_6_0_MatrixTest_v01.log &
#
# import ROOT in batch mode
import sys
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
ROOT.AutoLibraryLoader.enable()

# Create histograms, etc.
ROOT.gROOT.SetStyle('Plain') # white background
H_NPV = ROOT.TH1F ("NPV","NPV",101,-0.5,100.5)
H_PFH_HCAL = ROOT.TH1F ("PFH_HCAL", "Neutral hadron HCAL energy", 100, 0, 20)
H_PFH_CorECAL = ROOT.TH1F ("PFH_CorECAL", "Neutral hadron corrected ECAL energy", 100, 0, 20)
H_PFH_RawECAL = ROOT.TH1F ("PFH_RawECAL", "Neutral hadron raw ECAL energy", 100, 0, 20)
H_PFH_CorECAL_zoom = ROOT.TH1F ("PFH_CorECAL_zoom", "Neutral hadron corrected ECAL energy zoom", 60, -3., 3.)
H_PFH_RawECAL_zoom = ROOT.TH1F ("PFH_RawECAL_zoom", "Neutral hadron raw ECAL energy zoom", 60, -3., 3.)
H_PFH_withTrack = ROOT.TProfile ("H_PFH_withTrack", "H_PFH_withTrack", 80, -4.0, 4.0, -1.0, 2.0)
H_PFH_withoutTrack = ROOT.TProfile ("H_PFH_withoutTrack", "H_PFH_withoutTrack", 80, -4.0, 4.0, -1.0, 2.0)
H_PFH_HcalTrackdR = ROOT.TH1F ("H_PFH_HcalTrackdR", "H_PFH_HcalTrackdR", 80, 0., 4.0)
H_PFH_HcalClusterSize = ROOT.TH1F ("H_PFH_HcalClusterSize", "H_PFH_HcalClusterSize", 400, 0., 400.)
H_PFH_HcalClusterSizeVsEnergy = ROOT.TH2F ("H_PFH_HcalClusterSizeVsEnergy", "H_PFH_HcalClusterSizeVsEnergy", 50, 0., 100., 50, 0., 400.)

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

#muons, muonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"
#electrons, electronLabel = Handle("std::vector<pat::Electron>"), "slimmedElectrons"
#photons, photonLabel = Handle("std::vector<pat::Photon>"), "slimmedPhotons"
#taus, tauLabel = Handle("std::vector<pat::Tau>"), "slimmedTaus"
#jets, jetLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
#fatjets, fatjetLabel = Handle("std::vector<pat::Jet>"), "slimmedJetsAK8"
#mets, metLabel = Handle("std::vector<pat::MET>"), "slimmedMETs"
#vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlineSlimmedPrimaryVertices"
#offlinePrimaryVertices
vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlinePrimaryVertices"
pfcands, pfcandLabel = Handle("std::vector<reco::PFCandidate>"), "particleFlow"

pfcandPtScore = Handle("edm::ValueMap<float>")
verticesScore = Handle("edm::ValueMap<float>")

print len(sys.argv)

if len(sys.argv)>1:
    output=sys.argv[1]
else:
    output="PFH_CaloEnergy_JetHT_2018D-10_6_0_MatrixTest"
        
# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
if output == "PFH_CaloEnergy_FlatQCD_NoPU_2018":
    events = Events('root://cmsxrootd.fnal.gov//store/relval/CMSSW_10_6_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/106X_upgrade2018_realistic_v4-v1/10000/6C4012B3-14CA-8844-8301-D4B253D75644.root')
elif output == "PFH_CaloEnergy_FlatQCD_PU_2018":
    events = Events('root://cmsxrootd.fnal.gov//store/relval/CMSSW_10_6_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/PU25ns_106X_upgrade2018_realistic_v4-v1/10000/FDEFFEB4-107A-8B44-8E52-F84F8B3E158A.root')
elif output == "PFH_CaloEnergy_SingleMuon_2016C":
    events = Events('root://cmsxrootd.fnal.gov//store/data/Run2016C/SingleMuon/AOD/07Aug17-v1/110001/44B1F7D7-3F80-E711-AE45-001E67E6F4A9.root')
elif output == "PFH_CaloEnergy_SingleMuon_2017C":
    events = Events('root://cmsxrootd.fnal.gov//store/data/Run2017C/SingleMuon/AOD/17Nov2017-v1/70000/A073FB9C-6CDA-E711-9551-02163E01A465.root')
elif output == "PFH_CaloEnergy_SingleMuon_2018C":
    events = Events('root://cmsxrootd.fnal.gov//store/data/Run2018C/SingleMuon/AOD/17Sep2018-v1/00001/56DF33B4-0A6E-2747-99BC-3FD669B097FF.root')
elif output == "PFH_CaloEnergy_SingleMuon_2016H-10_6_0_Relval":
    events = Events('root://cmsxrootd.fnal.gov//store/relval/CMSSW_10_6_0/SingleMuon/RECO/106X_dataRun2_v10_RelVal_2016H-v1/10000/BC6B00A1-72AF-FF44-861E-87DD48C2F41C.root')
elif output == "PFH_CaloEnergy_SingleMuon_2017F-10_6_0_Relval":
    events = Events('root://cmsxrootd.fnal.gov//store/relval/CMSSW_10_6_0/SingleMuon/RECO/106X_dataRun2_v10_RelVal_2017F-v1/10000/515ED924-A207-2146-9C0F-6F99A20F2DB0.root')
elif output == "PFH_CaloEnergy_SingleMuon_2018C-10_6_0_Relval":
    events = Events('root://cmsxrootd.fnal.gov//store/relval/CMSSW_10_6_0/SingleMuon/RECO/106X_dataRun2_v10_resub_RelVal_2018C-v1/10000/66466716-600E-B540-BB79-6FC151D62972.root')
elif output == "PFH_CaloEnergy_JetHT_2018D-10_6_0_MatrixTest":
    events = Events('/home/hatake/ana_cms/PF/sl7/CMSSW_11_0_0_pre3_caloFracFix/src/136.888_RunJetHT2018D+RunJetHT2018D+HLTDR2_2018+RECODR2_2018reHLT_skimJetHT_Prompt+HARVEST2018_Prompt/backup3/step3.root')
elif output == "PFH_CaloEnergy_JetHT_2016H-10_6_0_MatrixTest":
    events = Events('/home/hatake/ana_cms/PF/sl7/CMSSW_11_0_0_pre3_caloFracFix/src/136.772_RunJetHT2016H+RunJetHT2016H+HLTDR2_2016+RECODR2_2016reHLT_skimJetHT_Prompt+HARVESTDR2/backup3/step3.root')

#
# Dataset memo
#
# /SingleMuon/CMSSW_10_6_0-106X_dataRun2_PromptLike_v6_RelVal_2018D-v1/RECO
# /SingleMuon/CMSSW_10_6_0-106X_dataRun2_v10_RelVal_2018C-v1/RECO
# /SingleMuon/CMSSW_10_6_0-106X_dataRun2_v10_resub_RelVal_2018C-v1/RECO
#  /store/relval/CMSSW_10_6_0/SingleMuon/RECO/106X_dataRun2_v10_resub_RelVal_2018C-v1/10000/66466716-600E-B540-BB79-6FC151D62972.root
# /SingleMuon/CMSSW_10_6_0-106X_dataRun2_v10_RelVal_2017F-v1/RECO
#  /store/relval/CMSSW_10_6_0/SingleMuon/RECO/106X_dataRun2_v10_RelVal_2017F-v1/10000/515ED924-A207-2146-9C0F-6F99A20F2DB0.root
# /SingleMuon/CMSSW_10_6_0-106X_dataRun2_v10_RelVal_2016H-v1/RECO
#  /store/relval/CMSSW_10_6_0/SingleMuon/RECO/106X_dataRun2_v10_RelVal_2016H-v1/10000/BC6B00A1-72AF-FF44-861E-87DD48C2F41C.root
#

print events

debug = False

for iev,event in enumerate(events):
    #if iev >= 100: break 
    #event.getByLabel(muonLabel, muons)
    #event.getByLabel(electronLabel, electrons)
    #event.getByLabel(photonLabel, photons)
    #event.getByLabel(tauLabel, taus)
    #event.getByLabel(jetLabel, jets)
    #event.getByLabel(fatjetLabel, fatjets)
    #event.getByLabel(metLabel, mets)
    event.getByLabel(vertexLabel, vertices)
    #event.getByLabel(vertexLabel, verticesScore)
    event.getByLabel(pfcandLabel,pfcands)
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

    # PF candidates
    for i,j in enumerate(pfcands.product()):  # loop over pf candidates
        #if (abs(j.pdgId())!=211): continue    # charged hadron
        if (j.pdgId()!=130): continue         # only neutral hadrons
        #if j.ecalEnergy() == 0: continue      # ecalEnergy>0
        #if j.pt() < 5: continue               # pt>5 GeV
        #if j.pt() < 20: continue              # pt>20 GeV
        #print iev,i,j,jelementsInBlocks().size()
        #print "iev,i,j.eleinblock.size: ",iev,i,j.elementsInBlocks().size(),

        if ( debug ):
            if (abs(j.pdgId())==211):
                print "PF: iev %3d pfcands %3d: pt %5.1f eta %5.2f pdgId %5d trackpt: %5.2f ecalE %5.2f (raw: %5.2f) hcalE %5.2f (raw: %5.2f) " % (
                    iev, i, j.pt(), j.eta(), j.pdgId(), j.trackRef().pt(), j.ecalEnergy(), j.rawEcalEnergy(), j.hcalEnergy(), j.rawHcalEnergy() )
            else: 
                print "PF: iev %3d pfcands %3d: pt %5.1f eta %5.2f pdgId %5d ecalE %5.2f (raw: %5.2f) hcalE %5.2f (raw: %5.2f) " % (
                    iev, i, j.pt(), j.eta(), j.pdgId(), j.ecalEnergy(), j.rawEcalEnergy(), j.hcalEnergy(), j.rawHcalEnergy() )

        ntrack=0
        nhcal=0
        Tlv_trk = ROOT.TLorentzVector(0.,0.,0.,0.)
        Tlv_trk_tmp = ROOT.TLorentzVector(0.,0.,0.,0.)
        Tlv_hcal = ROOT.TLorentzVector(0.,0.,0.,0.)
        
        for e in range(0, j.elementsInBlocks().size()):  # loop over elements
            for iEle in range(0, j.elementsInBlocks()[e].first.elements().size()):  
                if (j.elementsInBlocks()[e].first.elements()[iEle].index() == j.elementsInBlocks()[e].second):
                    element = j.elementsInBlocks()[e].first.elements()[iEle]
                    if (element.type() == 1):  # track element (1, reco::PFBlockElement::TRACK)
                        track = element.trackRef()
                        if (debug): print " element",e,element.type(),track.pt(),track.eta(),track.phi(),"track"
                        Tlv_trk_tmp.SetPtEtaPhiE(track.pt(),track.eta(),track.phi(),track.p())
                        Tlv_trk += Tlv_trk_tmp
                        ntrack += 1
                    elif (element.type() == 5):  # calo element (5 reco::PFBlockElement::HCAL)
                        cluster = element.clusterRef()
                        #fracs = cluster.recHitFractions()
                        hfracs = cluster.hitsAndFractions()
                        nhcal += 1
                        if (debug): print " element",e,element.type(),cluster.pt(),cluster.eta(),cluster.phi(),"hcal",len(hfracs)
                        H_PFH_HcalClusterSize.Fill(len(hfracs))
                        H_PFH_HcalClusterSizeVsEnergy.Fill(cluster.energy(),len(hfracs))
                        Tlv_hcal.SetPtEtaPhiE(cluster.pt(),cluster.eta(),cluster.phi(),cluster.energy())
                        for k in range(0, hfracs.size()):  # loop over rechits
                            #print "  ", fracs[k], hfracs[k], hfracs[k].first, hfracs[k].second
                            id = hfracs[k].first.rawId()
                            if (debug): print "  ", id, hfracs[k].second
                            #print "  ", id, cluster.recHitFractions()[k].recHitRef().detId()
                            #print HcalDetId(id).ieta(),HcalDetId(id).iphi(),HcalDetId(id).depth()
                            # extracting ieta,iphi from detid: not working
                    elif (element.type() == 4 or element.type() == 11):  # calo element (4,11, reco::PFBlockElement::ECAL,HO)
                        cluster = element.clusterRef()
                        if (debug): print " element",e,element.type(),cluster.pt(),"ecal or ho"
                    else:
                        if (debug): print " element",e,element.type()
        if (ntrack>0):
            H_PFH_withTrack.Fill(j.eta(),1.)
            H_PFH_withoutTrack.Fill(j.eta(),0.)
            if (nhcal>0):
                if (Tlv_hcal.DeltaR(Tlv_trk)>0.4):
                    if (debug): print "KH warning: dR>0.4"
                    #    print " hcal:",Tlv_hcal.Print()
                    #    print " track:",Tlv_trk.Print()
                H_PFH_HcalTrackdR.Fill(Tlv_hcal.DeltaR(Tlv_trk));            
        else:
            H_PFH_withTrack.Fill(j.eta(),0.)
            H_PFH_withoutTrack.Fill(j.eta(),1.)
            
        H_PFH_HCAL.Fill(j.hcalEnergy())
        H_PFH_CorECAL.Fill(j.ecalEnergy())
        H_PFH_RawECAL.Fill(j.rawEcalEnergy())

    # # pfcands (pt>20 GeV, ecal+hcal<>0., pdgid==130)
    # for i,j in enumerate(pfcands.product()):
    #     if j.pt() < 20: continue
    #     if ((j.ecalEnergy()+j.hcalEnergy())==0.): continue
    #     if (j.pdgId()!=130): continue
    #     print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
    #     print "PFH: pfcands %3d: pt %5.1f eta %5.2f pdgId %5d ecalE %5.2f hcalE %5.2f ecalFrac %5.2f" % ( i, j.pt(), j.eta(), j.pdgId(), j.ecalEnergy(), j.hcalEnergy(), j.ecalEnergy()/(j.ecalEnergy()+j.hcalEnergy()) )

    # # pfcands (pt>20 GeV, pdgid==211)
    # for i,j in enumerate(pfcands.product()):
    #     #if j.pt() >2: continue
    #     if j.pt() < 20: continue
    #     if (abs(j.pdgId())!=211): continue
    #     print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
    #     print "CH1: pfcands %3d: pt %5.1f eta %5.2f pdgId %5d ecalE %5.2f hcalE %5.2f " % ( i, j.pt(), j.eta(), j.pdgId(), j.ecalEnergy(), j.hcalEnergy() )

    # pfcands (pdgid==211, ecalEnergy<0.)
    for i,j in enumerate(pfcands.product()):
        if (abs(j.pdgId())!=211): continue
        H_PFH_CorECAL_zoom.Fill(j.ecalEnergy())
        H_PFH_RawECAL_zoom.Fill(j.rawEcalEnergy())
        if (j.ecalEnergy()>=0.): continue 
        print "\nEvent %d: run %6d, lumi %4d, event %12d" % (iev,event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())
        print "CH2: pfcands %3d: pt %5.1f eta %5.2f pdgId %5d ecalE %5.2f (rawEcalE %5.2f) hcalE %5.2f (rawHcalE %5.2f) " % ( i, j.pt(), j.eta(), j.pdgId(), j.ecalEnergy(), j.rawEcalEnergy(), j.hcalEnergy(), j.rawHcalEnergy() )

# Set up canvas : 
w = 1400 
h =  700
can  = ROOT.TCanvas("can", "histograms   ", w, h)

#####
ROOT.gPad.SetLogy()
H_PFH_RawECAL.SetLineColor(4)
H_PFH_RawECAL.SetMinimum(0.5)
H_PFH_RawECAL.Draw()
H_PFH_CorECAL.SetLineColor(2)
H_PFH_CorECAL.Draw("sames")
H_PFH_HCAL.Draw("sames")

legend = ROOT.TLegend(0.6, 0.6, 0.85, 0.85)
legend.AddEntry(H_PFH_HCAL, "HCAL", "l")
legend.AddEntry(H_PFH_RawECAL, "RawECAL", "l")
legend.AddEntry(H_PFH_CorECAL, "CorECAL", "l")
legend.Draw()

can.SaveAs(output+".pdf")
can.SaveAs(output+".png")
can.SaveAs(output+".jpg")
can.SaveAs(output+".root")

#####
ROOT.gPad.SetLogy()
H_PFH_CorECAL_zoom.SetLineColor(4)
H_PFH_CorECAL_zoom.SetMinimum(0.5)
H_PFH_CorECAL_zoom.Draw()

can.SaveAs(output+"_ch_zoom.pdf")
can.SaveAs(output+"_ch_zoom.png")
can.SaveAs(output+"_ch_zoom.jpg")
can.SaveAs(output+"_ch_zoom.root")

#####
ROOT.gPad.SetLogy()
H_PFH_RawECAL_zoom.SetLineColor(4)
H_PFH_RawECAL_zoom.SetMinimum(0.5)
H_PFH_RawECAL_zoom.Draw()

can.SaveAs(output+"_ch_zoom.pdf")
can.SaveAs(output+"_ch_zoom.png")
can.SaveAs(output+"_ch_zoom.jpg")
can.SaveAs(output+"_ch_zoom.root")

#####
ROOT.gPad.SetLogy(0)
H_PFH_withTrack.Draw()
can.SaveAs(output+"_track.pdf")
can.SaveAs(output+"_track.png")
can.SaveAs(output+"_track.jpg")
can.SaveAs(output+"_track.root")

#####
ROOT.gPad.SetLogy(0)
H_PFH_withoutTrack.SetMaximum(1.5)
H_PFH_withoutTrack.SetMinimum(0.0)
H_PFH_withoutTrack.Draw()
can.SaveAs(output+"_wotrack.pdf")
can.SaveAs(output+"_wotrack.png")
can.SaveAs(output+"_wotrack.jpg")
can.SaveAs(output+"_wotrack.root")

#####
ROOT.gPad.SetLogy(0)
H_NPV.Draw()
can.SaveAs(output+"_NPV.pdf")
can.SaveAs(output+"_NPV.png")
can.SaveAs(output+"_NPV.jpg")
can.SaveAs(output+"_NPV.root")

#####
ROOT.gPad.SetLogy(1)
H_PFH_HcalTrackdR.Draw()
can.SaveAs(output+"_HcalTrackdR.pdf")
can.SaveAs(output+"_HcalTrackdR.png")
can.SaveAs(output+"_HcalTrackdR.jpg")
can.SaveAs(output+"_HcalTrackdR.root")

#####
ROOT.gPad.SetLogy(1)
H_PFH_HcalClusterSize.Draw()
can.SaveAs(output+"_HcalClusterSize.pdf")
can.SaveAs(output+"_HcalClusterSize.png")
can.SaveAs(output+"_HcalClusterSize.jpg")
can.SaveAs(output+"_HcalClusterSize.root")

#####
ROOT.gPad.SetLogy(0)
H_PFH_HcalClusterSizeVsEnergy.Draw("COLZ")
can.SaveAs(output+"_HcalClusterSizeVsEnergy.pdf")
can.SaveAs(output+"_HcalClusterSizeVsEnergy.png")
can.SaveAs(output+"_HcalClusterSizeVsEnergy.jpg")
can.SaveAs(output+"_HcalClusterSizeVsEnergy.root")
