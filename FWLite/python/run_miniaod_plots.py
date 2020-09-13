import math
import copy

def deltaR2( e1, p1, e2=None, p2=None):
    """Take either 4 arguments (eta,phi, eta,phi) or two objects that have 'eta', 'phi' methods)"""
    if (e2 == None and p2 == None):
        return deltaR2(e1.eta(),e1.phi(), p1.eta(), p1.phi())
    de = e1 - e2
    dp = deltaPhi(p1, p2)
    return de*de + dp*dp

def deltaPhi( p1, p2):
    '''Computes delta phi, handling periodic limit conditions.'''
    res = p1 - p2
    while res > math.pi:
        res -= 2*math.pi
    while res < -math.pi:
        res += 2*math.pi
    return res

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

# Create histograms, etc.
ROOT.gROOT.SetStyle('Plain') # white background
H_ElectronEta = ROOT.TH1F ("ElectronEta","ElectronEta",60,-3.,+3.)
H_ElectronEta_Fake = ROOT.TH1F ("ElectronEta_Fake","ElectronEta_Fake",60,-3.,+3.)

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

genpars, genParLabel = Handle("std::vector<reco::GenParticle>"), "prunedGenParticles"
pgenpars, pgenParLabel = Handle("std::vector<pat::PackedGenParticle>"), "packedGenParticles"
muons, muonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"
electrons, electronLabel = Handle("std::vector<pat::Electron>"), "slimmedElectrons"
photons, photonLabel = Handle("std::vector<pat::Photon>"), "slimmedPhotons"
taus, tauLabel = Handle("std::vector<pat::Tau>"), "slimmedTaus"
jets, jetLabel = Handle("std::vector<pat::Jet>"), "slimmedJets"
fatjets, fatjetLabel = Handle("std::vector<pat::Jet>"), "slimmedJetsAK8"
mets, metLabel = Handle("std::vector<pat::MET>"), "slimmedMETs"
vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlineSlimmedPrimaryVertices"
pfcands, pfcandLabel = Handle("std::vector<pat::PackedCandidate>"), "packedPFCandidates"

verticesScore = Handle("edm::ValueMap<float>")

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
#events = Events('file:step4_inMINIAODSIM.root')
#events = Events('file:/cms/data/store/user/hatake/RelValQCD_FlatPt_15_3000HS_14/pfvalidation/200906_033704/0000/step3_inMINIAODSIM_1.root')
#events = Events('file:/cms/data/store/user/hatake/RelValQCD_FlatPt_15_3000HS_14/pfvalidation/200905_204656/0000/step3_inMINIAODSIM_1.root')
# test
#
listFiles=[
     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012317/0000/step3_inMINIAODSIM_1.root',
     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012317/0000/step3_inMINIAODSIM_2.root',
     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012317/0000/step3_inMINIAODSIM_3.root',
     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012317/0000/step3_inMINIAODSIM_5.root',
     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012317/0000/step3_inMINIAODSIM_6.root',
     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012317/0000/step3_inMINIAODSIM_7.root',
     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012317/0000/step3_inMINIAODSIM_8.root',
     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012317/0000/step3_inMINIAODSIM_9.root'
]
# ref
# listFiles=[
#     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012521/0000/step3_inMINIAODSIM_1.root',
#     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012521/0000/step3_inMINIAODSIM_2.root',
#     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012521/0000/step3_inMINIAODSIM_3.root',
#     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012521/0000/step3_inMINIAODSIM_5.root',
#     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012521/0000/step3_inMINIAODSIM_6.root',
#     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012521/0000/step3_inMINIAODSIM_7.root',
#     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012521/0000/step3_inMINIAODSIM_8.root',
#     'file:/cms/data/store/user/hatake/RelValZEE_14/pfvalidation/200908_012521/0000/step3_inMINIAODSIM_9.root'
# ]
events = Events(listFiles)
#events = Events('step3_inMINIAODSIM_rerereco_all.root')

for iev,event in enumerate(events):
    #if iev >= 100: break 
    event.getByLabel(genParLabel, genpars)
    event.getByLabel(pgenParLabel, pgenpars)
    event.getByLabel(muonLabel, muons)
    event.getByLabel(electronLabel, electrons)
    event.getByLabel(photonLabel, photons)
    event.getByLabel(tauLabel, taus)
    event.getByLabel(jetLabel, jets)
    event.getByLabel(fatjetLabel, fatjets)
    event.getByLabel(metLabel, mets)
    event.getByLabel(vertexLabel, vertices)
    event.getByLabel(vertexLabel, verticesScore)
    event.getByLabel(pfcandLabel,pfcands)
    print "\nEvent: run %6d, lumi %4d, event %12d" % (event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event())

    # Vertices 
    if len(vertices.product()) == 0 or vertices.product()[0].ndof() < 4:
        print "Event has no good primary vertex."
        continue
    else:
        PV = vertices.product()[0]
        print "PV at x,y,z = run %6d, event %10d, %+15.13f, %+15.13f, %+16.13f, ndof: %.1f, score: (pt2 of clustered objects) %.11f" % (event.eventAuxiliary().run(), event.eventAuxiliary().event(), PV.x(), PV.y(), PV.z(), PV.ndof(),verticesScore.product().get(0))

    # GenParticles
    for i,genp in enumerate(genpars.product()): 
        if genp.pt() < 5 : continue
        print "genpar: run %6d, event %10d, pt %4.1f, eta %5.2f, phi %5.2f, pdgId %d." % (
            event.eventAuxiliary().run(), event.eventAuxiliary().event(), genp.pt(), genp.eta(), genp.phi(), genp.pdgId())

    # PackedGenParticles
    for i,pgenp in enumerate(pgenpars.product()): 
        if pgenp.pt() < 5 : continue
        print "pgenpar: run %6d, event %10d, pt %4.1f, eta %5.2f, phi %5.2f, pdgId %d." % (
            event.eventAuxiliary().run(), event.eventAuxiliary().event(), pgenp.pt(), pgenp.eta(), pgenp.phi(), pgenp.pdgId())

    # Muons
    for i,mu in enumerate(muons.product()): 
        if mu.pt() < 5 or not mu.isLooseMuon(): continue
        print "muon: run %6d, event %10d, pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d." % (
            event.eventAuxiliary().run(), event.eventAuxiliary().event(), mu.pt(), mu.muonBestTrack().dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV))

    # Electrons
    for i,el in enumerate(electrons.product()):
        if el.pt() < 5: continue
        print "elec: run %6d, event %10d, pt %6.2f, supercluster eta %+5.3f, phi %+5.3f, energy %5.2f (raw %5.2f), gsf pt %+5.2f sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), lost hits %d, pass conv veto %d" % (
                    event.eventAuxiliary().run(), event.eventAuxiliary().event(), el.pt(), el.superCluster().eta(), el.superCluster().phi(), el.superCluster().correctedEnergy(), el.superCluster().energy(), el.gsfTrack().pt(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.gsfTrack().hitPattern().numberOfLostHits(ROOT.reco.HitPattern.MISSING_INNER_HITS), el.passConversionVeto())
        match = False
        for i,genp in enumerate(genpars.product()):
            if abs(genp.pdgId()) == 11:
                if deltaR2(genp,el)<0.01:
                    genmatch = genp
                    match = True
        if match :
            H_ElectronEta.Fill(el.eta())
        else :
            H_ElectronEta_Fake.Fill(el.eta())
        
    # Photon
    for i,pho in enumerate(photons.product()):
        if pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3: continue
        print "phot: run %6d, event %10d, pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)" % (
                    event.eventAuxiliary().run(), event.eventAuxiliary().event(), pho.pt(), pho.superCluster().eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta())

    # Tau
    for i,tau in enumerate(taus.product()):
        if tau.pt() < 20: continue
        print "tau: run %6d, event %10d, pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d " % (
                    event.eventAuxiliary().run(), event.eventAuxiliary().event(), tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand().pt(), tau.leadCand().pdgId()) 

    # Jets (standard AK4)
    for i,j in enumerate(jets.product()):
        if j.pt() < 20: continue
        print "jet: run %6d, event %10d, pt %5.1f (raw pt %5.1f, matched-calojet pt %5.1f), eta %+4.2f, btag run1(CSV) ) %.3f, run2(pfCSVIVFV2) %.3f, pileup mva disc %+.2f" % (
            event.eventAuxiliary().run(), event.eventAuxiliary().event(), j.pt(), j.pt()*j.jecFactor('Uncorrected'), j.userFloat("caloJetMap:pt"), j.eta(), max(0,j.bDiscriminator("combinedSecondaryVertexBJetTags")), max(0,j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")), j.userFloat("pileupJetId:fullDiscriminant"))
        # if i == 0: # for the first jet, let's print the leading constituents
        #     constituents = [ j.daughter(i2) for i2 in xrange(j.numberOfDaughters()) ]
        #     constituents.sort(key = lambda c:c.pt(), reverse=True)
        #     for i2, cand in enumerate(constituents):
        #         if i2 > 4: 
        #                 print "         ....."
        #                 break
        #         print "         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d" % (i2,cand.pt(),cand.dz(PV.position()),cand.pdgId()) 

    # pfcands
    for i,j in enumerate(pfcands.product()):
        if j.pt() < 0: continue
        print "pfcands: run %6d, event %10d, pt %5.1f eta %5.2f pdgId %5d %5.3f %5.3f " % ( event.eventAuxiliary().run(), event.eventAuxiliary().event(), j.pt(), j.eta(), j.pdgId(), j.rawCaloFraction(), j.rawHcalFraction())

    # Fat AK8 Jets
    # for i,j in enumerate(fatjets.product()):
    #     print "jetAK8 %3d: pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.1f softdrop, %5.1f pruned, %5.1f trimmed, %5.1f filtered. CMS TopTagger %.1f" % (
    #         i, j.pt(), j.pt()*j.jecFactor('Uncorrected'), j.eta(), j.mass(), j.userFloat('ak8PFJetsCHSSoftDropMass'), j.userFloat('ak8PFJetsCHSPrunedMass'), j.userFloat('ak8PFJetsCHSTrimmedMass'), j.userFloat('ak8PFJetsCHSFilteredMass'), j.userFloat("cmsTopTagPFJetsCHSMassAK8"))

    #     # To get the constituents of the AK8 jets, you have to loop over all of the
    #     # daughters recursively. To save space, the first two constituents are actually
    #     # the Soft Drop SUBJETS, which will then point to their daughters.
    #     # The remaining constituents are those constituents removed by soft drop but
    #     # still in the AK8 jet.
    #     constituents = []
    #     for ida in xrange( j.numberOfDaughters() ) :
    #         cand = j.daughter(ida)
    #         if cand.numberOfDaughters() == 0 :
    #             constituents.append( cand )
    #         else :
    #             for jda in xrange( cand.numberOfDaughters() ) :
    #                 cand2 = cand.daughter(jda)
    #                 constituents.append( cand2 )
    #     constituents.sort(key = lambda c:c.pt(), reverse=True)
    #     for i2, cand in enumerate(constituents):
    #         if i2 > 4: 
    #                     print "         ....."
    #                     break
    #         print "         constituent %3d: pt %6.2f, pdgId %+3d, #dau %+3d" % (i2,cand.pt(),cand.pdgId(), cand.numberOfDaughters()) 
                    
    #     wSubjets = j.subjets('SoftDrop')        
    #     for iw,wsub in enumerate( wSubjets ) :
    #         print "   w subjet %3d: pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f " % (
    #             iw, wsub.pt(), wsub.pt()*wsub.jecFactor('Uncorrected'), wsub.eta(), wsub.mass()
    #             )
    #     tSubjets = j.subjets('CMSTopTag')
    #     for it,tsub in enumerate( tSubjets ) :
    #         print "   t subjet %3d: pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f " % (
    #             it, tsub.pt(), tsub.pt()*tsub.jecFactor('Uncorrected'), tsub.eta(), tsub.mass()
    #             )

f = ROOT.TFile.Open("myfile.root","RECREATE")
H_ElectronEta.Write()
H_ElectronEta_Fake.Write()
#ROOT.TFile.Close(f)
f.Write()
f.Close()
