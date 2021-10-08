# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideCMSDataAnalysisSchoolPisa2019ParticleFlowShortExercise
# Load PyROOT
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

# Load PyFWLite
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()
from collections import OrderedDict
from DataFormats.FWLite import Handle, Events

# other utilities
import sys
from math import sqrt
from PhysicsTools.HeppyCore.utils.deltar import deltaR, bestMatch
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes

def isAncestor(a,p) :
    if a == p :
        return True
    for i in xrange(0,p.numberOfMothers()) :
        if isAncestor(a,p.mother(i)) :
            return True
    return False

# define handles to read things from the event
handles = OrderedDict()
handles['candsEle'] = [Handle("std::vector<reco::GsfElectron>"), ("gedGsfElectrons", "", "RECO")]
handles['candsPho'] = [Handle("std::vector<reco::Photon>"), ("gedPhotons", "", "RECO")]
handles['cands'  ] = [Handle("std::vector<reco::PFCandidate>"), ("particleFlow"    , "", "RECO")]
handles['gens'   ] = [Handle("std::vector<reco::GenParticle>"), ("genParticles"    , "", "RECO")]
handles['pftaus' ] = [Handle('std::vector<reco::PFTau>')      , ("hpsPFTauProducer", "", "RECO")]

fname = sys.argv[1] # take the name of the file to imput from the command line
events = Events(fname)
print("Reading %s " % fname)
for event in events:
    print("\n", ("--" * 75))
    print("Processing run %d, lumi %d, event %d" % (event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(), event.eventAuxiliary().event()))

    for k, v in handles.iteritems():
        event.getByLabel(v[1], v[0])
        setattr(events, k, list(v[0].product()))
        
    # gen taus
    taus = [gp for gp in event.gens if gp.isLastCopy() and abs(gp.pdgId())==15]

    print("Generated %d taus:" %len(taus))
    
    for itau in taus:
        itau.final_state_particles = []
        itau.visp4 = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')()
        
        for igen, gen in enumerate(event.gens):
            if gen.status()!=1: continue # only final state particles
            if abs(gen.pdgId())==15: continue # no taus (redundant)
            if gen.isDirectPromptTauDecayProductFinalState() and abs(gen.pdgId()) in [12, 14, 16]: continue # neutrinos we know
            if isAncestor(itau, gen):
                itau.final_state_particles.append(gen)
                itau.visp4 += gen.p4()
        
        # print only 1 prong + 1/2 pi zero taus
        itau.gendm = tauDecayModes.genDecayMode(itau.final_state_particles)
        if itau.gendm not in ['kOneProng1PiZero','kOneProng2PiZero']:
            continue

        # if there's a reco tau and it passed DecayModeFinding - with the correct decay mode - we're happy anyways
        pftau = None
        dr2best = 1000.
        (pftau, dr2best) = bestMatch(itau.visp4, event.pftaus)
        if pftau and dr2best<0.3**2 and pftau.decayMode() in [1,2]:
            continue
        if pftau and dr2best>0.3**2:
            pftau = None
                        
        # print the tau information only if fishy, that is there are no PF pions
        reco_pions = []
        for gen in itau.final_state_particles:            
            # find the best PF candidate that matches
            (cand, dr2best) = bestMatch(gen, event.cands)
            if cand and dr2best < 0.4**2:
                if abs(cand.pdgId())==211:
                    reco_pions.append(cand)
                    
        if len(reco_pions)==0:
            print('\n tau        pt  %7.2f  eta %+5.2f  phi %+5.2f  pdgId %+5d\tgen decay mode: %s' %(itau.pt(), itau.eta(), itau.phi(), itau.pdgId(), itau.gendm))
            print(  ' tau vis    pt  %7.2f  eta %+5.2f  phi %+5.2f' %(itau.visp4.pt(), itau.visp4.eta(), itau.visp4.phi()))
            if pftau:
                print(' reco pftau pt  %7.2f  eta %+5.2f  phi %+5.2f \treco decay mode: %d' %(pftau.pt(), pftau.eta(), pftau.phi(), pftau.decayMode()))
            
            for gen in itau.final_state_particles:            
                fakeElectron = False
                fakePhoton = False
                # find the best PF candidate that matches
                (cand, dr2best) = bestMatch(gen, event.cands)
                print("\t gen   pt %7.2f  eta %+5.2f  phi %+5.2f  energy %7.2f   pdgId % +5d  charge %+1d -- mother particle pt %7.2f  eta %+5.2f  phi %+5.2f  pdgId %+5d" % (gen.pt(), gen.eta(), gen.phi(), gen.energy(), gen.pdgId(), gen.charge(), gen.mother(0).pt(), gen.mother(0).eta(), gen.mother(0).phi(), gen.mother(0).pdgId())),
                if cand and dr2best < 0.3**2:
                    print("  -->  best match PFC (deltaR %.3f): cand #%4d   pt %7.2f  eta %+5.2f  phi %+5.2f  energy %7.2f  charge %+1d   pdgId % +5d " % (sqrt(dr2best), event.cands.index(cand), cand.pt(), cand.eta(), cand.phi(), cand.energy(), cand.charge(), cand.pdgId()))
                    if abs(cand.pdgId()) == 11: fakeElectron = True
                    if abs(cand.pdgId()) == 22: fakePhoton = True
                else:
                    print("  --X  no match")
                    
                # find the best PF candidate that matches -- electron --
                if fakeElectron:
                    (cand, dr2best) = bestMatch(gen, event.candsEle)
                    print("\t gen   pt %7.2f  eta %+5.2f  phi %+5.2f  energy %7.2f   pdgId % +5d  charge %+1d -- mother particle pt %7.2f  eta %+5.2f  phi %+5.2f  pdgId %+5d" % (gen.pt(), gen.eta(), gen.phi(), gen.energy(), gen.pdgId(), gen.charge(), gen.mother(0).pt(), gen.mother(0).eta(), gen.mother(0).phi(), gen.mother(0).pdgId())),
                    if cand and dr2best < 0.3**2:
                        print("  -->  best match Ele (deltaR %.3f): cand #%4d   pt %7.2f  eta %+5.2f  phi %+5.2f  energy %7.2f  charge %+1d  dr03 %6.1f  H/E %6.3f  mva_e_pi %6.3f  mva_iso %6.3f  sigma_ietaieta %6.3f  1/e-1/p %6.2f  deta %6.3f  dphi %6.3f" \
                            % (sqrt(dr2best), event.candsEle.index(cand), cand.pt(), cand.eta(), cand.phi(), cand.energy(), cand.charge(), cand.dr03TkSumPt()+cand.dr03EcalRecHitSumEt()+cand.dr03HcalTowerSumEt(), cand.hadronicOverEm(), \
                               cand.mva_e_pi(), \
                               cand.mva_Isolated(), \
                               cand.full5x5_sigmaIetaIeta(), \
                               (1.0 - cand.eSuperClusterOverP()) / cand.ecalEnergy(), \
                               abs(cand.deltaEtaSeedClusterTrackAtVtx()), \
                               abs(cand.deltaPhiSuperClusterTrackAtVtx())))
                    else:
                        print("  --X  no match")
                    
                # find the best PF candidate that matches -- photon --
                if fakePhoton:
                    (cand, dr2best) = bestMatch(gen, event.candsPho)
                    print("\t gen   pt %7.2f  eta %+5.2f  phi %+5.2f  energy %7.2f   pdgId % +5d  charge %+1d -- mother particle pt %7.2f  eta %+5.2f  phi %+5.2f  pdgId %+5d" % (gen.pt(), gen.eta(), gen.phi(), gen.energy(), gen.pdgId(), gen.charge(), gen.mother(0).pt(), gen.mother(0).eta(), gen.mother(0).phi(), gen.mother(0).pdgId())),
                    if cand and dr2best < 0.3**2:
                        print("  -->  best match Pho (deltaR %.3f): cand #%4d   pt %7.2f  eta %+5.2f  phi %+5.2f  energy %7.2f  charge %+1d  dr03 %6.1f  H/E %6.3f  sigma_ietaieta %6.3f  " \
                            % (sqrt(dr2best), event.candsPho.index(cand), cand.pt(), cand.eta(), cand.phi(), cand.energy(), cand.charge(), cand.trkSumPtHollowConeDR03()+cand.ecalRecHitSumEtConeDR03()+cand.hcalTowerSumEtConeDR03(), cand.hadTowOverEm(), \
                               cand.sigmaIetaIeta()))
                    else:
                        print("  --X  no match")
                
#                 if abs(cand.pdgId())==11: 
#                     print 'electron: cand #%4d   pt %7.2f  eta %+5.2f  phi %+5.2f  energy %7.2f   pdgId % +5d  charge %+1d  E/P %3.3f   ecal %5.2f   hcal %5.2f   mva_e_pi %3.2f'  %(event.cands.index(cand), cand.pt(), cand.eta(), cand.phi(), cand.energy(), cand.pdgId(), cand.charge(), cand.energy()/cand.p(), cand.ecalEnergy(), cand.hcalEnergy(), cand.mva_e_pi())
                    

