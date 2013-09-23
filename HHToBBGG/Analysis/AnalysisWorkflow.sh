####################################
# Make Electron To Photon Fake Rate
####################################
root -l CMSAna/ObjectStudies/Photons/ComputePhotonEfficiency.C+'("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/ElectronToPhotonFakeRate/ntuples/ElectronToPhotonFakeRateNtuple.ttjll-START53_V7A.root",11,-1,"ElectronToPhotonFakeRate")'


####################################
# Plot Efficiencies
####################################
root -l CMSAna/HHToBBGG/Analysis/PlotEfficiencies.C+'("CMSAna/HHToBBGG/data/PhotonEfficiency_PromptPhoton.root","Efficiency_PtEta","Photon p_{T} [GeV/c]","Photon #eta","PromptPhoton", 0.0, 1.0)'
root -l CMSAna/HHToBBGG/Analysis/PlotEfficiencies.C+'("CMSAna/HHToBBGG/data/BTaggingEfficiency_BJetEfficiency.root","Efficiency_PtEta","B Jet p_{T} [GeV/c]","B Jet #eta","BJet", 0.0, 1.0)'


root -l CMSAna/HHToBBGG/Analysis/PlotEfficiencies.C+'("CMSAna/HHToBBGG/data/PhotonEfficiency_GluonJetFakeRate.root","MistagRate_CSVMedium_PtEta","Photon p_{T} [GeV/c]","Photon #eta","GluonJetFakesPhoton", 0.00001, 0.1, true)'
root -l CMSAna/HHToBBGG/Analysis/PlotEfficiencies.C+'("CMSAna/HHToBBGG/data/PhotonEfficiency_QuarkJetFakeRate.root","MistagRate_CSVMedium_PtEta","Photon p_{T} [GeV/c]","Photon #eta","QuarkJetFakesPhoton", 0.0001, 0.1, true)'

root -l CMSAna/HHToBBGG/Analysis/PlotEfficiencies.C+'("CMSAna/HHToBBGG/data/BTaggingEfficiency_LightJetsMistagRate.root","MistagRate_CSVMedium_Pt_Eta","Jet p_{T} [GeV/c]","Jet #eta","LightJetMistag", 0.0, 0.05, false, 30, 150)'
root -l CMSAna/HHToBBGG/Analysis/PlotEfficiencies.C+'("CMSAna/HHToBBGG/data/BTaggingEfficiency_CharmJetsMistagRate.root","MistagRate_CSVMedium_Pt_Eta","Jet p_{T} [GeV/c]","Jet #eta","CharmJetMistag", 0.0, 0.5, false, 30, 150)'

root -l CMSAna/HHToBBGG/Analysis/PlotEfficiencies.C+'("CMSAna/HHToBBGG/data/PhotonEfficiency_ElectronToPhotonFakeRate.root","Efficiency_PtEta","Electron p_{T} [GeV/c]","Electron #eta","ElectronFakesPhoton", 0.001, 0.5, true)'



####################################
# ROOT DRAW String
####################################
##Cut Stage1
 HHToBBGGEvent->Draw("diphoton.M()","3000*weight*(pho1.Pt()>25 && pho2.Pt()>25 && max(pho1.Pt(),pho2.Pt())>40 && abs(pho1.Eta())<2.5 && abs(pho2.Eta())<2.5 && bjet1.Pt()>30 && bjet2.Pt() > 30 && abs(bjet1.Eta())<2.4 && abs(bjet2.Eta())<2.4 && diphoton.M() > 100 && diphoton.M() < 150 && dibjet.M() > 60 && dibjet.M() < 200)");

##Cut Stage2
 HHToBBGGEvent->Draw("diphoton.M()","3000*weight*(pho1.Pt()>25 && pho2.Pt()>25 && max(pho1.Pt(),pho2.Pt())>40 && abs(pho1.Eta())<2.5 && abs(pho2.Eta())<2.5 && bjet1.Pt()>30 && bjet2.Pt() > 30 && abs(bjet1.Eta())<2.4 && abs(bjet2.Eta())<2.4 && diphoton.M() > 120 && diphoton.M() < 130 && dibjet.M() > 105 && dibjet.M() < 145)");

##Cut Stage3
 HHToBBGGEvent->Draw("diphoton.Pt()","3000*weight*(pho1.Pt()>25 && pho2.Pt()>25 && max(pho1.Pt(),pho2.Pt())>40 && abs(pho1.Eta())<2.5 && abs(pho2.Eta())<2.5 && bjet1.Pt()>30 && bjet2.Pt() > 30 && abs(bjet1.Eta())<2.4 && abs(bjet2.Eta())<2.4 && diphoton.M() > 120 && diphoton.M() < 130 && dibjet.M() > 105 && dibjet.M() < 145 && DRgg < 2 && minDRgb>1 && ncentraljets < 4)");


 HHToBBGGEvent->Draw("diphoton.M()","3000*weight*(pho1.Pt()>25 && pho2.Pt()>25 && abs(pho1.Eta())<2.5 && abs(pho2.Eta())<2.5 && bjet1.Pt()>30 && bjet2.Pt() > 30 && abs(bjet1.Eta())<2.4 && abs(bjet2.Eta())<2.4 && DRgg < 2 && minDRgb>1)");


####################################
# Cut-Based Analysis
####################################

###bbgg
root -l CMSAna/HHToBBGG/Analysis/CutBasedAnalysis.C+'("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/normalizedNtuples/HHToBBGGNtuple.DiPhotonBB_M60To200_14TeV-BBGG.normalized.root")'

###ggjj
root -l CMSAna/HHToBBGG/Analysis/CutBasedAnalysis.C+'("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/normalizedNtuples/HHToBBGGNtuple.DiPhotonJJ_M60To200_14TeV.normalized.root")'
root -l CMSAna/HHToBBGG/Analysis/CutBasedAnalysis.C+'("/afs/cern.ch/work/v/vlambert/public/HHbbggBackground/new_normalizedNtuples/HHToBBGGNtuple.DiPhotonJJ_M60To200_14TeV.normalized.root")'

###bbjj
root -l CMSAna/HHToBBGG/Analysis/CutBasedAnalysis.C+'("/afs/cern.ch/work/v/vlambert/public/HHbbggBackground/new_normalizedNtuples/BACONNtuple_GenOnly_BBJJ_M110To140_14TeV.normalized.root")'

###ccjj
root -l CMSAna/HHToBBGG/Analysis/CutBasedAnalysis.C+'("/afs/cern.ch/work/v/vlambert/public/HHbbggBackground/new_normalizedNtuples/BACONNtuple_GenOnly_CCJJ_M110To140_14TeV.normalized.root")'

###jjjj
root -l CMSAna/HHToBBGG/Analysis/CutBasedAnalysis.C+'("/afs/cern.ch/work/v/vlambert/public/HHbbggBackground/new_normalizedNtuples/BACONNtuple_GenOnly_JJJJ_14TeV.root")'

