####################################
# Make Electron To Photon Fake Rate
####################################
root -l CMSAna/ObjectStudies/Photons/ComputePhotonEfficiency.C+'("/afs/cern.ch/work/s/sixie/public/Phase2Upgrade/ElectronToPhotonFakeRate/ntuples/ElectronToPhotonFakeRateNtuple.ttjll-START53_V7A.root",11,-1,"ElectronToPhotonFakeRate")'



####################################
# ROOT DRAW String
####################################
 HHToBBGGEvent->Draw("diphoton.Pt()","3000*weight*(pho1.Pt()>25 && pho2.Pt()>25 && abs(pho1.Eta())<2.5 && abs(pho2.Eta())<2.5 && bjet1.Pt()>30 && bjet2.Pt() > 30 && abs(bjet1.Eta())<2.4 && abs(bjet2.Eta())<2.4 && diphoton.M() > 120 && diphoton.M() < 130 && dibjet.M() > 105 && dibjet.M() < 145 && DRgg < 2 && minDRgb>1)");


 HHToBBGGEvent->Draw("diphoton.M()","3000*weight*(pho1.Pt()>25 && pho2.Pt()>25 && abs(pho1.Eta())<2.5 && abs(pho2.Eta())<2.5 && bjet1.Pt()>30 && bjet2.Pt() > 30 && abs(bjet1.Eta())<2.4 && abs(bjet2.Eta())<2.4 && DRgg < 2 && minDRgb>1)");


