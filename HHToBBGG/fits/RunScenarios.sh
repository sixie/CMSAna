#!/bin/bash

clear

for s in {0..12}
do
	root -l -b -q CMSAna/HHToBBGG/fits/FitScenariosTwoD.C+\(\"/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root\",\"/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root\",\"/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root\",\"pho\",5000,$s\) | tee PhotonEnergyResolutionScan_${s}.log
done


for s in {0..12}
do
	root -l -b -q CMSAna/HHToBBGG/fits/FitScenariosTwoD.C+\(\"/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root\",\"/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root\",\"/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root\",\"jet\",5000,$s\) | tee JetEnergyResolutionScan_${s}.log
done



for s in {0..12}
do
	root -l -b -q CMSAna/HHToBBGG/fits/FitScenariosTwoD.C+\(\"/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_diphotonMass.root\",\"/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_dibjetMass.root\",\"/afs/cern.ch/work/d/daan/public/releases/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/data/HHToBBGG_SignalBkgd_AfterCuts_twoDMass.root\",\"lum\",5000,$s\)
done
