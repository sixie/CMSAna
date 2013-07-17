#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================


################
# Photon Efficiency
################
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/diphjets-START53_V7A.txt /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/diphjets2-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonEfficiencyTrees_diphjets-START53_V7A_${i}.out -J Phase2Upgrade_PhotonEfficiencyTrees_diphjets-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonEfficiencyNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonEfficiencyNtuple.diphjets-START53_V7A.${i}.root\",10\) PhotonEfficiencyNtuple.diphjets-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonEfficiency/ntuples/jobs/
  @ i = $i + 1
end






################
# Summer12 QCD
################
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt30To50-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonFakeRateTrees_qcd_pt30To50-START53_V7A_${i}.out -J Phase2Upgrade_PhotonFakeRateTrees_qcd_pt30To50-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonFakeRateNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonFakeRateNtuple.qcd_pt30To50-START53_V7A.${i}.root\",0\) PhotonFakeRateNtuple.qcd_pt30To50-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonFakeRate/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt50To80-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonFakeRateTrees_qcd_pt50To80-START53_V7A_${i}.out -J Phase2Upgrade_PhotonFakeRateTrees_qcd_pt50To80-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonFakeRateNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonFakeRateNtuple.qcd_pt50To80-START53_V7A.${i}.root\",0\) PhotonFakeRateNtuple.qcd_pt50To80-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonFakeRate/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt80To120-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonFakeRateTrees_qcd_pt80To120-START53_V7A_${i}.out -J Phase2Upgrade_PhotonFakeRateTrees_qcd_pt80To120-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonFakeRateNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonFakeRateNtuple.qcd_pt80To120-START53_V7A.${i}.root\",0\) PhotonFakeRateNtuple.qcd_pt80To120-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonFakeRate/ntuples/jobs/
  @ i = $i + 1
end


################
# Phase1 TTHbb
################
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/ttHbb-125-Age0_STAR17_61_V1A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonFakeRateTrees_ttHbb-125_STAR17_61_V1A.root_${i}.out -J Phase2Upgrade_PhotonFakeRateTrees_ttHbb-125_STAR17_61_V1A.root_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonFakeRateNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonFakeRateNtuple.Phase1Age0.ttHbb-125_STAR17_61_V1A.root.${i}.root\",0\) PhotonFakeRateNtuple.Phase1Age0.ttHbb-125_STAR17_61_V1A.root.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonFakeRate/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/ttHbb-125-Age0DES_DES17_61_V5.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonFakeRateTrees_ttHbb-125-DES_DES17_61_V5_${i}.out -J Phase2Upgrade_PhotonFakeRateTrees_ttHbb-125-DES_DES17_61_V5_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonFakeRateNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonFakeRateNtuple.Phase1Age0.ttHbb-125-DES_DES17_61_V5.${i}.root\",0\) PhotonFakeRateNtuple.Phase1Age0.ttHbb-125-DES_DES17_61_V5.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonFakeRate/ntuples/jobs/
  @ i = $i + 1
end








################
# Electron -> Photon Fake Rate
################

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/ttjll-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonEfficiencyTrees_ttjll-START53_V7A_${i}.out -J Phase2Upgrade_PhotonEfficiencyTrees_ttjll-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonEfficiencyNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"ElectronToPhotonFakeRateNtuple.ttjll-START53_V7A.${i}.root\",20\) ElectronToPhotonFakeRateNtuple.ttjll-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/ElectronToPhotonFakeRate/ntuples/jobs/
  @ i = $i + 1
end

