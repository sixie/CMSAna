\#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================




####################################################
#
# High Pileup Samples
#
####################################################


#############
# WZHgg
#############

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/WZHgg-125-Age0_STAR17_61_V1A.txt`)
  echo $file 
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_WZHgg-125-Age0_STAR17_61_V1A_${i}.out -J Phase2Upgrade_PhotonTrees_WZHgg-125-Age0_STAR17_61_V1A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeJetNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"JetNtuple.WZHgg-125-Age0_STAR17_61_V1A.${i}.root\",-1,-1\) JetNtuple.WZHgg-125-Age0_STAR17_61_V1A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/JetStudies/ntuples/jobs/
  sleep 0.5
  @ i = $i + 1
end

#############
# DY
#############
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/DYToEEAge0START_STAR17_61_V1A.txt /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/DYToMMAge0START_STAR17_61_V1A.txt`)
  echo $file 
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_DYAge0_${i}.out -J Phase2Upgrade_PhotonTrees_DYAge0_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeJetNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"JetNtuple.DYAge0.${i}.root\",-1,-1\) JetNtuple.DYAge0.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/JetStudies/ntuples/jobs/
  sleep 0.5
  @ i = $i + 1
end

#############
# ttHbb
#############
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/ttHbb-125-Age0_STAR17_61_V1A.txt`)
  echo $file 
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/JetTrees_ttHbb-125-Age0_STAR17_61_V1A_${i}.out -J Phase2Upgrade_JetTrees_ttHbb-125-Age0_STAR17_61_V1A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeJetNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"JetNtuple.ttHbb-125-Age0_STAR17_61_V1A.${i}.root\",-1,-1\) JetNtuple.ttHbb-125-Age0_STAR17_61_V1A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/JetStudies/ntuples/jobs/
  sleep 0.5
  @ i = $i + 1
end

