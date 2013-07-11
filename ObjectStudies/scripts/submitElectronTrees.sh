#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================




####################################################
#
# HH -> BB Gamma Gamma Selection
#
####################################################

#for SL5 nodes only
bsub -R "type=SLC5_64"

#############
# TTBAR Summer12
#############
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/ttjll-START53_V7A.txt `)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/ElectronTrees_ttjll-START53_V7A_${i}.out -J Phase2Upgrade_ElectronTrees_ttjll-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/ MakeElectronNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"ElectronNtuple.ttjll-START53_V7A.${i}.root\"\) ElectronNtuple.ttjll-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/ElectronStudies/ntuples/jobs/
  @ i = $i + 1
end


#############
# Zee Phase1 Age0
#############

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/DYToEEAge0.txt `)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/ElectronTrees_DYToEEAge0_${i}.out -J Phase2Upgrade_ElectronTrees_DYToEEAge0_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/ MakeElectronNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"ElectronNtuple.DYToEEAge0.${i}.root\"\) ElectronNtuple.DYToEEAge0.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/ElectronStudies/ntuples/jobs/
  @ i = $i + 1
end


#############
# Zee Phase1 Age3H
#############

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/DYToEEAge3H.txt `)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/ElectronTrees_DYToEEAge3H_${i}.out -J Phase2Upgrade_ElectronTrees_DYToEEAge3H_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/ MakeElectronNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"ElectronNtuple.DYToEEAge3H.${i}.root\"\) ElectronNtuple.DYToEEAge3H.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/ElectronStudies/ntuples/jobs/
  @ i = $i + 1
end


