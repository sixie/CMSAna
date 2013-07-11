#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================




################
# Summer12 QCD
################
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt30To50-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/BJetMistagTrees_qcd_pt30To50-START53_V7A_${i}.out -J Phase2Upgrade_BJetMistagTrees_qcd_pt30To50-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeBJetMistagNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"BJetMistagNtuple.qcd_pt30To50-START53_V7A.${i}.root\",0\) BJetMistagNtuple.qcd_pt30To50-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetMistagRate/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt50To80-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/BJetMistagTrees_qcd_pt50To80-START53_V7A_${i}.out -J Phase2Upgrade_BJetMistagTrees_qcd_pt50To80-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeBJetMistagNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"BJetMistagNtuple.qcd_pt50To80-START53_V7A.${i}.root\",0\) BJetMistagNtuple.qcd_pt50To80-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetMistagRate/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt80To120-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/BJetMistagTrees_qcd_pt80To120-START53_V7A_${i}.out -J Phase2Upgrade_BJetMistagTrees_qcd_pt80To120-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeBJetMistagNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"BJetMistagNtuple.qcd_pt80To120-START53_V7A.${i}.root\",0\) BJetMistagNtuple.qcd_pt80To120-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetMistagRate/ntuples/jobs/
  @ i = $i + 1
end


