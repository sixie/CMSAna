#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================



################
# Higgs->bb
################
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/HHtoBBGG-14tev-START53_V7A.txt /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/HHtoBBGG-8tev-START53_V7A.txt `)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/BJetEfficiencyTrees_HHtoBBGG-14tev-START53_V7A_${i}.out -J Phase2Upgrade_BJetEfficiencyTrees_HHtoBBGG-14tev-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeBJetEfficiencyNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"BJetEfficiencyNtuple.HHtoBBGG-14tev-START53_V7A.${i}.root\",0\) BJetEfficiencyNtuple.HHtoBBGG-14tev-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetEfficiencyRate/ntuples/jobs/
  @ i = $i + 1
end

################
# ttbar
################
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/ttjll-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/BJetEfficiencyTrees_ttjll-START53_V7A_${i}.out -J Phase2Upgrade_BJetEfficiencyTrees_ttjll-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeBJetEfficiencyNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"BJetEfficiencyNtuple.ttjll-START53_V7A.${i}.root\",0\) BJetEfficiencyNtuple.ttjll-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetEfficiencyRate/ntuples/jobs/
  @ i = $i + 1
end


################
# Summer12 QCD
################
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt30To50-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/BJetEfficiencyTrees_qcd_pt30To50-START53_V7A_${i}.out -J Phase2Upgrade_BJetEfficiencyTrees_qcd_pt30To50-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeBJetEfficiencyNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"BJetEfficiencyNtuple.qcd_pt30To50-START53_V7A.${i}.root\",0\) BJetEfficiencyNtuple.qcd_pt30To50-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetEfficiencyRate/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt50To80-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/BJetEfficiencyTrees_qcd_pt50To80-START53_V7A_${i}.out -J Phase2Upgrade_BJetEfficiencyTrees_qcd_pt50To80-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeBJetEfficiencyNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"BJetEfficiencyNtuple.qcd_pt50To80-START53_V7A.${i}.root\",0\) BJetEfficiencyNtuple.qcd_pt50To80-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetEfficiencyRate/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt80To120-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/BJetEfficiencyTrees_qcd_pt80To120-START53_V7A_${i}.out -J Phase2Upgrade_BJetEfficiencyTrees_qcd_pt80To120-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeBJetEfficiencyNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"BJetEfficiencyNtuple.qcd_pt80To120-START53_V7A.${i}.root\",0\) BJetEfficiencyNtuple.qcd_pt80To120-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetEfficiencyRate/ntuples/jobs/
  @ i = $i + 1
end


################
# Photon Jets
################
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/diphjets-START53_V7A.txt /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/diphjets2-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/BJetEfficiencyTrees_diphjets-START53_V7A_${i}.out -J Phase2Upgrade_BJetEfficiencyTrees_diphjets-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeBJetEfficiencyNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"BJetEfficiencyNtuple.diphjets-START53_V7A.${i}.root\",0\) BJetEfficiencyNtuple.diphjets-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetEfficiencyRate/ntuples/jobs/
  @ i = $i + 1
end


set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/diphjetsSherpa-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/BJetEfficiencyTrees_diphjetsSherpa-START53_V7A_${i}.out -J Phase2Upgrade_BJetEfficiencyTrees_diphjetsSherpa-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/BTagging/ MakeBJetEfficiencyNtuple.C +\(\"root://eoscms//eos/cms/${file}\",\"BJetEfficiencyNtuple.diphjetsSherpa-START53_V7A.${i}.root\",0\) BJetEfficiencyNtuple.diphjetsSherpa-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/BJetEfficiencyRate/ntuples/jobs/
  @ i = $i + 1
end
