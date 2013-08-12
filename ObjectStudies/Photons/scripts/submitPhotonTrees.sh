#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================




####################################################
#
# Real Photons
#
####################################################


#############
# HH->bbgg Summer12
#############
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/HHtoBBGG-14tev-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_HHtoBBGG-14tev-START53_V7A_${i}.out -J Phase2Upgrade_PhotonTrees_HHtoBBGG-14tev-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtuple.HHtoBBGG-14tev-START53_V7A.${i}.root\",true\) PhotonNtuple.HHtoBBGG-14tev-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end

#############
# DiPhoton Born Summer12
#############
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/diphjets-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_diphjets-START53_V7A_${i}.out -J Phase2Upgrade_PhotonTrees_diphjets-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtuple.diphjets-START53_V7A.${i}.root\",true\) PhotonNtuple.diphjets-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end


#############
# DiPhoton Phase1 Age0
#############

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/DiPhotonBornPt25To250Age0.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_DiPhotonBornPt25To250Age0_${i}.out -J Phase2Upgrade_PhotonTrees_DiPhotonBornPt25To250Age0_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtuple.DiPhotonBornPt25To250Age0.${i}.root\",true\) PhotonNtuple.DiPhotonBornPt25To250Age0.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/DiPhotonBoxPt25To250Age0_STAR17_61_V1A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_DiPhotonBoxPt25To250Age0_STAR17_61_V1A_${i}.out -J Phase2Upgrade_PhotonTrees_DiPhotonBoxPt25To250Age0_STAR17_61_V1A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtuple.DiPhotonBoxPt25To250Age0.${i}.root\",true\) PhotonNtuple.DiPhotonBoxPt25To250Age0.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end


#############
# Hgg Phase1 Age0
#############

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/WZHgg-125-Age0DES_DES17_61_V5.txt `)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_WZHgg-125-Age0DES_DES17_61_V5.txt_${i}.out -J Phase2Upgrade_PhotonTrees_WZHgg-125-Age0DES_DES17_61_V5.txt_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtuple.WZHgg-125-Age0DES_DES17_61_V5.${i}.root\",true\) PhotonNtuple.WZHgg-125-Age0DES_DES17_61_V5.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end


set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/WZHgg-125-Age0_STAR17_61_V1A.txt`)
  foreach j(`seq 0 1 49`)
    echo $file " " $i " - " $j
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_WZHgg-125-Age0_STAR17_61_V1A_${i}.out -J Phase2Upgrade_PhotonTrees_WZHgg-125-Age0_STAR17_61_V1A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtuple.real.WZHgg-125-Age0_STAR17_61_V1A.${i}.${j}.root\",true,50,${j}\) PhotonNtuple.real.WZHgg-125-Age0_STAR17_61_V1A.${i}.${j}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
    sleep 0.5
  end
  @ i = $i + 1
end



####################################################
#
# Fake Photons
#
####################################################

##########################
# Summer12 QCD samples
##########################

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt30To50-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_qcd_pt30To50-START53_V7A_${i}.out -J Phase2Upgrade_PhotonTrees_qcd_pt30To50-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtupleFake.qcd_pt30To50-START53_V7A.${i}.root\",false\) PhotonNtupleFake.qcd_pt30To50-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt50To80-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_qcd_pt50To80-START53_V7A_${i}.out -J Phase2Upgrade_PhotonTrees_qcd_pt50To80-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtupleFake.qcd_pt50To80-START53_V7A.${i}.root\",false\) PhotonNtupleFake.qcd_pt50To80-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/HHToBBGG/catalog/qcd_pt80To120-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_qcd_pt80To120-START53_V7A_${i}.out -J Phase2Upgrade_PhotonTrees_qcd_pt80To120-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtupleFake.qcd_pt80To120-START53_V7A.${i}.root\",false\) PhotonNtupleFake.qcd_pt80To120-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end

##########################
# Summer12 ttHbb
##########################
set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/ttHbb-125-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 2nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_ttHbb-125-START53_${i}.out -J Phase2Upgrade_PhotonTrees_ttHbb-125-START53_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtupleFake.ttHbb-125-START53.${i}.root\",false\) PhotonNtupleFake.ttHbb-125-START53.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end


##########################
# Phase1 ttbar, ttHbb
##########################

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/ttHbb-125-Age0_STAR17_61_V1A.txt`)
  foreach j(`seq 0 1 49`)
    echo $file " " $i " - " $j
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_ttHbb-125-Age0_STAR17_61_V1A_${i}.out -J Phase2Upgrade_PhotonTrees_ttHbb-125-Age0_STAR17_61_V1A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtuple.fake.ttHbb-125-Age0_STAR17_61_V1A.${i}.${j}.root\",false,50,${j}\) PhotonNtuple.fake.ttHbb-125-Age0_STAR17_61_V1A.${i}.${j}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
    sleep 0.5
  end
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/ttHbb-125-Age0DES_DES17_61_V5.txt`)
  foreach j(`seq 0 1 49`)
    echo $file " " $i " - " $j
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_ttHbb-125-Age0DES_DES17_61_V5_${i}.out -J Phase2Upgrade_PhotonTrees_ttHbb-125-Age0DES_DES17_61_V5_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtuple.fake.ttHbb-125-Age0DES_DES17_61_V5.${i}.${j}.root\",false,50,${j}\) PhotonNtuple.fake.ttHbb-125-Age0DES_DES17_61_V5.${i}.${j}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  end
  @ i = $i + 1
end



set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/TTBARAge0.txt`)
  echo $file " " $i
  bsub -q 2nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_TTBARAge0_${i}.out -J Phase2Upgrade_PhotonTrees_TTBARAge0_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtupleFake.TTBARAge0.${i}.root\",false\) PhotonNtupleFake.TTBARAge0.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/DiPhotonBornPt25To250Age0_STAR17_61_V1A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_ttHbb-125-Age0DES_DES17_61_V5_${i}.out -J Phase2Upgrade_PhotonTrees_ttHbb-125-Age0DES_DES17_61_V5_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtupleFake.ttHbb-125-Age0DES_DES17_61_V5.${i}.root\",false\) PhotonNtupleFake.ttHbb-125-Age0DES_DES17_61_V5.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end

set i=0
foreach file(`cat /afs/cern.ch/work/s/sixie/public/releases/analysis/CMSSW_5_3_9_patch3/src/CMSAna/ObjectStudies/catalog/DiPhotonBoxPt25To250Age0_STAR17_61_V1A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Phase2Upgrade/PhotonTrees_DiPhotonBoxPt25To250Age0_STAR17_61_V1A_${i}.out -J Phase2Upgrade_PhotonTrees_DiPhotonBoxPt25To250Age0_STAR17_61_V1A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/ObjectStudies/Photons/ MakePhotonNtupleFromMC.C +\(\"root://eoscms//eos/cms/${file}\",\"PhotonNtupleFake.DiPhotonBoxPt25To250Age0_STAR17_61_V1A.${i}.root\",false\) PhotonNtupleFake.DiPhotonBoxPt25To250Age0_STAR17_61_V1A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/PhotonStudies/ntuples/jobs/
  @ i = $i + 1
end
