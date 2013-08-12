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
# Signal
#############
set i=0
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/catalog/HHtoBBGG-8tev-START53_V7A.txt `)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HHToBBGG/HHToBBGGSelection_HHtoBBGG-8tev-START53_V7A_${i}.out -J HHToBBGGSelection_HHtoBBGG-8tev-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/ HHToBBGGSelection.C +\(\"root://eoscms//eos/cms/${file}\",\"HHToBBGGNtuple.HHtoBBGG-8tev-START53_V7A.${i}.root\",1\) HHToBBGGNtuple.HHtoBBGG-8tev-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/ntuples/jobs/
  @ i = $i + 1
end


set i=0
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/catalog/HHtoBBGG-14tev-START53_V7A.txt `)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HHToBBGG/HHToBBGGSelection_HHtoBBGG-14tev-START53_V7A_${i}.out -J HHToBBGGSelection_HHtoBBGG-14tev-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/ HHToBBGGSelection.C +\(\"root://eoscms//eos/cms/${file}\",\"HHToBBGGNtuple.HHtoBBGG-14tev-START53_V7A.${i}.root\",1\) HHToBBGGNtuple.HHtoBBGG-14tev-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/ntuples/jobs/
  @ i = $i + 1
end

#############
# ttH
#############

set i=0
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/catalog/ttHgg-125-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HHToBBGG/HHToBBGGSelection_ttHgg-125-START53_V7A_${i}.out -J HHToBBGGSelection_ttHgg-125-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/ HHToBBGGSelection.C +\(\"root://eoscms//eos/cms/${file}\",\"HHToBBGGNtuple.ttHgg-125-START53_V7A.${i}.root\",2\) HHToBBGGNtuple.ttHgg-125-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/ntuples/jobs/
  @ i = $i + 1
end

#############
# ZHgg
#############
set i=0
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/catalog/ZHgg-125-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HHToBBGG/HHToBBGGSelection_ZHgg-125-START53_V7A_${i}.out -J HHToBBGGSelection_ZHgg-125-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/ HHToBBGGSelection.C +\(\"root://eoscms//eos/cms/${file}\",\"HHToBBGGNtuple.ZHgg-125-START53_V7A.${i}.root\",3\) HHToBBGGNtuple.ZHgg-125-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/ntuples/jobs/
  @ i = $i + 1
end


#############
# ggHgg
#############

set i=0
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/catalog/ggHgg-125-START53_V7A.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HHToBBGG/HHToBBGGSelection_ggHgg-125-START53_V7A_${i}.out -J HHToBBGGSelection_ggHgg-125-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/ HHToBBGGSelection.C +\(\"root://eoscms//eos/cms/${file}\",\"HHToBBGGNtuple.ggHgg-125-START53_V7A.${i}.root\",4\) HHToBBGGNtuple.ggHgg-125-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/ntuples/jobs/
  sleep 0.25
  @ i = $i + 1
end



#############
# Diphoton+BB : bb+gg
#############


set i=0
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/catalog/DiPhotonBB_M60To200_14TeV.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HHToBBGG/HHToBBGGSelection_BBGG_DiPhotonBB_M60To200_14TeV_${i}.out -J HHToBBGGSelection_BBGG_DiPhotonBB_M60To200_14TeV_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/ HHToBBGGSelection_BBGG.C +\(\"root://eoscms//eos/cms/${file}\",\"HHToBBGGNtuple.DiPhotonBB_M60To200_14TeV-BBGG.${i}.root\",6\) HHToBBGGNtuple.DiPhotonBB_M60To200_14TeV-BBGG.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/ntuples/jobs/
  sleep 0.25
  @ i = $i + 1
end

#############
# Diphoton+Jets : gg+jj
#############


set i=0
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/catalog/DiPhotonJJ_M60To200_14TeV.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HHToBBGG/HHToBBGGSelection_DiPhotonJJ_M60To200_14TeV_${i}.out -J HHToBBGGSelection_DiPhotonJJ_M60To200_14TeV_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/ HHToBBGGSelectionMistags.C +\(\"root://eoscms//eos/cms/${file}\",\"HHToBBGGNtuple.DiPhotonJJ_M60To200_14TeV.${i}.root\",7\) HHToBBGGNtuple.DiPhotonJJ_M60To200_14TeV.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/ntuples/jobs/
  sleep 0.25
  @ i = $i + 1
end


#############
# ttbar
#############

set i=0
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/catalog/ttjll-START53_V7A.txt `)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HHToBBGG/HHToBBGGSelection_ttjll-START53_V7A_${i}.out -J HHToBBGGSelection_ttjll-START53_V7A_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/ HHToBBGGSelection_TTBAR.C +\(\"root://eoscms//eos/cms/${file}\",\"HHToBBGGNtuple.ttjll-START53_V7A.${i}.root\",5\) HHToBBGGNtuple.ttjll-START53_V7A.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/ntuples/jobs/
  sleep 0.25
  @ i = $i + 1
end


#############
# BBH->GG
#############


set i=0
foreach file(`cat /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/catalog/bbH_M60To200_14TeV.txt`)
  echo $file " " $i
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HHToBBGG/HHToBBGGSelection_BBHGG_bbH_M60To200_14TeV_${i}.out -J HHToBBGGSelection_BBHGG_bbH_M60To200_14TeV_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_upgrade.csh /afs/cern.ch/user/s/sixie/CMSSW_upgrade/src/CMSAna/HHToBBGG/ HHToBBGGSelection_BBHGG.C +\(\"root://eoscms//eos/cms/${file}\",\"HHToBBGGNtuple.bbH_M60To200_14TeV.${i}.root\",20\) HHToBBGGNtuple.bbH_M60To200_14TeV.${i}.root /afs/cern.ch/work/s/sixie/public/Phase2Upgrade/HHToBBGG/ntuples/jobs/
  sleep 0.25
  @ i = $i + 1
end
