##########################################################################################
#Produce Smeared Gen-Level samples
##########################################################################################
foreach j(`seq 0 1 99`)
  echo  $j
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HZZ4L/FastSim_1125_${j}.out -J HZZ4L_FastSim_1125_${j} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_HZZ4L.csh /afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/ MakeHZZEventNtupleUsingFastsim.C +\(\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/BACON/V3/BACONNtuple_1125.root\",\"HZZEventNtuple_Fastsim_1125.${j}.root\",100,${j}\) HZZEventNtuple_Fastsim_1125.${j}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/ntuples/HZZEvent/jobs/
end

foreach j(`seq 0 1 99`)
  echo  $j
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HZZ4L/FastSim_102_${j}.out -J HZZ4L_FastSim_102_${j} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_HZZ4L.csh /afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/ MakeHZZEventNtupleUsingFastsim.C +\(\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/BACON/V3/BACONNtuple_102.root\",\"HZZEventNtuple_Fastsim_102.${j}.root\",100,${j}\) HZZEventNtuple_Fastsim_102.${j}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/ntuples/HZZEvent/jobs/
end

foreach j(`seq 0 1 99`)
  echo  $j
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HZZ4L/FastSim_103_${j}.out -J HZZ4L_FastSim_103_${j} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_HZZ4L.csh /afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/ MakeHZZEventNtupleUsingFastsim.C +\(\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/BACON/V3/BACONNtuple_103.root\",\"HZZEventNtuple_Fastsim_103.${j}.root\",100,${j}\) HZZEventNtuple_Fastsim_103.${j}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/ntuples/HZZEvent/jobs/
end

foreach j(`seq 0 1 99`)
  echo  $j
  bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HZZ4L/FastSim_105_${j}.out -J HZZ4L_FastSim_105_${j} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_HZZ4L.csh /afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/ MakeHZZEventNtupleUsingFastsim.C +\(\"/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/BACON/V3/BACONNtuple_105.root\",\"HZZEventNtuple_Fastsim_105.${j}.root\",100,${j}\) HZZEventNtuple_Fastsim_105.${j}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/ntuples/HZZEvent/jobs/
end


##########################################################################################
#Produce fullsim samples
##########################################################################################
root -l -b -q /afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/MakeHZZEventNtuple.C+'("1125",1125)'
root -l -b -q /afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/MakeHZZEventNtuple.C+'("102",102)'
root -l -b -q /afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/MakeHZZEventNtuple.C+'("103",103)'
root -l -b -q /afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/MakeHZZEventNtuple.C+'("105",105)'


##########################################################################################
#Do fullsim vs smeared comparison
##########################################################################################
##HZZ signal
root -l -b -q CMSAna/HZZ4L/ValidateFastSim.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_1125.V3.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_Fastsim_1125.V3.root","1125")'
root -l -b -q CMSAna/HZZ4L/ValidateFastSim.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_1125.V3.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_Fastsim_1125.V3.root","1125_Normalized",true)'

##ZZ bkg
root -l -b -q CMSAna/HZZ4L/ValidateFastSim.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_102.V3.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_Fastsim_102.V3.root","102")'
root -l -b -q CMSAna/HZZ4L/ValidateFastSim.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_103.V3.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_Fastsim_103.V3.root","103")'
root -l -b -q CMSAna/HZZ4L/ValidateFastSim.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_105.V3.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_Fastsim_105.V3.root","105")'



##Gen-Level comparisons
root -l -b -q CMSAna/HZZ4L/ValidateFastSim_GenLevelMomenta.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_1125.V3.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_Fastsim_1125.V3.root","1125_GenLevelMomenta")'
root -l -b -q CMSAna/HZZ4L/ValidateFastSim_GenLevelMomenta.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_102.V3.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_Fastsim_102.V3.root","102_GenLevelMomenta")'
root -l -b -q CMSAna/HZZ4L/ValidateFastSim_GenLevelMomenta.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_103.V3.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_Fastsim_103.V3.root","103_GenLevelMomenta")'
root -l -b -q CMSAna/HZZ4L/ValidateFastSim_GenLevelMomenta.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_105.V3.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/HZZEvent/HZZEventNtuple_Fastsim_105.V3.root","105_GenLevelMomenta")'


