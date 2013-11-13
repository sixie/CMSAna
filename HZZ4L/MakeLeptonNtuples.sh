##########################################################################################
#Electrons
##########################################################################################

root -l -b -q CMSAna/HZZ4L/MakeHZZEleNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V3/BACONNtuple_1125.root","ElectronNtuple.HZZ4L.1125.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZEleNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V3/BACONNtuple_100.root","ElectronNtuple.HZZ4L.100.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZEleNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V3/BACONNtuple_101.root","ElectronNtuple.HZZ4L.101.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZEleNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V3/BACONNtuple_102.root","ElectronNtuple.HZZ4L.102.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZEleNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V3/BACONNtuple_105.root","ElectronNtuple.HZZ4L.105.root")'

foreach i(0 1 2 3 4 5 6 7 8 9)
   bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HZZ4L/MakeHZZFakeEleNtuple_${i}.out -J MakeHZZFakeEleNtuple_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_HZZ4L.csh /afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/ MakeHZZFakeEleNtuple.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/BACON/V4/111/BACONNtuple_111.${i}.root\",\"FakeElectronNtuple.HZZ4L.111.${i}.root\"\) FakeElectronNtuple.HZZ4L.111.${i}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/ntuples/LeptonResponse/V4/jobs/
end
foreach i(0 1 2 3 4 5 6 7 8 9)
   bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/HZZ4L/MakeHZZFakeMuonNtuple_${i}.out -J MakeHZZFakeMuonNtuple_${i} /afs/cern.ch/work/s/sixie/public/condor/bin/runRootJob_HZZ4L.csh /afs/cern.ch/work/s/sixie/public/releases/syncForRun1LegacyPaper/CMSSW_5_3_9_patch3/src/CMSAna/HZZ4L/ MakeHZZFakeMuonNtuple.C +\(\"root://eoscms//eos/cms/store/group/phys_higgs/cmshzz4l/sixie/BACON/V4/111/BACONNtuple_111.${i}.root\",\"FakeMuonNtuple.HZZ4L.111.${i}.root\"\) FakeMuonNtuple.HZZ4L.111.${i}.root /afs/cern.ch/work/s/sixie/public/HZZ4l/ntuples/LeptonResponse/V4/jobs/
end


##########################################################################################
#Muons
##########################################################################################
root -l -b -q CMSAna/HZZ4L/MakeHZZMuonNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V3/BACONNtuple_1125.root","MuonNtuple.HZZ4L.1125.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZMuonNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V3/BACONNtuple_100.root","MuonNtuple.HZZ4L.100.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZMuonNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V3/BACONNtuple_101.root","MuonNtuple.HZZ4L.101.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZMuonNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V3/BACONNtuple_103.root","MuonNtuple.HZZ4L.103.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZMuonNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V3/BACONNtuple_105.root","MuonNtuple.HZZ4L.105.root")'

##########################################################################################
#Make Efficiency Ntuple
##########################################################################################
root -l -b -q CMSAna/HZZ4L/MakeEfficiencyMapNtuple.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/LeptonResponse/V3/ElectronNtuple.HZZ4L.combined.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/LeptonResponse/V3/MuonNtuple.HZZ4L.combined.root", "combined")'
root -l -b -q CMSAna/HZZ4L/MakeEfficiencyMapNtuple.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/LeptonResponse/V3/ElectronNtuple.HZZ4L.1125.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/LeptonResponse/V3/MuonNtuple.HZZ4L.1125.root", "1125")'
root -l -b -q CMSAna/HZZ4L/MakeEfficiencyMapNtuple.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/LeptonResponse/V3/ElectronNtuple.HZZ4L.ZZ.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/LeptonResponse/V3/MuonNtuple.HZZ4L.ZZ.root", "ZZ")'


##########################################################################################
#Make HZZ Efficiency Maps and Lepton Resolution Maps
##########################################################################################
root -l -b -q CMSAna/HZZ4L/HLL/CreateEfficiencyMap.C+'("HZZEfficiencyMap_combined.root","combined",0)'
root -l -b -q CMSAna/HZZ4L/HLL/CreateEfficiencyMap.C+'("HZZEfficiencyMap_1125.root","1125",0)'
root -l -b -q CMSAna/HZZ4L/HLL/CreateEfficiencyMap.C+'("HZZEfficiencyMap_ZZ.root","ZZ",0)'

root -l -b -q CMSAna/HZZ4L/HLL/CreateLeptonResponseMap.C+'("HZZEfficiencyMap_combined.root","combined",0)'


##########################################################################################
#Make Fake rates
##########################################################################################
root -l -b -q CMSAna/HZZ4L/HLL/CreateElectronFakeRateMap.C+'("/afs/cern.ch/work/s/sixie/public/HZZ4l/ntuples/LeptonResponse/V4/FakeElectronNtuple.HZZ4L.111.root","ele",0)'
root -l -b -q CMSAna/HZZ4L/HLL/CreateElectronFakeRateMap.C+'("/afs/cern.ch/work/s/sixie/public/HZZ4l/ntuples/LeptonResponse/V4/FakeElectronNtuple.HZZ4L.111.root","ele",1)'

root -l -b -q CMSAna/HZZ4L/HLL/CreateMuonFakeRateMap.C+'("/afs/cern.ch/work/s/sixie/public/HZZ4l/ntuples/LeptonResponse/V4/FakeMuonNtuple.HZZ4L.111.root","mu",0)'
root -l -b -q CMSAna/HZZ4L/HLL/CreateMuonFakeRateMap.C+'("/afs/cern.ch/work/s/sixie/public/HZZ4l/ntuples/LeptonResponse/V4/FakeMuonNtuple.HZZ4L.111.root","mu",1)'
