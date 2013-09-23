##########################################################################################
#Electrons
##########################################################################################

root -l -b -q CMSAna/HZZ4L/MakeHZZEleNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V2/BACONNtuple_1125.root","ElectronNtuple.HZZ4L.1125.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZEleNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V2/BACONNtuple_100.root","ElectronNtuple.HZZ4L.100.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZEleNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V2/BACONNtuple_101.root","ElectronNtuple.HZZ4L.101.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZEleNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V2/BACONNtuple_102.root","ElectronNtuple.HZZ4L.102.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZEleNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V2/BACONNtuple_105.root","ElectronNtuple.HZZ4L.105.root")'


##########################################################################################
#Muons
##########################################################################################
root -l -b -q CMSAna/HZZ4L/MakeHZZMuonNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V2/BACONNtuple_1125.root","MuonNtuple.HZZ4L.1125.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZMuonNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V2/BACONNtuple_100.root","MuonNtuple.HZZ4L.100.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZMuonNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V2/BACONNtuple_101.root","MuonNtuple.HZZ4L.101.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZMuonNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V2/BACONNtuple_103.root","MuonNtuple.HZZ4L.103.root")'
root -l -b -q CMSAna/HZZ4L/MakeHZZMuonNtuple.C+'("/data1/sixie/ntuples/BACON/HZZ4L/V2/BACONNtuple_105.root","MuonNtuple.HZZ4L.105.root")'

##########################################################################################
#Make Efficiency Ntuple
##########################################################################################
root -l -b -q CMSAna/HZZ4L/MakeEfficiencyMapNtuple.C+'("/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/LeptonResponse/ElectronNtuple.HZZ4L.combined.root","/afs/cern.ch/user/s/sixie/work/public/HZZ4l/ntuples/LeptonResponse/MuonNtuple.HZZ4L.combined.root", "combined")'


##########################################################################################
#Make HZZ Efficiency Maps and Lepton Resolution Maps
##########################################################################################
root -l -b -q CMSAna/HZZ4L/HLL/CreateEfficiencyMap.C+'("HZZEfficiencyMap_combined.root","combined",0)'
root -l -b -q CMSAna/HZZ4L/HLL/CreateLeptonResponseMap.C+'("HZZEfficiencyMap_combined.root","combined",0)'
