//


#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "CMSAna/HHToBBGG/interface/HHToBBGGEventTree.hh"
#include "CMSAna/Utils/interface/SimpleTable.h"


//--------------------------------------------------------------------------------------------------
// Get Total Number of Events in the sample
//--------------------------------------------------------------------------------------------------
Double_t getNormalizationWeight(string filename, string datasetName) {
  // Get Normalization Weight Factor

  //Get Number of Events in the Sample
  TFile *file = new TFile(filename.c_str(),"READ");
  if (!file) {
    cout << "Could not open file " << filename << endl;
    return 0;
  }

  //TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  //if (!dir) {
  //  cout << "Could not find directory AnaFwkMod"
  //       << " in file " << filename << endl;
  //  delete file;
  //  return 0;
  //}

  TH1F *hist = (TH1F*) file->Get("NEvents");
  if (!hist) {
    cout << "Could not find histogram hDAllEvents in directory AnaFwkMod"
         << " in file " << filename << endl;
    //delete dir;
    file->Close();
    delete file;
    return 0;
  }
  Double_t NEvents = hist->Integral();
  cout << "Original events in the sample: " << NEvents << endl;

  //Get CrossSection
  cmsana::SimpleTable xstab("$CMSSW_BASE/src/CMSAna/HHToBBGG/data/xs.dat");
  Double_t CrossSection = xstab.Get(datasetName.c_str());  
  Double_t Weight = CrossSection / NEvents;
  // weight for data is always 1 (-1 to make a trick for fakes)
  if(CrossSection < 0) Weight = -1.0;

  cout << "Cross Section = " << CrossSection << " for dataset " << datasetName << "\n";
  cout << "Events get weight: " << Weight << "\n";

  //delete dir;
  file->Close();
  delete file;
  return Weight;

}

//*************************************************************************************************
//Main part of the macro
//*************************************************************************************************
void NormalizeNtuple(const string InputFilename, 
                     const string datasetName, 
                     const string OutputFilename) {
  
  cmsana::HHToBBGGEventTree event;
  event.LoadTree(InputFilename.c_str());
  event.InitTree();

  Double_t normalizationWeight = getNormalizationWeight(InputFilename, datasetName);

  //*************************************************************************************************
  //Create new normalized tree
  //*************************************************************************************************
  TFile *outputFile = new TFile(OutputFilename.c_str(), "RECREATE");
  outputFile->cd();

  TTree *normalizedTree = event.tree_->CloneTree(0);  
  cout << "Events in the ntuple: " << event.tree_->GetEntries() << endl;
  
  for (int n=0;n<event.tree_->GetEntries();n++) { 
    event.tree_->GetEntry(n);
    if     (normalizationWeight < 0){
      event.weight = event.weight * 1.0;
    } else {
      event.weight = event.weight * normalizationWeight *  1.0;
    }
    normalizedTree->Fill(); 
  }

  normalizedTree->Write();
  outputFile->Close();
}
