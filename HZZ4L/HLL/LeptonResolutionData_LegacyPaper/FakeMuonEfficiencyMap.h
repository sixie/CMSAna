UInt_t FindFakeMuonEfficiencyBin( double value, double bins[], UInt_t nbins) {
  UInt_t nbinboundaries = nbins+1;
  UInt_t bin = 0;
  for (uint i=0; i < nbinboundaries; ++i) {
    if (i < nbinboundaries-1) {
      if (value >= bins[i] && value < bins[i+1]) {
        bin = i+1;
        break;
      }
    } else if (i == nbinboundaries-1) {
      if (value >= bins[i]) {
        bin = nbinboundaries;
        break;
      }
    }    
  }
  return bin;
}


Double_t GetFakeMuonEfficiencyPtEta(Double_t Pt, Double_t Eta) {

  Double_t ptBins[8] = {5,7,10,15,20,25,30,35};
  Double_t etaBins[4] = {0,1.2,2.2,2.4};


  Double_t Efficiency[9][5] = {
    {0,6.34246e-05,0.000379899,0.000366359,4.41215e-05    },
    {0,0.000179326,0.000424109,0.000375324,7.84913e-05    },
    {0,0.000256191,0.00044567,0.000351584,7.53574e-05    },
    {0,0.000327362,0.00050925,0.000468659,7.41818e-05    },
    {0,0.000336716,0.00048003,0.000518476,8.56091e-05    },
    {0,0.00046874,0.000605364,0.000471846,0.000175346    },
    {0,0.000470612,0.000679046,0.000577182,0.000144009    },
    {0,0.000601961,0.000583475,0.00133537,0.000215332    },
    {0,0.000680583,0.000971012,0.00134601,0.000269505}
  };


  Int_t tmpPtBin = FindFakeMuonEfficiencyBin( Pt , ptBins, 7);
  Int_t tmpEtaBin = FindFakeMuonEfficiencyBin( fabs(Eta) , etaBins, 3);
  return Efficiency[tmpPtBin][tmpEtaBin];
}
