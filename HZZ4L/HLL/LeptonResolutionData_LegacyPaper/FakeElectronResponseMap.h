UInt_t FindFakeElectronResponseBin( double value, double bins[], UInt_t nbins) {
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


Double_t GetElectronResponseMeanPtEta(Double_t Pt, Double_t Eta) {

  Double_t ptBins[5] = {5,7,10,15,25};
  Double_t etaBins[3] = {0,1.5,2.5};


  Double_t ResponseMean[6][4] = {
    {-5.01215,-5.01215,-5.01215,-5.01215    },
    {-5.01215,-0.267489,-0.208188,-5.01215    },
    {-5.01215,-0.237655,-0.255802,-5.01215    },
    {-5.01215,-0.241834,-0.245421,-5.01215    },
    {-5.01215,-0.247748,-0.198943,-5.01215    },
    {-5.01215,-0.24667,-0.199304,-5.01215}
  };


  Int_t tmpPtBin = FindFakeElectronResponseBin( Pt , ptBins, 4);
  Int_t tmpEtaBin = FindFakeElectronResponseBin( fabs(Eta) , etaBins, 2);
  return ResponseMean[tmpPtBin][tmpEtaBin];
}


Double_t GetElectronResponseSigmaPtEta(Double_t Pt, Double_t Eta) {

  Double_t ptBins[5] = {5,7,10,15,25};
  Double_t etaBins[3] = {0,1.5,2.5};


  Double_t ResponseSigma[6][4] = {
    {0.05,0.05,0.05,0.05    },
    {0.05,0.166146,0.187806,0.05    },
    {0.05,0.190038,0.205906,0.05    },
    {0.05,0.165298,0.185579,0.05    },
    {0.05,0.170588,0.192198,0.05    },
    {0.05,0.156883,0.156483,0.05}
  };


  Int_t tmpPtBin = FindFakeElectronResponseBin( Pt , ptBins, 4);
  Int_t tmpEtaBin = FindFakeElectronResponseBin( fabs(Eta) , etaBins, 2);
  return ResponseSigma[tmpPtBin][tmpEtaBin];
}
