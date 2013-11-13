UInt_t FindFakeMuonResponseBin( double value, double bins[], UInt_t nbins) {
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


Double_t GetMuonResponseMeanPtEta(Double_t Pt, Double_t Eta) {

  Double_t ptBins[5] = {5,7,10,15,25};
  Double_t etaBins[4] = {0,1.2,2.2,2.4};


  Double_t ResponseMean[6][5] = {
    {0,0,0,0,0    },
    {0,-0.249732,-0.363671,-0.36232,0    },
    {0,-0.346771,-0.439966,-0.413834,0    },
    {0,-0.402514,-0.482566,-0.44711,0    },
    {0,-0.424105,-0.490755,-0.48053,0    },
    {0,-0.388474,-0.478037,-0.454262,0}
  };


  Int_t tmpPtBin = FindFakeMuonResponseBin( Pt , ptBins, 4);
  Int_t tmpEtaBin = FindFakeMuonResponseBin( fabs(Eta) , etaBins, 3);
  return ResponseMean[tmpPtBin][tmpEtaBin];
}


Double_t GetMuonResponseSigmaPtEta(Double_t Pt, Double_t Eta) {

  Double_t ptBins[5] = {5,7,10,15,25};
  Double_t etaBins[4] = {0,1.2,2.2,2.4};


  Double_t ResponseSigma[6][5] = {
    {0,0,0,0,0    },
    {0,0.11903,0.167462,0.186057,0    },
    {0,0.143104,0.172082,0.165416,0    },
    {0,0.155361,0.183543,0.194661,0    },
    {0,0.156805,0.182135,0.17721,0    },
    {0,0.194042,0.201404,0.24297,0}
  };


  Int_t tmpPtBin = FindFakeMuonResponseBin( Pt , ptBins, 4);
  Int_t tmpEtaBin = FindFakeMuonResponseBin( fabs(Eta) , etaBins, 3);
  return ResponseSigma[tmpPtBin][tmpEtaBin];
}
