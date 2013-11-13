UInt_t FindFakeElectronEfficiencyBin( double value, double bins[], UInt_t nbins) {
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


Double_t GetFakeElectronEfficiencyPtEta(Double_t Pt, Double_t Eta) {

  Double_t ptBins[8] = {5,7,10,15,20,25,30,35};
  Double_t etaBins[4] = {0,1,2,2.5};


  Double_t Efficiency[9][5] = {
    {0,0,0,0,0    },
    {0,1.68221e-05,5.841e-05,8.95084e-05,0    },
    {0,5.10577e-05,0.00014004,0.000143058,0    },
    {0,0.000243361,0.000461522,0.000433131,0    },
    {0,0.000725843,0.00104842,0.00090709,0    },
    {0,0.00100456,0.00165935,0.00175351,0    },
    {0,0.00124005,0.00211075,0.00176218,0    },
    {0,0.00132031,0.00222611,0.00256162,0    },
    {0,0.00142783,0.00266223,0.0024952,0}
  };


  //find nearest bins
  Int_t nearestPtBin = FindFakeElectronEfficiencyBin( Pt , ptBins, 7);
  Int_t nearestEtaBin = FindFakeElectronEfficiencyBin( fabs(Eta) , etaBins, 3);

  return Efficiency[nearestPtBin][nearestEtaBin];

}



Double_t GetFakeElectronEfficiencyInterpolatedPtEta(Double_t Pt, Double_t Eta) {

  Double_t ptBins[8] = {5,7,10,15,20,25,30,35};
  Double_t etaBins[4] = {0,1,2,2.5};


  Double_t Efficiency[9][5] = {
    {0,0,0,0,0    },
    {0,1.68221e-05,5.841e-05,8.95084e-05,0    },
    {0,5.10577e-05,0.00014004,0.000143058,0    },
    {0,0.000243361,0.000461522,0.000433131,0    },
    {0,0.000725843,0.00104842,0.00090709,0    },
    {0,0.00100456,0.00165935,0.00175351,0    },
    {0,0.00124005,0.00211075,0.00176218,0    },
    {0,0.00132031,0.00222611,0.00256162,0    },
    {0,0.00142783,0.00266223,0.0024952,0}
  };


  //find nearest bins
  Int_t nearestPtBin = FindFakeElectronEfficiencyBin( Pt , ptBins, 7);
  Int_t nearestEtaBin = FindFakeElectronEfficiencyBin( fabs(Eta) , etaBins, 3);

  return Efficiency[nearestPtBin][nearestEtaBin];

  //find the nearest bin centers
  double NearestBinCenter_Pt = (ptBins[nearestPtBin-1] + ptBins[nearestPtBin])/2.0;
  double NearestBinCenter_Eta = (etaBins[nearestEtaBin-1] + etaBins[nearestEtaBin])/2.0;

  //find next to nearest bins
  Int_t nextToNearestPtBin = -2;
  Int_t nextToNearestEtaBin = -2;
  if (Pt > NearestBinCenter_Pt) nextToNearestPtBin = nearestPtBin+1;
  else nextToNearestPtBin = nearestPtBin-1;
  if (fabs(Eta) > NearestBinCenter_Eta) nextToNearestEtaBin = nearestEtaBin+1;
  else nextToNearestEtaBin = nearestEtaBin-1;
  
  //find next to nearest bin centers
  double NextToNearestBinCenter_Pt = -999;(ptBins[nextToNearestPtBin-1] + ptBins[nextToNearestPtBin])/2.0;
  double NextToNearestBinCenter_Eta = -999; (etaBins[nextToNearestEtaBin-1] + etaBins[nextToNearestEtaBin])/2.0;
  if (nextToNearestPtBin >= 1) {
    NextToNearestBinCenter_Pt = (ptBins[nextToNearestPtBin-1] + ptBins[nextToNearestPtBin])/2.0;
  }
  if (nextToNearestEtaBin >= 1) {
    NextToNearestBinCenter_Eta = (etaBins[nextToNearestEtaBin-1] + etaBins[nextToNearestEtaBin])/2.0;  
  }

  //4 bin centers are:
  //(NearestBinCenter_Pt, NearestBinCenter_Eta)
  //(NearestBinCenter_Pt, NextToNearestBinCenter_Eta)
  //(NextToNearestBinCenter_Pt, NearestBinCenter_Eta)
  //(NextToNearestBinCenter_Pt, NextToNearestBinCenter_Eta)

  double EffInterpolated = 0;

//   cout << "\n\n";
//   cout << "Electron: " << Pt << " " << Eta << "\n";
//   cout << "nearestPtBin = " << nearestPtBin << " , nearestEtaBin = " << nearestEtaBin << "\n";
//   cout << "NearestBinCenter_Pt = " << NearestBinCenter_Pt << " , NearestBinCenter_Eta = " << NearestBinCenter_Eta << "\n";
//   cout << "nextToNearestPtBin = " << nextToNearestPtBin << " , nextToNearestEtaBin = " << nextToNearestEtaBin << "\n";
//   cout << "NextToNearestBinCenter_Pt = " << NextToNearestBinCenter_Pt << " , NextToNearestBinCenter_Eta = " << NextToNearestBinCenter_Eta << "\n";


  //*************************************************
  //Overflow Eta bins
  //*************************************************
  if (nearestEtaBin == 0 || nearestEtaBin == 4) {
    EffInterpolated = 0; 
  } else if (nearestEtaBin == 1) {


    //*************************************************
    //On left-most Eta bin
    //*************************************************
//     cout << "left-most eta bin\n";

    //#################################################
    // underflow Pt bin
    //#################################################
    if (nearestPtBin == 0) {
      EffInterpolated = Efficiency[nearestPtBin][nearestEtaBin]; //it's 0
    } else if (nearestPtBin == 1) {
      //#################################################
      // left-most Pt Bin
      //#################################################
//      cout << "left most pt bin\n";
      if (fabs(Eta) < NearestBinCenter_Eta && Pt < NearestBinCenter_Pt ) {
        //a constant
        EffInterpolated =Efficiency[nearestPtBin][nearestEtaBin]; 
//         cout << "left-left: constant \n";
      } else if ( fabs(Eta) < NearestBinCenter_Eta && Pt >= NearestBinCenter_Pt) {
        //do 1D linear interpolation in pt direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Pt-NextToNearestBinCenter_Pt))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(Pt - NextToNearestBinCenter_Pt)
          + Efficiency[nextToNearestPtBin][nearestEtaBin]*fabs(Pt - NearestBinCenter_Pt)
          );      
//         cout << "left-right: 1D interpolation in pT : " << "(1.0/fabs(" << NearestBinCenter_Pt << " - " << NextToNearestBinCenter_Pt << "))*(" << Efficiency[nearestPtBin][nearestEtaBin] << " * fabs(" << Pt << " - " <<  NextToNearestBinCenter_Pt << ") + " << Efficiency[nextToNearestPtBin][nearestEtaBin] << "*fabs(" << Pt << " - " <<  NearestBinCenter_Pt << "))" << "\n";
      } else if ( fabs(Eta) >= NearestBinCenter_Eta && Pt < NearestBinCenter_Pt) {
        //do 1D linear interpolation in eta direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Eta-NextToNearestBinCenter_Eta))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(fabs(Eta) - NextToNearestBinCenter_Eta)
          + Efficiency[nearestPtBin][nextToNearestEtaBin]*fabs(fabs(Eta) - NearestBinCenter_Eta)
          );
      } else {
        //do 2D interpolation
        EffInterpolated = 
          (1.0 / (fabs(NextToNearestBinCenter_Pt - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );
      }
    } else if (nearestPtBin == 7) {
      //#################################################
      // right-most Pt Bin
      //#################################################
      if (fabs(Eta) < NearestBinCenter_Eta && Pt >= NearestBinCenter_Pt ) {
        //do 1D linear interpolation in pt direction up to Pt bin edge
        EffInterpolated = (1.0/fabs(NearestBinCenter_Pt-ptBins[nearestPtBin]))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(Pt - ptBins[nearestPtBin])
          + Efficiency[nextToNearestPtBin][nearestEtaBin]*fabs(Pt - NearestBinCenter_Pt)
          );
      } else if ( fabs(Eta) < NearestBinCenter_Eta && Pt < NearestBinCenter_Pt) {
        //do 1D linear interpolation in pt direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Pt-NextToNearestBinCenter_Pt))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(Pt - NextToNearestBinCenter_Pt)
          + Efficiency[nextToNearestPtBin][nearestEtaBin]*fabs(Pt - NearestBinCenter_Pt)
          );      
      } else if ( fabs(Eta) >= NearestBinCenter_Eta && Pt >= NearestBinCenter_Pt) {
        //do 2D interpolation up to Pt bin edge 
        EffInterpolated = 
          (1.0 / (fabs(ptBins[nearestPtBin] - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(ptBins[nearestPtBin] - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(ptBins[nearestPtBin] - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );
      } else {
        //do 2D interpolation
        EffInterpolated = 
          (1.0 / (fabs(NextToNearestBinCenter_Pt - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );
      }
    } else if (nearestPtBin == 8) {
      //#################################################
      // overflow Pt Bin
      //#################################################
      if (fabs(Eta) < NearestBinCenter_Eta) {
        //a constant
        EffInterpolated = Efficiency[nearestPtBin][nearestEtaBin]; 
      } else {
        //do 1D linear interpolation in eta direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Eta-NextToNearestBinCenter_Eta))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(fabs(Eta) - NextToNearestBinCenter_Eta)
          + Efficiency[nearestPtBin][nextToNearestEtaBin]*fabs(fabs(Eta) - NearestBinCenter_Eta)
          );        
      }      
    } else {
      //#################################################
      // regular Pt Bin
      //#################################################
      if (fabs(Eta) < NearestBinCenter_Eta) {
        //do 1D linear interpolation in pt direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Pt-NextToNearestBinCenter_Pt))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(Pt - NextToNearestBinCenter_Pt)
          + Efficiency[nextToNearestPtBin][nearestEtaBin]*fabs(Pt - NearestBinCenter_Pt)
          );       
      } else {
        //do 2D interpolation
        EffInterpolated = 
          (1.0 / (fabs(NextToNearestBinCenter_Pt - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );        
      }
    }

  } else if (nearestEtaBin == 3) {
    //*************************************************
    //On right-most Eta bin
    //*************************************************

    //#################################################
    // underflow Pt bin
    //#################################################
    if (nearestPtBin == 0) {
      EffInterpolated = Efficiency[nearestPtBin][nearestEtaBin]; //it's 0
    } else if (nearestPtBin == 1) {
      //#################################################
      // left-most Pt Bin
      //#################################################
      if (fabs(Eta) >= NearestBinCenter_Eta && Pt < NearestBinCenter_Pt ) {
        //a constant
        EffInterpolated =Efficiency[nearestPtBin][nearestEtaBin]; 
      } else if ( fabs(Eta) >= NearestBinCenter_Eta && Pt >= NearestBinCenter_Pt) {
        //do 1D linear interpolation in pt direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Pt-NextToNearestBinCenter_Pt))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(Pt - NextToNearestBinCenter_Pt)
          + Efficiency[nextToNearestPtBin][nearestEtaBin]*fabs(Pt - NearestBinCenter_Pt)
          );      
      } else if ( fabs(Eta) < NearestBinCenter_Eta && Pt < NearestBinCenter_Pt) {
        //do 1D linear interpolation in eta direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Eta-NextToNearestBinCenter_Eta))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(fabs(Eta) - NextToNearestBinCenter_Eta)
          + Efficiency[nearestPtBin][nextToNearestEtaBin]*fabs(fabs(Eta) - NearestBinCenter_Eta)
          );
      } else {
        //do 2D interpolation
        EffInterpolated = 
          (1.0 / (fabs(NextToNearestBinCenter_Pt - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );
      }
    } else if (nearestPtBin == 7) {
      //#################################################
      // right-most Pt Bin
      //#################################################
      if (fabs(Eta) >= NearestBinCenter_Eta && Pt >= NearestBinCenter_Pt ) {
        //do 1D linear interpolation in pt direction up to Pt bin edge
        EffInterpolated = (1.0/fabs(NearestBinCenter_Pt-ptBins[nearestPtBin]))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(Pt - ptBins[nearestPtBin])
          + Efficiency[nextToNearestPtBin][nearestEtaBin]*fabs(Pt - NearestBinCenter_Pt)
          );
      } else if ( fabs(Eta) >= NearestBinCenter_Eta && Pt < NearestBinCenter_Pt) {
        //do 1D linear interpolation in pt direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Pt-NextToNearestBinCenter_Pt))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(Pt - NextToNearestBinCenter_Pt)
          + Efficiency[nextToNearestPtBin][nearestEtaBin]*fabs(Pt - NearestBinCenter_Pt)
          );      
      } else if ( fabs(Eta) < NearestBinCenter_Eta && Pt >= NearestBinCenter_Pt) {
        //do 2D interpolation up to Pt bin edge 
        EffInterpolated = 
          (1.0 / (fabs(ptBins[nearestPtBin] - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(ptBins[nearestPtBin] - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(ptBins[nearestPtBin] - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );
      } else {
        //do 2D interpolation
        EffInterpolated = 
          (1.0 / (fabs(NextToNearestBinCenter_Pt - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );
      }
    } else if (nearestPtBin == 8) {
      //#################################################
      // overflow Pt Bin
      //#################################################
      if (fabs(Eta) >= NearestBinCenter_Eta) {
        //a constant
        EffInterpolated = Efficiency[nearestPtBin][nearestEtaBin]; 
      } else {
        //do 1D linear interpolation in eta direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Eta-NextToNearestBinCenter_Eta))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(fabs(Eta) - NextToNearestBinCenter_Eta)
          + Efficiency[nearestPtBin][nextToNearestEtaBin]*fabs(fabs(Eta) - NearestBinCenter_Eta)
          );        
      }      
    } else {
      //#################################################
      // regular Pt Bin
      //#################################################
      if (fabs(Eta) >= NearestBinCenter_Eta) {
        //do 1D linear interpolation in pt direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Pt-NextToNearestBinCenter_Pt))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(Pt - NextToNearestBinCenter_Pt)
          + Efficiency[nextToNearestPtBin][nearestEtaBin]*fabs(Pt - NearestBinCenter_Pt)
          );       
      } else {
        //do 2D interpolation
        EffInterpolated = 
          (1.0 / (fabs(NextToNearestBinCenter_Pt - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );        
      }
    }
    
  } else {
    //*************************************************
    //Regular Eta bins
    //*************************************************

    //#################################################
    // underflow Pt bin
    //#################################################
    if (nearestPtBin == 0) {
      EffInterpolated = Efficiency[nearestPtBin][nearestEtaBin]; //it's 0
    } else if (nearestPtBin == 1) {

      //#################################################
      // left-most Pt Bin
      //#################################################
      if (Pt < NearestBinCenter_Pt ) {
        //do 1D linear interpolation in eta direction
        EffInterpolated = (1.0/fabs(NearestBinCenter_Eta-NextToNearestBinCenter_Eta))*(
          Efficiency[nearestPtBin][nearestEtaBin]*fabs(fabs(Eta) - NextToNearestBinCenter_Eta)
          + Efficiency[nearestPtBin][nextToNearestEtaBin]*fabs(fabs(Eta) - NearestBinCenter_Eta)
          );
      } else {
        //do 2D interpolation
        EffInterpolated = 
          (1.0 / (fabs(NextToNearestBinCenter_Pt - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );
      }
      
    } else if (nearestPtBin == 7) {

      //#################################################
      // right-most Pt Bin
      //#################################################
      if (Pt >= NearestBinCenter_Pt ) {
        //do 2D interpolation up to Pt bin edge 
        EffInterpolated = 
          (1.0 / (fabs(ptBins[nearestPtBin] - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(ptBins[nearestPtBin] - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(ptBins[nearestPtBin] - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );
      } else {
        //do 2D interpolation
        EffInterpolated = 
          (1.0 / (fabs(NextToNearestBinCenter_Pt - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
          (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
           + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
            );
      }

    } else if (nearestPtBin == 8) {

      //#################################################
      // overflow Pt Bin
      //#################################################
      //do 1D linear interpolation in eta direction
      EffInterpolated = (1.0/fabs(NearestBinCenter_Eta-NextToNearestBinCenter_Eta))*(
        Efficiency[nearestPtBin][nearestEtaBin]*fabs(fabs(Eta) - NextToNearestBinCenter_Eta)
        + Efficiency[nearestPtBin][nextToNearestEtaBin]*fabs(fabs(Eta) - NearestBinCenter_Eta)
        );        

    } else {

      //#################################################
      // regular Pt Bin
      //#################################################
      //do 2D interpolation
      EffInterpolated = 
        (1.0 / (fabs(NextToNearestBinCenter_Pt - NearestBinCenter_Pt)*fabs(NextToNearestBinCenter_Eta-NearestBinCenter_Eta)))*
        (Efficiency[nearestPtBin][nearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
         + Efficiency[nearestPtBin][nextToNearestEtaBin] * (fabs(NextToNearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
         + Efficiency[nextToNearestPtBin][nearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NextToNearestBinCenter_Eta-fabs(Eta)) )
         + Efficiency[nextToNearestPtBin][nextToNearestEtaBin] * (fabs(NearestBinCenter_Pt - Pt)*fabs(NearestBinCenter_Eta-fabs(Eta)) )
          );        
    }
    

  }
    
  //cout << "EffInterpolated = " << EffInterpolated << "\n";

  if (EffInterpolated != EffInterpolated) {
    cout << "\n\n";
    cout << "Electron: " << Pt << " " << Eta << "\n";
    cout << "nearestPtBin = " << nearestPtBin << " , nearestEtaBin = " << nearestEtaBin << "\n";
    cout << "NearestBinCenter_Pt = " << NearestBinCenter_Pt << " , NearestBinCenter_Eta = " << NearestBinCenter_Eta << "\n";
    cout << "nextToNearestPtBin = " << nextToNearestPtBin << " , nextToNearestEtaBin = " << nextToNearestEtaBin << "\n";
    cout << "NextToNearestBinCenter_Pt = " << NextToNearestBinCenter_Pt << " , NextToNearestBinCenter_Eta = " << NextToNearestBinCenter_Eta << "\n";
    cout << "EffInterpolated = " << EffInterpolated << "\n";  
  }


  return EffInterpolated;
}
