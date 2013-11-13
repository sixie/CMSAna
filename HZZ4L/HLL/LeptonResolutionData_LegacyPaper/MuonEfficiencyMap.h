UInt_t FindMuonEfficiencyBin( double value, double bins[], UInt_t nbins) {
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


Double_t GetMuonEfficiencyPtEtaPhi(Double_t Pt, Double_t Eta, Double_t Phi) {

  Double_t ptBins[16] = {5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,50};
  Double_t etaBins[17] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4442,1.566,1.8,2,2.1,2.2,2.3,2.4,2.5,2.6};
  Double_t phiBins[13] = {-3.2,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3.2};


  Double_t Efficiency[17][18][14]  = {
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.0041833,0.00169164,0.00219941,0.00323625,0.0016568,0.00340053,0.00246002,0.00284225,0.0026712,0.00189934,0.00243962,0.00323071,0},
      {0,0.00305227,0.0029725,0.00219191,0.004003,0.00319018,0.00219191,0.00343643,0.0043881,0.00165916,0.00390911,0.00574282,0.00286096,0},
      {0,0.00364683,0.00474644,0.00219405,0.00579491,0.00365943,0.00475475,0.0035515,0.0035824,0.00383436,0.0039555,0.00395159,0.00291149,0},
      {0,0.0034034,0.00446077,0.00309039,0.00608365,0.00417101,0.00384813,0.00324919,0.00233584,0.00412903,0.00431801,0.00453629,0.00442211,0},
      {0,0.00685643,0.0057128,0.00475059,0.00531208,0.00647075,0.00555409,0.00561347,0.00397878,0.0046133,0.00392876,0.00560598,0.00492914,0},
      {0,0.00689195,0.00370054,0.00655831,0.00612643,0.00601547,0.0047196,0.0044906,0.0044906,0.00511073,0.00453515,0.00514433,0.00663717,0},
      {0,0.00717116,0.006917,0.00434678,0.00508598,0.00749625,0.00496032,0.00544379,0.00721214,0.00469252,0.00559067,0.00636475,0.00707051,0},
      {0,0.00616776,0.00459184,0.00363636,0.00689655,0.00313152,0.00716846,0.00538503,0.00655738,0.00363259,0.00382305,0.00401003,0.00499584,0},
      {0,0.00448113,0.00440917,0.00512975,0.00398284,0.00702165,0.00443131,0.0060024,0.00533333,0.00567334,0.00449102,0.00658092,0.00516796,0},
      {0,0.00452216,0.0071885,0.00537222,0.00763942,0.00620636,0.00628931,0.00550964,0.00626714,0.00461007,0.00784314,0.00540123,0.00461113,0},
      {0,0.00392157,0.00958466,0.00594732,0.0100587,0.00642792,0.00428816,0.0085034,0.00740132,0.0062167,0.0066778,0.0111972,0.00480769,0},
      {0,0.0034965,0.00453309,0.00777454,0.00633484,0.0117223,0.0100457,0.000910747,0.00973451,0.00353045,0.0107914,0.00184502,0.00839748,0},
      {0,0.00846805,0.00605449,0.00585366,0.00867052,0.00622407,0.0127326,0.00691017,0.0056872,0.0077821,0.00690335,0.00673725,0.00306279,0},
      {0,0.00738916,0.00806452,0.0040404,0.00919305,0.00652174,0.00827301,0.00408998,0.0042508,0.004995,0.00616016,0.00523013,0.00661704,0},
      {0,0.000896861,0,0.00110988,0.00118203,0,0.00114548,0,0,0.00113766,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.692804,0.67823,0.672297,0.706832,0.701923,0.669056,0.678788,0.693154,0.681373,0.694994,0.679952,0.699334,0},
      {0,0.697697,0.696576,0.691748,0.694836,0.679245,0.688586,0.726221,0.685882,0.703518,0.684729,0.698856,0.698351,0},
      {0,0.710425,0.707809,0.689655,0.732573,0.696765,0.722222,0.712362,0.723636,0.717224,0.733978,0.724782,0.706511,0},
      {0,0.715757,0.722295,0.708877,0.749064,0.73868,0.75443,0.725845,0.740741,0.701351,0.72346,0.745198,0.736354,0},
      {0,0.725449,0.727038,0.730769,0.70227,0.715847,0.727891,0.720506,0.745536,0.727156,0.738571,0.745257,0.751351,0},
      {0,0.732162,0.72863,0.694521,0.740841,0.702128,0.725705,0.744741,0.715288,0.707493,0.724234,0.742657,0.758511,0},
      {0,0.735736,0.735936,0.753713,0.715385,0.752741,0.764554,0.726918,0.735149,0.762953,0.735409,0.734127,0.737512,0},
      {0,0.791845,0.770235,0.73224,0.731383,0.762058,0.750725,0.712846,0.753316,0.712766,0.761658,0.780702,0.697826,0},
      {0,0.762452,0.781056,0.757377,0.760976,0.804878,0.767267,0.771518,0.768519,0.748115,0.760976,0.778434,0.768029,0},
      {0,0.803863,0.746094,0.779443,0.74375,0.762846,0.770961,0.779048,0.770115,0.788538,0.775391,0.794118,0.7648,0},
      {0,0.735915,0.71308,0.746606,0.714286,0.752212,0.716495,0.725118,0.797235,0.763158,0.715026,0.741784,0.738462,0},
      {0,0.670251,0.729885,0.649485,0.716763,0.746411,0.640212,0.708738,0.715789,0.646226,0.758242,0.697115,0.706349,0},
      {0,0.65,0.606383,0.604061,0.560209,0.630542,0.656977,0.661376,0.534884,0.677249,0.636364,0.611111,0.651064,0},
      {0,0.473029,0.439024,0.461538,0.484076,0.502994,0.576923,0.489655,0.55625,0.5,0.564972,0.516854,0.486486,0},
      {0,0.0049505,0,0,0,0.00609756,0,0,0,0,0.0125,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.773526,0.74937,0.737637,0.767352,0.777778,0.735294,0.746851,0.749296,0.780856,0.788535,0.748705,0.761905,0},
      {0,0.778656,0.772368,0.765092,0.754248,0.811295,0.761905,0.761598,0.779661,0.754098,0.776099,0.797011,0.78821,0},
      {0,0.778243,0.792135,0.764946,0.784553,0.792857,0.7867,0.780612,0.777328,0.765799,0.784211,0.798621,0.782435,0},
      {0,0.784292,0.785319,0.820793,0.822626,0.797386,0.819783,0.81456,0.812256,0.787599,0.776159,0.781534,0.810811,0},
      {0,0.808036,0.804878,0.77003,0.805746,0.798551,0.772263,0.785607,0.809192,0.838757,0.804809,0.798592,0.801624,0},
      {0,0.817527,0.826953,0.787736,0.816641,0.805901,0.832386,0.839813,0.838168,0.828528,0.828571,0.797806,0.814516,0},
      {0,0.841206,0.824828,0.806319,0.831625,0.813472,0.824147,0.82133,0.832418,0.83355,0.828146,0.805369,0.819979,0},
      {0,0.841772,0.836676,0.824561,0.841499,0.863372,0.838323,0.814493,0.824513,0.837912,0.837349,0.839763,0.837037,0},
      {0,0.833744,0.840426,0.85173,0.856419,0.8304,0.815152,0.825796,0.847377,0.866889,0.846281,0.843023,0.835634,0},
      {0,0.839939,0.807792,0.82167,0.809091,0.87473,0.853712,0.845433,0.858696,0.852679,0.851582,0.843891,0.844365,0},
      {0,0.819923,0.829493,0.771084,0.777251,0.811594,0.884058,0.79703,0.789216,0.843318,0.872449,0.857759,0.838951,0},
      {0,0.788136,0.764045,0.818182,0.738095,0.79096,0.769231,0.752941,0.789744,0.75,0.814815,0.795699,0.787879,0},
      {0,0.738938,0.664835,0.719298,0.698324,0.69375,0.736527,0.714286,0.688623,0.788043,0.756757,0.699454,0.726368,0},
      {0,0.592417,0.524096,0.642424,0.624204,0.62,0.691176,0.568047,0.636943,0.636364,0.61039,0.622754,0.595745,0},
      {0,0.00534759,0.00833333,0,0,0,0,0.0142857,0,0,0,0,0.00581395,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.836066,0.801471,0.80083,0.796424,0.781122,0.815612,0.806638,0.796919,0.802557,0.794366,0.818057,0.816435,0},
      {0,0.822839,0.808424,0.819885,0.801144,0.821683,0.811533,0.79805,0.813953,0.795362,0.816901,0.821162,0.804878,0},
      {0,0.841866,0.811852,0.789404,0.822314,0.821522,0.835374,0.817923,0.83115,0.818841,0.793696,0.814815,0.810268,0},
      {0,0.831498,0.851003,0.850227,0.847076,0.854467,0.868182,0.815668,0.862573,0.836283,0.822767,0.85782,0.832134,0},
      {0,0.8478,0.843195,0.796267,0.848576,0.836767,0.840121,0.827922,0.868502,0.843849,0.833608,0.822294,0.872749,0},
      {0,0.849367,0.833061,0.838957,0.849057,0.846761,0.835355,0.801303,0.83982,0.865574,0.853698,0.835526,0.883333,0},
      {0,0.856981,0.829023,0.867521,0.847305,0.841667,0.867069,0.844729,0.88169,0.868239,0.83795,0.878561,0.877483,0},
      {0,0.867299,0.87619,0.848387,0.858896,0.879888,0.870091,0.832237,0.879257,0.863777,0.875723,0.858025,0.884804,0},
      {0,0.856164,0.88676,0.842199,0.81993,0.862385,0.882459,0.831283,0.861301,0.868571,0.85689,0.865922,0.864407,0},
      {0,0.858657,0.875,0.863109,0.882353,0.878109,0.87239,0.880841,0.857482,0.862981,0.863753,0.851415,0.865248,0},
      {0,0.865079,0.815,0.826816,0.817204,0.857868,0.882682,0.835106,0.879581,0.814815,0.827027,0.863905,0.845528,0},
      {0,0.765217,0.795455,0.833333,0.776042,0.779874,0.805263,0.804734,0.8125,0.863354,0.820988,0.772222,0.834711,0},
      {0,0.778947,0.805714,0.719512,0.706587,0.708333,0.781421,0.748428,0.722973,0.8375,0.835366,0.777778,0.752525,0},
      {0,0.633721,0.666667,0.622951,0.613139,0.590909,0.597484,0.70229,0.645161,0.710843,0.689873,0.695364,0.606061,0},
      {0,0.00558659,0,0.013986,0.00787402,0.00763359,0.016,0.00775194,0,0,0,0.0172414,0.00555556,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.830108,0.811953,0.812336,0.820961,0.837762,0.828986,0.859748,0.853247,0.832215,0.840941,0.829371,0.829621,0},
      {0,0.848786,0.825269,0.848821,0.826277,0.854676,0.817391,0.822129,0.851182,0.866571,0.844286,0.856354,0.858521,0},
      {0,0.84307,0.853625,0.830585,0.867486,0.847737,0.874063,0.832196,0.848787,0.856932,0.839827,0.869873,0.836639,0},
      {0,0.849633,0.875179,0.850679,0.859813,0.847262,0.86119,0.882521,0.858209,0.861516,0.841958,0.830882,0.87067,0},
      {0,0.875949,0.843271,0.852751,0.885759,0.880259,0.854071,0.897227,0.877061,0.861194,0.881064,0.880551,0.855006,0},
      {0,0.887793,0.881739,0.857143,0.852041,0.881739,0.861066,0.867153,0.860504,0.877193,0.871753,0.883257,0.879548,0},
      {0,0.874272,0.860031,0.841705,0.872894,0.872024,0.885075,0.879687,0.899701,0.869499,0.864266,0.851744,0.865654,0},
      {0,0.84399,0.884273,0.891089,0.848665,0.913313,0.886792,0.877358,0.86478,0.890323,0.863924,0.892256,0.890819,0},
      {0,0.899593,0.863378,0.888689,0.876307,0.886824,0.902111,0.868613,0.86116,0.883142,0.862069,0.913121,0.894737,0},
      {0,0.879377,0.90709,0.9022,0.872123,0.890306,0.889401,0.888626,0.909794,0.880899,0.888361,0.87381,0.876712,0},
      {0,0.839552,0.860963,0.863158,0.809524,0.898148,0.892473,0.844444,0.834225,0.845455,0.831579,0.8867,0.873518,0},
      {0,0.844444,0.8,0.81592,0.881988,0.796296,0.835227,0.849673,0.819209,0.886364,0.859756,0.855263,0.819742,0},
      {0,0.766169,0.729167,0.77439,0.739437,0.791367,0.802326,0.75,0.80226,0.80625,0.8,0.753333,0.822581,0},
      {0,0.65,0.670886,0.695946,0.60625,0.604651,0.646259,0.642857,0.6875,0.742857,0.657718,0.741667,0.65896,0},
      {0,0,0,0.00724638,0,0,0,0,0.00877193,0,0.00869565,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.834731,0.842814,0.848739,0.84816,0.860741,0.818182,0.872464,0.885794,0.851689,0.861919,0.857351,0.845815,0},
      {0,0.87037,0.864431,0.811782,0.846377,0.874608,0.853659,0.841808,0.841615,0.861968,0.864706,0.88006,0.854606,0},
      {0,0.876586,0.866473,0.880368,0.86,0.868463,0.862661,0.884615,0.853411,0.882448,0.882526,0.875367,0.879433,0},
      {0,0.87696,0.855974,0.862595,0.868878,0.880911,0.872621,0.836041,0.864865,0.884679,0.86385,0.869775,0.86231,0},
      {0,0.896071,0.905344,0.861027,0.86901,0.855462,0.875981,0.891339,0.854063,0.885329,0.897727,0.8848,0.872801,0},
      {0,0.869727,0.890017,0.875828,0.899497,0.898223,0.890485,0.889081,0.871287,0.898305,0.870337,0.905336,0.901547,0},
      {0,0.897377,0.905145,0.882979,0.868976,0.887755,0.880184,0.866562,0.894984,0.871355,0.895522,0.889065,0.901408,0},
      {0,0.899263,0.911538,0.896875,0.901316,0.92163,0.912338,0.891892,0.895683,0.871711,0.845614,0.901587,0.879717,0},
      {0,0.891854,0.897839,0.9,0.886691,0.897921,0.908463,0.883534,0.886239,0.895753,0.902672,0.905838,0.89955,0},
      {0,0.909091,0.886842,0.887019,0.906634,0.892768,0.905612,0.868812,0.901554,0.892694,0.9,0.887399,0.913556,0},
      {0,0.902128,0.879121,0.858696,0.871134,0.9,0.91018,0.853933,0.870466,0.874286,0.892473,0.911458,0.915094,0},
      {0,0.869565,0.838323,0.843243,0.850932,0.829268,0.850932,0.866667,0.90184,0.891566,0.883041,0.817073,0.875598,0},
      {0,0.785714,0.786667,0.810219,0.830601,0.791367,0.829268,0.873239,0.839744,0.834532,0.802395,0.80663,0.773196,0},
      {0,0.744318,0.671429,0.679389,0.681818,0.73913,0.75,0.766423,0.771429,0.694915,0.803797,0.75969,0.686391,0},
      {0,0,0.00813008,0.00877193,0,0.00925926,0,0,0,0,0.00862069,0,0.0136986,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.885289,0.881874,0.867847,0.882883,0.901038,0.894958,0.891137,0.887973,0.886396,0.892401,0.882546,0.884533,0},
      {0,0.892485,0.887973,0.881974,0.88132,0.887295,0.897122,0.890738,0.892659,0.884501,0.888657,0.88931,0.89557,0},
      {0,0.891464,0.900808,0.889288,0.870448,0.892679,0.874145,0.888356,0.89338,0.895478,0.891048,0.878592,0.891456,0},
      {0,0.9134,0.904048,0.90717,0.910663,0.905415,0.894509,0.897644,0.886888,0.908541,0.911411,0.903226,0.889632,0},
      {0,0.900117,0.910162,0.910104,0.909366,0.913012,0.902821,0.895652,0.894813,0.900234,0.885565,0.90462,0.917586,0},
      {0,0.912552,0.905321,0.900314,0.910053,0.894108,0.906107,0.924663,0.911111,0.925175,0.898387,0.906793,0.905549,0},
      {0,0.913323,0.905144,0.89508,0.903034,0.90411,0.908616,0.898199,0.898438,0.922603,0.906317,0.927322,0.906617,0},
      {0,0.909091,0.915629,0.893466,0.928264,0.926205,0.897038,0.914286,0.920659,0.901709,0.927469,0.912463,0.928655,0},
      {0,0.908856,0.925681,0.909657,0.933836,0.912309,0.915584,0.906504,0.913478,0.921939,0.927101,0.92635,0.922236,0},
      {0,0.918941,0.930886,0.921892,0.916392,0.912556,0.918182,0.909774,0.931697,0.916149,0.913811,0.916314,0.919118,0},
      {0,0.919065,0.913043,0.911111,0.913043,0.922018,0.923788,0.914798,0.92,0.899533,0.89485,0.919903,0.935368,0},
      {0,0.897839,0.909774,0.9125,0.858247,0.871795,0.862981,0.887139,0.92654,0.884236,0.890274,0.895,0.905482,0},
      {0,0.794926,0.817949,0.832041,0.816216,0.805755,0.847594,0.881137,0.870523,0.848649,0.845178,0.846966,0.859756,0},
      {0,0.753555,0.755495,0.77135,0.751678,0.75969,0.770718,0.781818,0.709581,0.781421,0.784062,0.782723,0.770531,0},
      {0,0.015748,0.00471698,0,0.0135747,0,0,0.005,0.00534759,0.011976,0.0046729,0.00921659,0.00393701,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.905131,0.901799,0.899925,0.901138,0.894537,0.899644,0.885221,0.913562,0.90176,0.888975,0.902493,0.899886,0},
      {0,0.916282,0.913176,0.902788,0.901311,0.911786,0.902883,0.914243,0.904,0.913043,0.909091,0.911179,0.914729,0},
      {0,0.917015,0.913995,0.894537,0.915985,0.914115,0.904586,0.904726,0.912481,0.914877,0.911442,0.917882,0.903432,0},
      {0,0.936347,0.921175,0.91245,0.9297,0.917708,0.918524,0.89572,0.920712,0.926667,0.937453,0.92459,0.918626,0},
      {0,0.927781,0.913149,0.920034,0.919184,0.9181,0.911483,0.919032,0.925627,0.916272,0.916933,0.919969,0.917877,0},
      {0,0.931803,0.924179,0.920775,0.927512,0.919932,0.90996,0.917181,0.928999,0.933891,0.908155,0.920465,0.925852,0},
      {0,0.930023,0.917638,0.913702,0.9274,0.924981,0.919943,0.920973,0.937408,0.9041,0.934181,0.920623,0.918854,0},
      {0,0.919949,0.926791,0.933649,0.943615,0.925806,0.929134,0.927803,0.941363,0.934621,0.935229,0.934718,0.930292,0},
      {0,0.951236,0.938959,0.938865,0.939421,0.936481,0.928634,0.924702,0.928637,0.929266,0.933102,0.925926,0.930055,0},
      {0,0.93231,0.934905,0.934633,0.921546,0.920725,0.934884,0.927252,0.938845,0.947239,0.937729,0.947778,0.936522,0},
      {0,0.933194,0.899244,0.9075,0.896552,0.914286,0.931759,0.914286,0.917476,0.928747,0.953125,0.937028,0.939148,0},
      {0,0.902196,0.917112,0.910941,0.912371,0.915789,0.888268,0.90113,0.91129,0.911688,0.9175,0.91092,0.925553,0},
      {0,0.857143,0.848214,0.874286,0.885117,0.87055,0.880435,0.875371,0.869159,0.880519,0.886305,0.863222,0.87257,0},
      {0,0.801471,0.761755,0.81388,0.788079,0.795527,0.83642,0.765079,0.819767,0.860534,0.828746,0.766272,0.792941,0},
      {0,0.00796813,0.00518135,0.00653595,0.0119048,0.00966184,0,0.0175439,0.0126582,0.0052356,0,0.0111732,0.0168776,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.923839,0.93297,0.927877,0.924267,0.934528,0.920823,0.913738,0.92428,0.904577,0.920572,0.929688,0.932639,0},
      {0,0.9325,0.924937,0.904485,0.910941,0.916944,0.915117,0.93066,0.912458,0.917923,0.919275,0.922638,0.937657,0},
      {0,0.927733,0.93263,0.916399,0.930364,0.930788,0.91885,0.915493,0.923204,0.945217,0.936239,0.940887,0.922372,0},
      {0,0.931701,0.939241,0.938655,0.936833,0.932287,0.928159,0.934095,0.916807,0.934783,0.916739,0.932594,0.932712,0},
      {0,0.918768,0.918573,0.917749,0.916353,0.933954,0.925957,0.918584,0.931878,0.92767,0.933099,0.919311,0.934754,0},
      {0,0.936893,0.935614,0.92,0.921409,0.9375,0.936626,0.928372,0.92885,0.937557,0.93007,0.93231,0.930348,0},
      {0,0.943653,0.931692,0.937746,0.933871,0.937903,0.939469,0.931596,0.939123,0.921584,0.93029,0.92823,0.932915,0},
      {0,0.934783,0.940375,0.938038,0.957983,0.942078,0.94086,0.943428,0.937086,0.944171,0.928333,0.926789,0.924264,0},
      {0,0.948276,0.954163,0.948718,0.938462,0.943044,0.934132,0.917166,0.924386,0.944551,0.938192,0.939061,0.935557,0},
      {0,0.940312,0.94898,0.942602,0.954887,0.922208,0.922597,0.955204,0.951673,0.959647,0.963183,0.959792,0.932814,0},
      {0,0.934156,0.940104,0.931122,0.934959,0.930412,0.925208,0.941011,0.936813,0.925558,0.945455,0.949495,0.929006,0},
      {0,0.93932,0.919643,0.919643,0.896875,0.933962,0.914454,0.896848,0.902913,0.923913,0.933333,0.952113,0.93617,0},
      {0,0.898618,0.862295,0.89899,0.877358,0.901899,0.889571,0.890323,0.904762,0.912121,0.901587,0.900958,0.898345,0},
      {0,0.842377,0.847458,0.817276,0.85,0.801444,0.875421,0.841912,0.825215,0.846645,0.871186,0.851852,0.861878,0},
      {0,0.0100503,0.0127389,0,0.00588235,0.00746269,0.00571429,0.0171429,0.00555556,0.00662252,0.00653595,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.928709,0.931677,0.916307,0.934946,0.937713,0.925311,0.915212,0.928139,0.936546,0.940747,0.92916,0.934405,0},
      {0,0.928667,0.934002,0.9408,0.946281,0.922111,0.938038,0.924481,0.926494,0.921348,0.940972,0.932331,0.92743,0},
      {0,0.932718,0.940135,0.930556,0.932088,0.935648,0.936688,0.926307,0.93617,0.944043,0.934426,0.937288,0.948113,0},
      {0,0.93456,0.943925,0.933102,0.937133,0.938811,0.940273,0.945346,0.936073,0.933447,0.926981,0.934932,0.937457,0},
      {0,0.933866,0.949087,0.939123,0.937095,0.934105,0.928758,0.944231,0.935077,0.942959,0.938326,0.932829,0.947368,0},
      {0,0.938776,0.93432,0.923579,0.95015,0.936087,0.939703,0.942884,0.925263,0.942507,0.933075,0.950249,0.950856,0},
      {0,0.950502,0.945968,0.92721,0.937604,0.927001,0.949033,0.940348,0.943316,0.939035,0.944444,0.949367,0.939009,0},
      {0,0.939355,0.952118,0.953043,0.932971,0.953358,0.936731,0.942857,0.952984,0.930018,0.945946,0.945255,0.943089,0},
      {0,0.940554,0.947585,0.946626,0.934985,0.948913,0.951659,0.937113,0.924037,0.956389,0.960847,0.964611,0.948638,0},
      {0,0.939555,0.950798,0.949807,0.959514,0.945431,0.947503,0.945652,0.941748,0.957474,0.941799,0.946615,0.960163,0},
      {0,0.946067,0.92545,0.958824,0.913747,0.934426,0.94051,0.936963,0.946176,0.951429,0.951952,0.962963,0.964059,0},
      {0,0.929504,0.930303,0.963855,0.905172,0.908571,0.949405,0.940439,0.952113,0.934328,0.939683,0.950867,0.920398,0},
      {0,0.918159,0.896774,0.885196,0.917722,0.932515,0.89547,0.909091,0.916418,0.909357,0.902017,0.915152,0.903073,0},
      {0,0.844985,0.855422,0.837545,0.824818,0.820789,0.859712,0.821705,0.857678,0.922509,0.878229,0.837302,0.865889,0},
      {0,0.00952381,0.00591716,0,0,0,0.0212766,0,0,0.00689655,0.00653595,0.00684932,0.0251256,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.941447,0.937008,0.929276,0.931204,0.938272,0.932701,0.939597,0.939108,0.932546,0.943709,0.937807,0.94075,0},
      {0,0.949208,0.945874,0.94103,0.936644,0.950085,0.943994,0.937965,0.935426,0.935965,0.941176,0.943445,0.936725,0},
      {0,0.950306,0.936393,0.950545,0.933002,0.942342,0.936306,0.932465,0.948582,0.9375,0.942758,0.947707,0.949935,0},
      {0,0.952006,0.94411,0.939014,0.939736,0.934876,0.944395,0.939231,0.943447,0.942105,0.949731,0.932337,0.934631,0},
      {0,0.948699,0.950185,0.948004,0.934545,0.93934,0.940789,0.930683,0.941942,0.939221,0.945005,0.949541,0.936264,0},
      {0,0.948407,0.940576,0.938776,0.947516,0.951265,0.947419,0.940442,0.938017,0.943284,0.945328,0.949393,0.942192,0},
      {0,0.934856,0.931034,0.954426,0.949097,0.948427,0.945378,0.947644,0.950607,0.939292,0.955274,0.942831,0.955182,0},
      {0,0.962293,0.953704,0.954128,0.926136,0.963462,0.952381,0.947896,0.951036,0.946524,0.942639,0.954635,0.94835,0},
      {0,0.951954,0.950166,0.967775,0.957001,0.967249,0.950604,0.938559,0.949772,0.946114,0.953511,0.947312,0.95436,0},
      {0,0.940362,0.95,0.953782,0.960437,0.949934,0.944595,0.962008,0.949176,0.951157,0.955357,0.952,0.962012,0},
      {0,0.954106,0.948795,0.930168,0.938028,0.947214,0.930556,0.961111,0.959064,0.966387,0.962366,0.964179,0.950104,0},
      {0,0.930233,0.92381,0.919444,0.949254,0.923077,0.946588,0.913333,0.916667,0.93617,0.942693,0.957377,0.960591,0},
      {0,0.927461,0.921824,0.948529,0.906516,0.898246,0.930091,0.938776,0.913333,0.932692,0.925651,0.925466,0.921519,0},
      {0,0.900302,0.868132,0.888889,0.862816,0.90681,0.88716,0.839695,0.858779,0.87108,0.881188,0.878676,0.900277,0},
      {0,0.0168539,0.013986,0,0.0150376,0,0,0.00719424,0.0075188,0,0.00675676,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.952291,0.953026,0.948429,0.948869,0.949231,0.95744,0.94364,0.951087,0.949075,0.949985,0.95122,0.948425,0},
      {0,0.952767,0.952116,0.953061,0.953155,0.949021,0.952038,0.947079,0.947433,0.950352,0.95403,0.945831,0.95569,0},
      {0,0.954746,0.947815,0.956347,0.951969,0.954445,0.958513,0.949937,0.955022,0.952995,0.95148,0.952605,0.95449,0},
      {0,0.955711,0.966209,0.958743,0.959783,0.951899,0.949735,0.956228,0.954368,0.958414,0.955266,0.954964,0.953311,0},
      {0,0.956032,0.948649,0.955873,0.955877,0.950086,0.94993,0.956671,0.961706,0.950786,0.958291,0.954148,0.957509,0},
      {0,0.958634,0.955862,0.953902,0.964693,0.952312,0.955604,0.949257,0.955068,0.955625,0.9525,0.953471,0.953444,0},
      {0,0.960543,0.954949,0.948329,0.950158,0.960655,0.95636,0.958753,0.954915,0.955651,0.96124,0.952786,0.95869,0},
      {0,0.957961,0.96477,0.958859,0.958253,0.967363,0.954274,0.963365,0.962636,0.961115,0.966473,0.966837,0.963596,0},
      {0,0.964575,0.957236,0.964973,0.963719,0.961713,0.962064,0.950122,0.954014,0.957584,0.963057,0.956422,0.963204,0},
      {0,0.958828,0.970302,0.971145,0.96117,0.966045,0.964861,0.967679,0.96063,0.964706,0.970273,0.96854,0.968402,0},
      {0,0.959968,0.94917,0.960499,0.942768,0.957627,0.962475,0.957778,0.972458,0.959826,0.964356,0.957292,0.97551,0},
      {0,0.943089,0.938664,0.961332,0.94214,0.949721,0.955044,0.94311,0.944191,0.948124,0.950474,0.955882,0.956835,0},
      {0,0.920737,0.934195,0.958525,0.927145,0.953434,0.936905,0.927928,0.936959,0.950363,0.941514,0.942377,0.945455,0},
      {0,0.895075,0.891746,0.915646,0.901408,0.920912,0.911323,0.897474,0.898148,0.926923,0.931559,0.889043,0.920668,0},
      {0,0.0125313,0,0.0141844,0.00320513,0.0131148,0.00983607,0.00884956,0.00724638,0.00955414,0.0118343,0.00307692,0.00797872,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.965887,0.960875,0.958934,0.954475,0.962769,0.953037,0.954487,0.954404,0.960977,0.956588,0.952396,0.957529,0},
      {0,0.960943,0.95804,0.965245,0.957013,0.954897,0.955983,0.957938,0.955799,0.95296,0.956183,0.956508,0.959596,0},
      {0,0.959817,0.958587,0.964984,0.968023,0.95717,0.953686,0.95921,0.959725,0.957663,0.963535,0.958428,0.960376,0},
      {0,0.96027,0.961473,0.963523,0.958447,0.965778,0.967322,0.960169,0.96142,0.964015,0.957475,0.955224,0.961226,0},
      {0,0.966713,0.963826,0.959717,0.956614,0.965917,0.964596,0.960187,0.95785,0.959677,0.962624,0.959382,0.969173,0},
      {0,0.961211,0.966654,0.964273,0.960212,0.96365,0.959581,0.962257,0.962039,0.964299,0.959262,0.956124,0.964923,0},
      {0,0.964156,0.954855,0.95416,0.953061,0.964167,0.960081,0.964129,0.968658,0.956479,0.961178,0.960733,0.961832,0},
      {0,0.968386,0.963519,0.954386,0.961783,0.964489,0.963925,0.961644,0.961871,0.963493,0.968622,0.963889,0.962863,0},
      {0,0.964907,0.967301,0.969591,0.964842,0.96354,0.961295,0.96042,0.958904,0.962104,0.970451,0.965006,0.972205,0},
      {0,0.961391,0.967146,0.96611,0.966032,0.973926,0.967271,0.96628,0.972344,0.967181,0.973863,0.961929,0.964444,0},
      {0,0.969724,0.971061,0.97,0.959137,0.966316,0.973539,0.966427,0.966555,0.961024,0.97377,0.952874,0.963466,0},
      {0,0.955268,0.961313,0.965027,0.967941,0.949878,0.96368,0.953865,0.967459,0.970964,0.95321,0.951878,0.953405,0},
      {0,0.942058,0.95092,0.954436,0.948622,0.946809,0.950303,0.93734,0.955451,0.948212,0.960957,0.950556,0.956478,0},
      {0,0.909404,0.943262,0.921379,0.941748,0.913917,0.940256,0.935681,0.941718,0.93918,0.927864,0.941512,0.933619,0},
      {0,0.013624,0.0299625,0.00378788,0.013986,0.0104895,0.018315,0.0200669,0.0133779,0,0.00357143,0.00701754,0.0123077,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.964999,0.967514,0.961408,0.966985,0.957624,0.9625,0.959098,0.964945,0.959962,0.963146,0.964754,0.960982,0},
      {0,0.968159,0.964585,0.959519,0.95851,0.96532,0.96919,0.96554,0.962319,0.957704,0.964356,0.967782,0.961567,0},
      {0,0.968639,0.970846,0.962822,0.965177,0.969279,0.971438,0.96554,0.96126,0.963489,0.963183,0.97032,0.966982,0},
      {0,0.961549,0.97185,0.964121,0.966779,0.968354,0.964762,0.965835,0.969215,0.970338,0.97009,0.966135,0.967369,0},
      {0,0.969402,0.965062,0.968935,0.967434,0.96088,0.968017,0.967025,0.967951,0.969046,0.964184,0.968451,0.971947,0},
      {0,0.969226,0.969136,0.968177,0.970313,0.971883,0.961595,0.966568,0.964193,0.966704,0.97158,0.968808,0.964413,0},
      {0,0.966957,0.961398,0.956975,0.958361,0.970179,0.962554,0.967066,0.969536,0.970244,0.96658,0.964653,0.968703,0},
      {0,0.969835,0.965768,0.977512,0.971935,0.970051,0.955619,0.963483,0.967972,0.964481,0.971705,0.975212,0.966377,0},
      {0,0.969378,0.974949,0.970576,0.979253,0.969607,0.971336,0.961231,0.964214,0.969833,0.971625,0.970185,0.965464,0},
      {0,0.970454,0.969239,0.975367,0.971029,0.969964,0.973456,0.966249,0.966171,0.970483,0.971007,0.978437,0.970445,0},
      {0,0.968577,0.975664,0.972789,0.967356,0.962195,0.968198,0.970623,0.978648,0.980222,0.973804,0.971864,0.968326,0},
      {0,0.960486,0.966357,0.965473,0.962092,0.967941,0.959427,0.964912,0.962042,0.973138,0.973964,0.97598,0.960674,0},
      {0,0.940412,0.969974,0.953401,0.94709,0.955243,0.966501,0.951524,0.967785,0.965612,0.95478,0.939142,0.969419,0},
      {0,0.923451,0.944954,0.944993,0.945827,0.937097,0.946853,0.92623,0.948966,0.954745,0.938659,0.943114,0.947778,0},
      {0,0.0183486,0.0150943,0.0127119,0.00813008,0.0107527,0.00408163,0.0110294,0.0181159,0.0229008,0.0115385,0.021645,0.00961538,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.974671,0.967316,0.966333,0.971878,0.974787,0.971268,0.964946,0.966409,0.972297,0.967646,0.969434,0.963078,0},
      {0,0.971568,0.969628,0.970874,0.971245,0.971466,0.971391,0.973303,0.969069,0.972864,0.973081,0.970696,0.971098,0},
      {0,0.971624,0.969409,0.968403,0.969307,0.972909,0.963704,0.965505,0.975321,0.969458,0.970447,0.971593,0.970677,0},
      {0,0.967604,0.973393,0.961687,0.970262,0.971746,0.972696,0.965517,0.977234,0.96203,0.97079,0.972909,0.977224,0},
      {0,0.972262,0.975467,0.97201,0.964529,0.973731,0.969652,0.966667,0.969925,0.96787,0.973012,0.969842,0.971831,0},
      {0,0.969305,0.965918,0.970067,0.963636,0.967359,0.962264,0.968474,0.966043,0.972582,0.974865,0.970462,0.970685,0},
      {0,0.972289,0.966014,0.964706,0.966177,0.964578,0.970054,0.969107,0.966136,0.968168,0.974701,0.969533,0.970113,0},
      {0,0.976852,0.973786,0.978478,0.962052,0.981871,0.970939,0.967189,0.959298,0.972561,0.976848,0.975753,0.968595,0},
      {0,0.973702,0.975453,0.977236,0.973362,0.977556,0.96867,0.971014,0.962629,0.970662,0.973229,0.972138,0.973919,0},
      {0,0.975254,0.974632,0.974947,0.972785,0.975877,0.972441,0.969898,0.98172,0.972771,0.974526,0.973813,0.974196,0},
      {0,0.96738,0.971461,0.972696,0.959158,0.971178,0.975933,0.973398,0.96988,0.971631,0.971257,0.977053,0.972469,0},
      {0,0.965921,0.969169,0.977424,0.957254,0.958904,0.97037,0.974497,0.972261,0.976131,0.959658,0.960265,0.969,0},
      {0,0.926012,0.949333,0.962687,0.966245,0.956585,0.959128,0.946071,0.954733,0.976712,0.961433,0.964286,0.970522,0},
      {0,0.920755,0.945161,0.957746,0.946372,0.946367,0.955556,0.929766,0.95624,0.959459,0.951104,0.962963,0.955975,0},
      {0,0.0140351,0.0141509,0.00454545,0.0134529,0,0.0166667,0.00961538,0.012605,0.0078125,0,0.00431034,0.00704225,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.96768,0.972488,0.97107,0.964972,0.975951,0.97079,0.969448,0.974059,0.973197,0.975552,0.97037,0.973041,0},
      {0,0.970432,0.972979,0.976157,0.968368,0.974813,0.970121,0.972507,0.96737,0.966298,0.971138,0.976119,0.96858,0},
      {0,0.972318,0.976081,0.975873,0.971784,0.975523,0.973472,0.975137,0.972106,0.970037,0.974063,0.973533,0.973993,0},
      {0,0.974986,0.976551,0.97495,0.970546,0.972229,0.971018,0.97245,0.975701,0.974957,0.969644,0.972103,0.974277,0},
      {0,0.971143,0.973244,0.968482,0.972271,0.978102,0.968096,0.972049,0.973558,0.971802,0.968821,0.97468,0.973207,0},
      {0,0.971303,0.972563,0.970125,0.972798,0.974574,0.967709,0.974257,0.97291,0.971831,0.972603,0.973664,0.976429,0},
      {0,0.97026,0.96835,0.963356,0.96689,0.969122,0.968876,0.969854,0.971709,0.968669,0.972311,0.976127,0.966348,0},
      {0,0.970843,0.973575,0.972772,0.969652,0.975223,0.97398,0.972596,0.972892,0.977471,0.970326,0.975873,0.974704,0},
      {0,0.974657,0.970622,0.973684,0.97327,0.97975,0.972032,0.968664,0.97142,0.971396,0.973609,0.976453,0.97692,0},
      {0,0.979641,0.979159,0.973694,0.974002,0.973291,0.975313,0.976569,0.97702,0.979442,0.973378,0.979465,0.975797,0},
      {0,0.971641,0.957322,0.974957,0.970991,0.970563,0.971843,0.969512,0.976151,0.968929,0.974316,0.97515,0.980959,0},
      {0,0.970739,0.966241,0.976895,0.961326,0.971197,0.970534,0.95942,0.973321,0.962142,0.976534,0.982309,0.966381,0},
      {0,0.954798,0.970356,0.959204,0.967285,0.967513,0.974464,0.955466,0.974012,0.969548,0.978988,0.967342,0.96313,0},
      {0,0.940171,0.94898,0.955504,0.950243,0.95881,0.956977,0.940024,0.952719,0.955422,0.946619,0.959491,0.969557,0},
      {0,0.00665188,0.0031348,0.00617284,0.00298507,0.0112994,0.0141243,0.0027933,0.015625,0.00287356,0.00630915,0.012945,0.0133929,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    },
    {
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0.974885,0.971802,0.974969,0.970474,0.974129,0.974528,0.971452,0.971834,0.972393,0.974328,0.971953,0.96995,0},
      {0,0.975459,0.977315,0.97163,0.973722,0.976492,0.975889,0.974634,0.968578,0.970308,0.974898,0.973755,0.978576,0},
      {0,0.972947,0.977496,0.977982,0.979369,0.978124,0.977731,0.97658,0.974559,0.974806,0.971737,0.97625,0.973835,0},
      {0,0.97629,0.97684,0.974686,0.975646,0.978304,0.974812,0.976774,0.974989,0.97602,0.972343,0.976373,0.977416,0},
      {0,0.976764,0.973971,0.976367,0.974828,0.981619,0.972582,0.975314,0.974335,0.975078,0.974491,0.975129,0.97474,0},
      {0,0.977386,0.973411,0.97355,0.977815,0.975475,0.975502,0.973276,0.96899,0.972848,0.973578,0.975223,0.978876,0},
      {0,0.971509,0.971786,0.966537,0.963773,0.973926,0.971768,0.970993,0.976471,0.971795,0.969057,0.973723,0.97042,0},
      {0,0.97851,0.976437,0.973267,0.973721,0.975488,0.977317,0.97546,0.979721,0.971984,0.977242,0.976056,0.972045,0},
      {0,0.976754,0.974929,0.97199,0.974581,0.98108,0.975862,0.965369,0.960434,0.975137,0.980122,0.973272,0.978732,0},
      {0,0.974374,0.97867,0.975783,0.976854,0.981648,0.981342,0.974636,0.975137,0.978089,0.976642,0.978059,0.982274,0},
      {0,0.978015,0.96933,0.972354,0.965682,0.978435,0.970564,0.973064,0.976461,0.97343,0.977011,0.973579,0.970923,0},
      {0,0.969346,0.967401,0.976577,0.965642,0.973921,0.966087,0.971377,0.956679,0.971891,0.973357,0.971093,0.978112,0},
      {0,0.948198,0.96084,0.973659,0.959104,0.96371,0.972736,0.96976,0.971857,0.963366,0.970443,0.964358,0.970896,0},
      {0,0.935909,0.937143,0.957731,0.941935,0.961538,0.958797,0.944255,0.952278,0.969457,0.95811,0.961283,0.951639,0},
      {0,0.00981997,0.00646552,0.00671141,0.0121212,0.00976562,0.0105708,0.00625,0.0100806,0,0.0142857,0.0043956,0.00822368,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
};


  Int_t tmpPtBin = FindMuonEfficiencyBin( Pt , ptBins, 15);
  Int_t tmpEtaBin = FindMuonEfficiencyBin( Eta , etaBins, 16);
  Int_t tmpPhiBin = FindMuonEfficiencyBin( Phi , phiBins, 12);
  return Efficiency[tmpPtBin][tmpEtaBin][tmpPhiBin];
}




Double_t GetMuonEfficiencyPtEta(Double_t Pt, Double_t Eta) {

  Double_t ptBins[16] = {5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,50};
  Double_t etaBins[17] = {0,0.2,0.4,0.6,0.8,1,1.2,1.4442,1.566,1.8,2,2.1,2.2,2.3,2.4,2.5,2.6};


  Double_t Efficiency[17][18] = {
    {0,0.00310655,0.00376616,0.00442773,0.00460216,0.00609777,0.00624254,0.00698353,0.00577092,0.00599372,0.0067733,0.00803596,0.00755214,0.00809036,0.00720825,0.000518065,0,0    },
    {0,0.79081,0.799727,0.821869,0.839758,0.839475,0.836364,0.851779,0.857797,0.884295,0.889957,0.849263,0.801925,0.718383,0.578852,0.00235174,0,0    },
    {0,0.798185,0.813868,0.821194,0.84086,0.840193,0.860549,0.865077,0.878393,0.882582,0.884576,0.86695,0.819139,0.755767,0.642021,0.00297114,0,0    },
    {0,0.818603,0.823544,0.829957,0.856189,0.853911,0.857929,0.871391,0.879998,0.87099,0.880473,0.857226,0.816761,0.77707,0.657194,0.00689198,0,0    },
    {0,0.844785,0.85652,0.862426,0.87051,0.883154,0.885467,0.882724,0.891544,0.89731,0.901019,0.869883,0.850329,0.791095,0.675151,0.00202191,0,0    },
    {0,0.864884,0.868423,0.885853,0.879464,0.892678,0.901323,0.900814,0.907389,0.909827,0.910153,0.900489,0.873202,0.824297,0.74095,0.00427069,0,0    },
    {0,0.899748,0.90273,0.901473,0.916077,0.917461,0.922043,0.921247,0.928266,0.932186,0.932596,0.929714,0.905869,0.850771,0.77638,0.00635865,0,0    },
    {0,0.913176,0.923517,0.924026,0.935773,0.932712,0.936499,0.936406,0.945602,0.948164,0.948519,0.936907,0.924398,0.885065,0.814694,0.00880694,0,0    },
    {0,0.938205,0.935037,0.942054,0.945271,0.938845,0.945369,0.94818,0.951941,0.953122,0.959868,0.949132,0.937199,0.908581,0.857262,0.00651605,0,0    },
    {0,0.943991,0.945746,0.950374,0.950843,0.952598,0.953318,0.955402,0.957979,0.961,0.963173,0.958725,0.94898,0.922178,0.865236,0.00733988,0,0    },
    {0,0.951185,0.955679,0.956763,0.95545,0.956241,0.958642,0.960313,0.964707,0.96718,0.966832,0.965252,0.949272,0.93776,0.893001,0.00580332,0,0    },
    {0,0.950236,0.951173,0.953432,0.956132,0.954751,0.955091,0.956232,0.962074,0.960141,0.965906,0.960336,0.949169,0.939457,0.908669,0.00852713,0,0    },
    {0,0.957911,0.95772,0.960108,0.96135,0.96242,0.962072,0.960435,0.963349,0.965235,0.967116,0.966238,0.959255,0.950209,0.931981,0.0122542,0,0    },
    {0,0.962815,0.963788,0.96662,0.967082,0.967496,0.967654,0.965408,0.96831,0.969698,0.970974,0.971376,0.965849,0.956838,0.941952,0.0137029,0,0    },
    {0,0.969177,0.971359,0.969926,0.970385,0.970631,0.968544,0.968582,0.971928,0.972583,0.974459,0.971142,0.967553,0.957586,0.94859,0.00880902,0,0    },
    {0,0.956886,0.956622,0.959032,0.95875,0.957539,0.958023,0.954711,0.958679,0.959045,0.961807,0.957473,0.955189,0.95201,0.938766,0.00813665,0,0    },
    {0,0.958114,0.959787,0.961188,0.961276,0.960797,0.960194,0.956426,0.960943,0.959582,0.96318,0.958689,0.955747,0.950968,0.937852,0.00817563,0,0}
  };


  Int_t tmpPtBin = FindMuonEfficiencyBin( Pt , ptBins, 15);
  Int_t tmpEtaBin = FindMuonEfficiencyBin( Eta , etaBins, 16);
  return Efficiency[tmpPtBin][tmpEtaBin];
}
