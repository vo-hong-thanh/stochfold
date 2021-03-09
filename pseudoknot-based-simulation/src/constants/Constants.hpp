
#pragma once

#include <vector>
#include <limits>




namespace EnergyConstants {

   constexpr double INF =  std::numeric_limits<double>::max();
   //stochfold
   constexpr double  R = 1.98717E-3;   // kcal/(K*mol)
   constexpr double T = 310.5;
   constexpr double k0 = 1e4;
   constexpr int maxLoopSize = 30;
   constexpr double lxc37 = 1.07856; // R*T*1.75 at 37 C

  //turner 2004 and 1999
  //AU end
  constexpr double AU_end = 0.45;
  //GU end
  constexpr double GU_end = 0.45;

  //turner 2004 and 1999
  inline const std::vector<std::vector<double>> stem {              /* CG      GC      GU      UG      AU     UA   */  
                        /*(i, j)*/      /*CG*/ std::vector<double> { -3.26,  -2.36,  -1.41,  -2.11,  -2.11,  -2.08  },
                        /* |  */        /*GC*/ std::vector<double> { -3.42,  -3.26,  -1.53,  -2.51,  -2.35,  -2.24  },  
                        /* v  */        /*GU*/ std::vector<double> { -2.51,  -2.11,   0.47,   1.29,  -1.27,  -1.36  },
                                        /*UG*/ std::vector<double> { -1.53,  -1.41,   0.30,   0.47,  -1.00,  -0.55  },
                                        /*AU*/ std::vector<double> { -2.24,  -2.08,  -0.55,  -1.36,  -0.93,  -1.10  },
                                        /*UA*/ std::vector<double> { -2.35,  -2.11,  -1.00,  -1.27,  -1.33,  -0.93  }
                                        };
  //Turener 1999 
  inline const std::vector<double> multiloop { 10.1, -0.3, -0.3 };

  // idx = 0 special C bulge (2004) | idx > 0 idx = nmbr of unpaired bases
  // turner 2004 and 1999
  inline const std::vector<double> bulge { -0.9, 3.80, 2.80, 3.20, 3.60, 4.00, 4.40, 4.60, 4.70, 4.80, 4.90, 
                                             5.00, 5.10, 5.20, 5.3, 5.4, 5.4, 5.5, 5.5, 5.60, 5.7, 5.7,
                                             5.8, 5.80, 5.8, 5.9, 5.9, 6.00, 6.00, 6.0, 6.1
                                           };
  //hairpin 
  //turner 2004
  inline const std::vector<double> hairpin04 = { INF, INF, INF, 5.40, 5.60, 5.70, 5.40, 6.00, 5.50, 6.40, 6.5,
                                               6.6, 6.7, 6.8, 6.9, 6.9, 7.0, 7.1, 7.1, 7.2, 7.2,
                                               7.3, 7.3, 7.4, 7.4, 7.5, 7.5, 7.5, 7.6, 7.6, 7.7
                                             };
  //turner 1999                                         
  inline const std::vector<double> hairpin = { INF, INF, INF, 5.70, 5.60, 5.60, 5.40, 5.90, 5.60, 6.40, 6.5,
                                               6.6, 6.7, 6.8, 6.9, 6.9, 7.0, 7.1, 7.1, 7.2, 7.2,
                                               7.3, 7.3, 7.4, 7.4, 7.5, 7.5, 7.5, 7.6, 7.6, 7.7
                                             };


   
  
    
      


  
  // UU or GA first mismatch 
  //turner 2004
  constexpr double hairpinUUGAmm04 = -0.9;
  // turner 1999
  constexpr double hairpinUUGAmm = -0.8;

  // GG first mismatch 
   //turner 2004
  constexpr double hairpinGGmm04 = -0.8;

  // special GU closure 2004 and 1999
  constexpr double hairpinGUclosure = -2.2;

  // c3 loop
  //turner 2004
  constexpr double hairpinC304 = 1.5;
  // turner 1999
  constexpr double hairpinC3 = 1.4;
  //all-C loop A    2004 and 1999
  constexpr double hairpinAllCA = 0.3;
  //all-C loop B  2004 and 1999
  constexpr double hairpinAllCB = 1.6;

  //special hairpin loops
  //turner 2004
  inline const std::vector<double>  hairpinSpecial04 = {
                                          /*CAACG*/	    6.8
                                          /*GUUAC*/,	  6.9
                                          /*CUACGG*/,	  2.8
                                          /*CUCCGG*/,	  2.7
                                          /*CUUCGG*/,	  3.7
                                          /*CUUUGG*/,   3.7
                                          /*CCAAGG*/,	  3.3
                                          /*CCCAGG*/,	  3.4
                                          /*CCGAGG*/,	  3.5
                                          /*CCUAGG*/,	  3.7
                                          /*CCACGG*/,	  3.7
                                          /*CCGCGG*/,	  3.6
                                          /*CCUCGG*/,	  2.5
                                          /*CUAAGG*/,	  3.6
                                          /*CUCAGG*/,	  3.7
                                          /*CUUAGG*/,	  3.5
                                          /*CUGCGG*/,	  2.8
                                          /*CAACGG*/,	  5.5
                                          /*ACAGUGCU*/,	2.9
                                          /*ACAGUGAU*/,	3.6
                                          /*ACAGUGUU*/,	1.8
                                          /*ACAGUACU*/,	2.8
 

  };
// turner 1999
 inline const  std::vector<double> hairpinTetraloop = {
                                                /*CGGAAG*/ -3.0, /*CUUCGG*/ -3.0, /*CGUGAG*/ -3.0,/*CGAAAG*/ -3.0, /*CGCAAG*/ -3.0, /*CGAAGG*/ -2.5, 
                                                /*CUACGG*/ -2.5, /*CGCGAG*/ -2.5, /*CGUAAG*/ -2.0, /*CUAACG*/ -2.0, /*CGAGAG*/ -2.0,/*CGGGAG*/ -1.5,
                                                
                                                /*GGGGAC*/ -3.0, /*GGUGAC*/ -3.0,  /*GGAGAC*/ -3.0,  /*GGAAAC*/ -3.0,  /*GGCAAC*/ -2.5,	/*GGAAGC*/ -1.5,
                                                /*GGCGAC*/ -1.5, /*GGGAGC*/ -1.5, /*GUGAAC*/ -1.5, /*GGGAAC*/ -1.5,
    
                                                /*UGAGAG*/ -2.5,/*UGAAAG*/ -2.0,  /*UGAAAA*/ -1.5, /*UGGAAA*/ -1.5,

                                                /*AGAAAU*/ -2.0,/*AGCAAU*/ -1.5, /*AGUAAU*/ -1.5,	/*AGUGAU*/ -1.5

                                            };

                                          //  /*GGGGAC*/ -3.0, /*GGUGAC*/ -3.0, /*CGAAAG*/ -3.0, /*GGAGAC*/ -3.0, /*CGCAAG*/ -3.0, /*GGAAAC*/ -3.0,		
                                           //     /*CGGAAG*/ -3.0, /*CUUCGG*/ -3.0, /*CGUGAG*/ -3.0, /*CGAAGG*/ -2.5, /*CUACGG*/ -2.5, /*GGCAAC*/ -2.5,	
                                            //    /*CGCGAG*/ -2.5, /*UGAGAG*/ -2.5, /*CGAGAG*/ -2.0, /*AGAAAU*/ -2.0, /*CGUAAG*/ -2.0, /*CUAACG*/ -2.0,		
                                             //   /*UGAAAG*/ -2.0, /*GGAAGC*/ -1.5, /*GGGAAC*/ -1.5, /*UGAAAA*/ -1.5, /*AGCAAU*/ -1.5, /*AGUAAU*/ -1.5,	
                                               // /*CGGGAG*/ -1.5, /*AGUGAU*/ -1.5, /*GGCGAC*/ -1.5, /*GGGAGC*/ -1.5, /*GUGAAC*/ -1.5, /*UGGAAA*/ -1.5

//turner 2004 terminal mismatch
  inline const std::vector<std::vector<double>> terminal_mismatch04 = {
               /*XY -> */
                                    /* AA    AC    AG    AU    CA    CC    CG    CU    GA    GC    GG    GU    UA    UC    UG    UU   */  
       /*AU*/  std::vector<double> { -0.8, -1.0, -0.8, -1.0, -0.6, -0.7, -0.6, -0.7, -0.8, -1.0, -0.8, -0.1, -0.6, -0.8, -0.6, -0.8  },
       /*CG*/  std::vector<double> { -1.5, -1.5, -1.4, -1.5, -1.0, -1.1, -1.0, -0.8, -1.4, -1.5, -1.6, -1.5, -1.0, -1.4, -1.0, -1.2  },  
       /*GC*/  std::vector<double> { -1.1, -1.5, -1.3, -1.5, -1.1, -0.7, -1.1, -0.5, -1.6, -1.5, -1.4, -1.5, -1.1, -1.0, -1.1, -0.7  },
       /*GU*/  std::vector<double> { -0.3, -1.0, -0.8, -1.0, -0.6, -0.7, -0.6, -0.7, -0.6, -1.0, -0.8, -1.0, -0.6, -0.8, -0.6, -0.6  },
       /*UA*/  std::vector<double> { -1.0, -0.8, -1.1, -0.8, -0.7, -0.6, -0.7, -0.5, -1.1, -0.8, -1.2, -0.8, -0.7, -0.6, -0.7, -0.5  },
       /*UG*/  std::vector<double> { -1.0, -0.8, -1.1, -0.8, -0.7, -0.6, -0.7, -0.8, -0.5, -0.8, -0.8, -0.8, -0.7, -0.6, -0.7, -0.5  }
                                        };
  //1999                                     
  inline const std::vector<std::vector<double>> terminal_mismatch = {
               /*XY -> */
                                    /* AA    AC    AG    AU    CA    CC    CG    CU    GA    GC    GG    GU    UA    UC    UG    UU   */  
       /*AU*/  std::vector<double> { -0.8, -1.0, -0.8, -0.4, -0.6, -0.7, -1.7, -0.7, -0.8, -1.6, -0.8, -0.1, -0.6, -0.8, -0.9, -0.8  },
       /*CG*/  std::vector<double> { -1.5, -1.5, -1.4, -1.6, -1.0, -1.1, -2.8, -0.8, -1.4, -1.9, -1.6, -0.9, -1.6, -1.4, -1.6, -1.2  },  
       /*GC*/  std::vector<double> { -1.1, -1.5, -1.3, -1.9, -1.1, -0.7, -2.9, -0.5, -1.6, -2.8, -1.4, -1.0, -1.7, -1.0, -2.0, -0.7  },
       /*GU*/  std::vector<double> { -0.3, -1.0, -0.8, -0.8, -0.6, -0.7, -2.0, -0.7, -0.6, -1.6, -0.8,  0.0, -0.9, -0.8, -1.8, -0.8  },
       /*UA*/  std::vector<double> { -1.0, -0.8, -1.1, -0.8, -0.7, -0.6, -1.9, -0.5, -1.1, -1.6, -1.2, -0.5, -0.4, -0.6, -0.8, -0.5  },
       /*UG*/  std::vector<double> { -1.0, -0.8, -1.1, -0.5, -0.7, -0.6, -1.0, -0.5, -0.5, -0.9, -0.8, -0.8, -0.1, -0.6, -0.0, -0.5  }
                                        };                                
//turner 1999 
  inline const std::vector<std::vector<double>>  internal_loop_11 = {
                                     /* AA    AC    AG    AU    CA    CC    CG    CU    GA    GC    GG    GU    UA    UC    UG    UU   */  
      /*AU-AU*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.5  },
      /*AU-CG*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.0  },
      /*AU-GC*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*AU-UA*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.2  },
      /*AU-GU*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      /*AU-UG*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      
      /*CG-AU*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*CG-CG*/ std::vector<double> {  0.4, -0.4,  0.4,  0.4,  0.3,  0.5,  0.4,  0.5, -0.1,  0.4, -1.7,  0.4,  0.4,  0.0,  0.4, -0.3  },
      /*CG-GC*/ std::vector<double> {  1.1,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4, -1.4,  0.4,  0.4,  0.4,  0.4,  0.4  },
      /*CG-UA*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*CG-GU*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*CG-UG*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      
      /*GC-AU*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*GC-CG*/ std::vector<double> {  0.8,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4, -2.1,  0.4,  0.4,  0.4,  0.4, -0.7  },
      /*GC-GC*/ std::vector<double> {  0.4,  0.3, -0.1,  0.4, -0.4,  0.5,  0.4,  0.0,  0.4,  0.4, -1.7,  0.4,  0.4,  0.5,  0.4, -0.3  },
      /*GC-UA*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.0  },
      /*GC-GU*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*GC-UG*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      
      /*UA-AU*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.8  },
      /*UA-CG*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*UA-GC*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*UA-UA*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.5  },
      /*UA-GU*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      /*UA-UG*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      
      /*GU-AU*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      /*GU-CG*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*GU-GC*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*GU-UA*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      /*GU-GU*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      /*GU-UG*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      
      /*UG-AU*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      /*UG-CG*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*UG-GC*/ std::vector<double> {  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1,  1.1, -1.0,  1.1,  1.1,  1.1,  1.1,  1.1  },
      /*UG-UA*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      /*UG-GU*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  },
      /*UG-UG*/ std::vector<double> {  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7,  1.7, -0.4,  1.7,  1.7,  1.7,  1.7,  1.7  }  
    };

        //internal loop  1999 idx 2 and 3 formula 
     inline const std::vector<double> internal_loop = { INF, INF, 4.10, 5.10, 1.70, 1.80, 2.00, 2.20, 2.30, 2.40, 2.50,
                                                     2.60, 2.70, 2.80, 2.9, 3.00, 3.00, 3.10, 3.10, 3.20, 3.30,
                                                     3.30, 3.40, 3.40, 3.40, 3.5, 3.50, 3.6, 3.60, 3.60, 3.70
                                                  };
    // turner 2004
   constexpr double internalAsymmetry04 =0.6;   
   constexpr double internalAUGUClosure04 = 0.7;
   //                                            AG    GA    GG    UU
   inline const std::vector<double> internalmm04{ -0.8, -1.0, -1.2, -0.7 };    

  // truener 1999
   constexpr double internalAsymmetry = 0.48;
    constexpr double internalAUGUClosure = 0.65;
     //                                            AG    GA    GG    UU
    inline const std::vector<double> internalmm { -1.1, -1.1,  0.0, -0.7};
                       
  
  //dangling end 
  //turner 1999 and 2004
                                                  /* 5' jX 3'
                                                     3' i  5'  */
     inline const std::vector<std::vector<double>> dangling_3  = {/*  X=   A     C     G     U   */  
                                            //ij
                                            /*AU*/ std::vector<double> { -0.7, -0.1, -0.7, -0.1  },
                                            /*CG*/ std::vector<double> { -1.1, -0.4, -1.3, -0.6  },  
                                            /*GC*/ std::vector<double> { -1.7, -0.8, -1.7, -1.2  },
                                            /*GU*/ std::vector<double> { -0.7, -0.1, -0.7, -0.1  },
                                            /*UA*/ std::vector<double> { -0.8, -0.5, -0.8, -0.6  },
                                            /*UG*/ std::vector<double> { -0.8, -0.5, -0.8, -0.6  }
                                           };
    
                                                  /*X -> */
    inline const std::vector<std::vector<double>> dangling_5  = {      /*  A     C     G     U   */  
                                           /*AU*/ std::vector<double> { -0.3, -0.3, -0.4, -0.2  },
                                           /*CG*/ std::vector<double> { -0.5, -0.3, -0.2, -0.1  },
                                           /*GC*/ std::vector<double> { -0.2, -0.3,  0.0,  0.0  },
                                           /*GU*/ std::vector<double> { -0.3, -0.3, -0.4, -0.2  },
                                           /*UA*/ std::vector<double> { -0.3, -0.1, -0.2, -0.2  },
                                           /*UG*/ std::vector<double> { -0.3, -0.1, -0.2, -0.2  }
                                           };
//used in Multiloop energy and pkenergy
    inline const std::vector<std::vector<double>> pseudoknot {
                                     //                              βX       βmul      βpsudo    β2        β3
                                    /*H-type*/  std::vector<double> {9.6,     15.0,     15.0,     0.1,      0.1 },
                                    /*K-type*/  std::vector<double> {12.6,    18.0,     18.0,     0.1,      0.1 },
                                    /*L-type*/  std::vector<double> {14.6,    20.0,     20.0,     0.1,      0.1 },
                                    /*M-type*/  std::vector<double> {17.6,    23.0,     23.0,     0.1,      0.1 }
    };


  //DP09 Params  IN USE 
     inline const std::vector<double> pseudoknotDP09 {
       
 /* initiating external */  -1.38
 /* initiating in multiloop*/,  10.07
 /* initiating in pseudoknot */,15.00
 /* initiating a band */ ,2.46
 /* unpaired base in a pseudoknot*/, 0.06
 /* nested substructure in a pseudoknot */,  0.96
 /* initiating a multiloop that spans a band */, 3.41
 /* branch in a multiloop that spans a band*/,  0.56
 /* for an unpaired base in a multiloop that spans a band */,  0.12
 /* Multiplicative for a stacked pair that spans a band */, 0.89
/* Multiplicative for an internal loop that spans a band */, 0.74
     };
    
    inline const std::vector<double> pseudoknotDP03 {
       
 /* initiating external  */ 9.60,
 /*  initiating in multiloop */15.00, 
 /* initiating  in pseudoknot */15.00,  
 /* initiating a bandY */0.20, 
 /* unpaired base in a pseudoknot */0.10,
 /* nested substructure in a pseudoknot */0.10,  
 /* initiating a multiloop that spans a band */3.40, 
 /* branch in a multiloop that spans a band*/0.40,  
 /* for an unpaired base in a multiloop that spans a band */0.00,  
 /* Multiplicative for a stacked pair that spans a band */0.83, 
/* Multiplicative for an internal loop that spans a bandN*/0.83
     };
};