#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <cstdio>
#include "math.h"
//#include "genAtomicState.cpp"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_coupling.h>
#include "DiagMatWithB.cpp"


using namespace std;

void printArray2D_complex(gsl_complex *array, int col, int row);
void printArray2D_real(gsl_complex *array, int col, int row); 
void printArray2D_imag(gsl_complex *array, int col, int row);
//double calc_gf(int F2, int J2, int I2,int L2, int S2,double g_I);
//double calc_gj(int J2, int L2, int S2);
bool isHermitian(gsl_complex *array, int size);

int main(int argc, char* argv[]) {

  bool production = true;
  double epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY * pow(10,27); // A^2 ns^4 / cm^3 g
  const double speed_of_light = GSL_CONST_CGSM_SPEED_OF_LIGHT * pow(10,-9); // cm/ns
  const double hbar = GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR    * pow(10,-9); // g cm^2 / ns
  printf("hbar = %8.6G\n",hbar);

  DiagMatWithB decomp;

  string isotope = "37K";
  int I2 = 3;
  double Aj[2];
  double g_I = 0.39094;
  int Je2 = 1;
  double tmax = 1.0; //ns
  double dt = 0.1; //ns
  double B_z = 2.0; //G
  double omega_Di = 2*M_PI*(speed_of_light / (769.9*pow(10,-7))) ; //ns^-1 omega = c/(2pi*lambda) and for Potassium the D1 laser has lambda = 769.9 nm.  

  if(argc > 1) {  
    if(strcmp(argv[1],"-h")==0) {
      cout << "First paramter is isotope (K37 or K38) [K38]" << endl;
      cout << "Second parameter is OP transition: 1 for D1 line 3 for D2 line [1]" << endl;
      cout << "Third parameter is tmax in ns [1]" << endl;
      cout << "Fourth paramter is B_ext in G [2]" << endl;
      cout << "Fifth paramter is dt in ns [0.1]" << endl;
      exit(0);
    }
    isotope = argv[1];
    if(argc > 2) {
      Je2 = atoi(argv[2]);
      if(argc > 3) {
	tmax = atof(argv[3]);
	if(argc > 4) {
	  B_z = atof(argv[4]);
	  if(argc > 5) { 
	    dt = atof(argv[5]);
	  }
	}
      }
    }
  } //End checking for input
  
    //Aj[i] and omega_Di are both in ns^-1 !!!
    //g_I is dimensionless
  if( strcmp(isotope.c_str(),"K37")  == 0 || strcmp(isotope.c_str(),"37K")  ==0 ) { 
    I2 = 3; 
    Aj[0] = 0.1201336 ; 
    Aj[1] = 0.01445  ; 
    g_I = .39094 ; 
    if(Je2 == 3) {
      omega_Di = 2*M_PI*(speed_of_light / (766.7 *pow(10,-7)));
      Aj[1] = 0.00805;
    }
  } //see Besch 1968 and Dan's thesis

  if( strcmp(isotope.c_str(),"K38")  == 0 || strcmp(isotope.c_str(),"38K")  ==0 ) { 
    I2 = 6; 
    Aj[0] = 0.70765   ; 
    Aj[1] = 0.0853   ; 
    g_I = 0.2029 ; 
    if(Je2 == 3) {
      omega_Di = 2*M_PI*(speed_of_light / (766.7 *pow(10,-7)));
      Aj[1] = 0.05575;
    }
  } //see Besch 1968 and Dan's thesis

  if( strcmp(isotope.c_str(),"K47")  == 0 || strcmp(isotope.c_str(),"47K")  ==0 ) { 
    I2 = 1; Aj[0] = 0.70765   ;
    Aj[1] = 0.0853   ; 
    g_I = 0.2029 ; 
    if(Je2 == 3) {
      omega_Di = 2*M_PI*(speed_of_light / (766.7 *pow(10,-7)));
      Aj[1] = 0.00805;
    }
  } //see Besch 1968 and Dan's thesis THIS IS NOT REAL NUMBRES, JUST NEEDED AN ISOTOPE WITH I = 1/2 FOR TESTING!
  
  
  //The idea is to adjust these parameters even though they will later be switched to better units
  //*********Inputs, constants*******
  double laser_I_fe = 0.100; //mW/cm^2
  double laser_I_ge = 0.100; //mW/cm^2
  double laser_pol[3] = {0.0,0.0,1.0};
  double detune = 0.0; //MHz
  double gamma_1 = 0.0; //Mhz = Linewidth of laser on the g->e transition omega_1 in T&J
  double gamma_2 = 0.0; //Mhz = Linewidth of laser on the f->e transition omega_2 in T&J
  //*********************************
  
  
  
  double I = (double)I2/2.0; 
  double Je = (double)Je2/2.0;
  double excitedShift = ((I + Je) * (I + Je + 1.0)) - ( I * (I+1.0)) - ( Je * (Je + 1.0));
  excitedShift *= 0.5 * Aj[1];
  detune *= pow(10,-3); //puts detune in ns^-1
  double omega_2 = omega_Di - (0.5*Aj[0]* ( ((I+0.5)*(I+1.5)) - ( I * (I+1.0)) - 0.75)) + excitedShift + detune;
  //double omega_1 = omega_Di  - (0.5*Aj[0]* ( ((I-0.5)*(I+0.5)) - ( I * (I+1.0)) - 0.75)) + excitedShift + detune;  //Frequency of laser tuned to the g->e transition (ns^-1)
  double omega_1 = omega_2 - 0.23746474;
  
  double tau = 26.2; // ns
  double gamma = 1.0/tau; //ns^-1
  double mu_B = 1.39962417998; //MHz/G = Bohr magneton
  
  
  
  //****Now switch all units from more natural units to ns, g, cm, A, G
  laser_I_fe *= pow(10,-23); // g/ns^3
  laser_I_ge *= pow(10,-23); // g/ns^3
  detune     *= pow(10,-3 ); //ns^-1
  gamma_1   *= pow(10,-3 ); //ns^-1
  gamma_2   *= pow(10,-3 ); //ns^-1
  mu_B       *= pow(10,-3 ); //ns^-1/G
  //*************************************************************


  double laser_E_fe = sqrt(laser_I_fe/(epsilon_0*speed_of_light)); // g cm^2 / A ns^3
  double laser_E_ge = sqrt(laser_I_ge/(epsilon_0*speed_of_light)); // g cm^2 / A ns^3
  
  
  //e States are any excited state either with J = 1/2, or 3/2 (D1 or D2)
  //f States are ground states with F = I + J
  //g States are ground states with F = I - J
  const int numEStates = (I2+1)*(Je2+1);
  const int numFStates = I2 + 2;
  const int numGStates = I2;


  if(!production) {
    cout << "************** OPTICAL PUMPING WITH OPTICAL-BLOCH EQUATIONS **************" << endl << endl;
    cout << "Pumping " << isotope << "\t As = " << Aj[0] << " ns^1 \t Ap = " << Aj[1] << " ns^-1 \t g-factor: " << g_I << endl;
    printf("omega_Di = %14.10G ns^-1\n",omega_Di);
    //This value is the D1 line in K37
    if(I2%2 == 0) { cout << endl << "Nuclear spin: " << I2/2 << endl; } 
    else{           cout << endl << "Nuclear spin: " << I2 << "/2" << endl;}
    cout << "tau: " <<  tau << " ns \t gamma = " << gamma << "n s^-1" << endl;
    cout << "Tmax: " << tmax << " ns.  Time step: " << dt << " ns." << endl;
    
    cout << "Magnetic field: " << B_z << " G" << endl;
    
    printf("Laser 1 frequency %12.10G ns^-1 \t linewidth = %4.2G ns^-1 \t Intensity = %G g/ns^3 \t E-field = %10.8G g*cm/(A*ns^3) \n", omega_1, gamma_1, laser_I_ge, laser_E_ge);
    printf("Laser 2 frequency %12.10G ns^-1 \t linewidth = %4.2G ns^-1 \t Intensity = %G g/ns^3 \t E-field = %10.8G g*cm/(A*ns^3) \n", omega_2, gamma_2, laser_I_fe, laser_E_fe);

    cout << "Number of excited states: " << numEStates << endl;
    cout << "f States with 2F = " << I2 + 1 << " -->  " << numFStates << endl;
    cout << "g States with 2F = " << I2 - 1 << " -->  " << numGStates << endl;

    //cout << "Epsilon_0 [MKS] = " << GSL_CONST_MKSA_VACUUM_PERMITTIVITY << endl; //A^2 s^4/(m^3 kg)
    //cout << "Speed of light [CGS] = " << GSL_CONST_CGSM_SPEED_OF_LIGHT << endl; //cm/s
    //cout << "Hbar [cgs] = " << GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR << endl; //g cm^2 /s
  } //Print out parameters  
  
    
  //Stuff to get the I,J decomposition
  double groundDecomp[numFStates+numGStates][numFStates+numGStates];
  decomp.diagH(I2, 0, 1 , mu_B, B_z, 0.0, g_I, Aj[0], &groundDecomp[0][0]);

  double excitedDecomp[numEStates][numEStates];
  decomp.diagH(I2, 1, Je2,mu_B, B_z, 0.0, g_I, Aj[1], &excitedDecomp[0][0]);
  

  bool debugDecomp = false;
  
  if(!production && debugDecomp) {
    printf("\n\n IN OPTICALBLOCH:\n");
    for(int i = 0; i < numEStates; i++) {
      for(int j = 0; j < numEStates; j++) {
	printf("%10.6G\t ",excitedDecomp[i][j]);
      }
      printf("\n");
    }
    printf("\n");
    for(int i = 0; i < numFStates+numGStates; i++) {
      for(int j = 0; j < numFStates+numGStates; j++) {
	printf("%10.6G\t ",groundDecomp[i][j]);
      }
      printf("\n");
    }
  }
  //END DECOMPOSITION


  //Getting an ordered list of Iz (same order that the helper class uses)
  double Izground[numFStates+numGStates];
  double Izexcited[numFStates+numGStates];
  for(int i = 0; i < numFStates + numGStates; i++)  Izground[i]  = -((double)I2/2.0) + (double)(i%(I2+1));
  for(int i = 0; i < numEStates             ; i++)  Izexcited[i] = -((double)I2/2.0) + (double)(i%(I2+1)); 


  //Set up output
  const int maxTOut = 1000; //Maxiumum number of time steps to write to a file
  int printFreq; //Number of ns between print statements
  int totalPrint; //Total number of print statements
  if(tmax <= maxTOut) {
    printFreq = 1;
    totalPrint = tmax;
  }
  else {
    printFreq = tmax / maxTOut;
    printf("printFreq = %i\n",printFreq);
    if(printFreq == 0) printFreq = 1;
    totalPrint = maxTOut;
  }
  //printFreq = 0.1;
  FILE * file;
  printf("printFreq = %i\n",printFreq);

  if(production) {
    file = fopen("opData.dat","w");
    stdout = file;
    fprintf(file,"%i \t %i \t %6.2G \t %i \t %i \t %i 4.0\n",numEStates+numFStates+numGStates,numEStates,1.0,totalPrint,I2,Je2);
   
  }
  else file = stdout;
  double updateFreq = max(dt * 1000.0,1000.0);

  FILE * eeFile = fopen("eeData.dat","w");
  FILE * ffFile = fopen("ffData.dat","w");
  FILE * ggFile = fopen("ggData.dat","w");
  FILE * efFile = fopen("efData.dat","w");
  FILE * egFile = fopen("egData.dat","w");
  FILE * fgFile = fopen("fgData.dat","w");
  //End set up output

  gsl_complex rhoee[numEStates][numEStates];
  gsl_complex rhoff[numFStates][numFStates];
  gsl_complex rhogg[numGStates][numGStates];
  
  gsl_complex zrhoef[numEStates][numFStates];
  gsl_complex zrhoeg[numEStates][numGStates];
  gsl_complex zrhofg[numFStates][numGStates];
  
  gsl_complex drhoee[numEStates][numEStates];
  gsl_complex drhoff[numFStates][numFStates];
  gsl_complex drhogg[numGStates][numGStates];

  gsl_complex dzrhoef[numEStates][numFStates];
  gsl_complex dzrhoeg[numEStates][numGStates];
  gsl_complex dzrhofg[numFStates][numGStates];

  //Fill populations to starting vales
  for(int i = 0; i < numEStates; i++) {
      for(int j = 0; j < numEStates; j++) {
	double fill = 0.0;
	double fillAll = 1.0/ ( (double)numFStates+numGStates);
	if( i == j) {
	  //fill = 1.0/( (double)(numFStates+numGStates));
	  fill = 1.0/ ( (double)numFStates+numGStates);
	}
	//rhoee[i][j] = gsl_complex_rect(0.0,0.0);
	rhoee[i][j] = gsl_complex_rect(0.0, 0.0);
	drhoee[i][j] = gsl_complex_rect(0.0,0.0);
	
	if(i < numFStates && j < numFStates) {
	  //rhoff[i][j] = gsl_complex_rect(fill,0.0);
	  rhoff[i][j] = gsl_complex_rect(fill,0.0);
	  drhoff[i][j] = gsl_complex_rect(0.0,0.0);
	}
	if(j < numFStates) {
	  zrhoef[i][j] = gsl_complex_rect(0.0,0.0);
	  //zrhoef[i][j] = gsl_complex_rect(fillAll,fillAll);
	  dzrhoef[i][j] = gsl_complex_rect(0.0,0.0);
	}
	if(i < numGStates && j < numGStates) {
	  //rhogg[i][j] = gsl_complex_rect(fill,0.0);
	  rhogg[i][j] = gsl_complex_rect(fill,0.0);
	  drhogg[i][j] = gsl_complex_rect(0.0,0.0);
	} 
	if(j < numGStates) {
	  //zrhoeg[i][j] = gsl_complex_rect(fill,fill);
	  zrhoeg[i][j] = gsl_complex_rect(0.0,0.0);
	  dzrhoeg[i][j] = gsl_complex_rect(0.0,0.0);
	}	
	if(i < numFStates && j < numGStates) {
	  zrhofg[i][j] = gsl_complex_rect(0.0,0.0);
	  dzrhofg[i][j] = gsl_complex_rect(0.0,0.0);
	}


      }
   }
  
  
  //Setup ordered list of eigenvalues to access by index
  int Fe2Array[numEStates];
  int MFe2Array[numEStates];
  
  int Ff2Array[numFStates];
  int MFf2Array[numFStates];

  int Fg2Array[numGStates];
  int MFg2Array[numGStates];

  int Fe2 = abs(I2 - Je2);
  int MFe2 = -Fe2;
  for(int e = 0; e < numEStates; e++) {
    Fe2Array[e] = Fe2;
    MFe2Array[e] = MFe2;
    //printf("e = %i \t Fe2 = %i \t Mfe2 = %i\n",e,Fe2Array[e],MFe2Array[e]);
    MFe2 += 2;
    if(MFe2 > Fe2) {
      Fe2 += 2;
      MFe2 = -Fe2;
    }
  }
  
  int Ff2 = (I2 + 1);
  int MFf2 = -Ff2;
  for(int f = 0; f < numFStates; f++) {
    Ff2Array[f] = Ff2;
    MFf2Array[f] = MFf2;
    //printf("f = %i \t Ff2 = %i \t Mff2 = %i\n",f,Ff2Array[f],MFf2Array[f]);
    MFf2 += 2;
  }

  int Fg2 = (I2 - 1);
  int MFg2 = -Fg2;
  for(int g = 0; g < numGStates; g++) {
    Fg2Array[g] = Fg2;
    MFg2Array[g] = MFg2;
    //printf("g = %i \t Fg2 = %i \t Mfg2 = %i\n",g,Fg2Array[g],MFg2Array[g]);
    MFg2 += 2;
  }


 
  //Array to hold the relative energy of each state E = 0 is at the S1/2 level before the hyperfine splitting!!!
  double omegaE[numEStates]; //Ordering of states: F = I - J --> F = I + J ;;; Mf = -F --> Mf = F
  double omegaF[numFStates]; 
  double omegaG[numGStates];
  
  

  //Get all energies for all states
  bool debugEnergy = true;
  int fStep = 0;
  int gStep = 0;
  for(int i = 0; i < numEStates; i++)   {
    
    double I = ((double)I2)/2.0;
    double J_excited = ((double)Je2/2.0);
    double F_excited = ((double)Fe2Array[i])/2.0;
    
    //Excited states
    omegaE[i] = omega_Di;//Fine splitting [ns^-1]
    omegaE[i] += 0.5*Aj[1]*( (F_excited * (F_excited+1.0)) - (I * (I+1.0)) - (J_excited * (J_excited+1.0))); //Hyperfine [ns^-1]
    double g_Fe = decomp.calc_gf(Fe2, Je2, I2, 2, 1, g_I);
    omegaE[i] += g_Fe * mu_B * B_z * (double)MFe2Array[i] / 2.0; //[ns-1]
    //if(!production && debugEnergy) printf(" F = %3.1f \t Mf = %i/2 \t omegaE = %12.8f ns^-1 \n",F_excited,MFe2Array[i],omegaE[i]);
    
    //F states
    if( i < numFStates) {

      double F_f = (double)Ff2Array[fStep]/2.0;
      omegaF[fStep] = 0.5*Aj[0]*( (F_f * (F_f + 1.0)) - (I * (I+1.0)) - 0.75); //[ns^-1]
      double g_Ff = decomp.calc_gf(Ff2Array[fStep], 1, I2, 0, 1,g_I);
      omegaF[fStep]+= g_Ff * mu_B * B_z * (double)MFf2Array[fStep] / 2.0; // ns^-1
      if(!production && debugEnergy) printf(" F = %3.1f \t Mf = %i/2 \t omegaE = %12.8f ns^-1 \n",F_f,MFf2Array[fStep],omegaF[fStep]);
      fStep++;
      
    }
    //G states
    if( i < numGStates) {
      double F_g = (double)Fg2Array[gStep]/2.0;
      omegaG[gStep] = 0.5*Aj[0]*( (F_g * (F_g + 1.0)) - (I * (I+1.0)) - 0.75); //ns^-1
      double g_Fg = decomp.calc_gf(Fg2Array[gStep], 1, I2, 0, 1, g_I);
      omegaG[gStep] += g_Fg * mu_B * B_z * (double)MFg2Array[gStep] /2.0;  //ns^-1
      if(!production && debugEnergy) printf(" F = %3.1f \t Mf = %i/2 \t omegaE = %12.8f ns^-1 \n",F_g,MFg2Array[gStep],omegaG[gStep]);
      gStep++;
    }
    
  }



  //Coupling constants:  Equation 17 of T&J
  //Third index is polarization of the light
  const double emitProb = 3*M_PI*epsilon_0*hbar*pow(speed_of_light,3)/tau; //cm^2 A^2 /ns
  //printf("Constant factor in Eq 31 = %12.5G cm^2 A^2 /ns \n",emitProb);
  double aef[numEStates][numFStates][3]; //[]
  double aeg[numEStates][numGStates][3]; //[]
  
  double Def[numEStates][numFStates][3]; //[A * cm * s ]
  double Deg[numEStates][numGStates][3]; //[A * cm * s ]

  bool debugCoupling = false;
  for(int i = 0; i < numEStates; i++) {  
    for(int j = 0; j < numFStates; j++) {
      for(int q = 0; q < 3; q++) {
	double I = ((double)I2) / 2.0;
	double J_e = ((double)Je2) / 2.0;
	double F_e = ((double)Fe2Array[i]) / 2.0;
	double F_f = ((double)Ff2Array[j]) / 2.0;
	double MFe = ((double)MFe2Array[i]) / 2.0;
	
	//***************Trembaly and Jacques Eq 17 ********************
	aef[i][j][q] = pow(-1.0, 1.0 + I + J_e + F_e + F_f - MFe);
	aef[i][j][q] *= (sqrt( (2*F_f) + 1.0) * sqrt( (2*F_e) + 1.0) * sqrt( (2*J_e) + 1.0));
	aef[i][j][q] *= gsl_sf_coupling_3j( Fe2Array[i] , 2 , Ff2Array[j] , -MFe2Array[i], (q-1)*2 , MFf2Array[j]);
	aef[i][j][q] *= gsl_sf_coupling_6j( Fe2Array[i] , 2 , Ff2Array[j] , 1            , I2      , Je2         );
	//****************************************************************


	double omega_ef = fabs(omegaE[i] - omegaF[j]); //ns^-1

	Def[i][j][q] = aef[i][j][q] * sqrt(emitProb/pow(omega_ef,3)); // ns*cm*A
	if(!production && debugCoupling && fabs(aef[i][j][q]) > 0.0) printf(" |%3.1f, %3.1f > ---> |%3.1f, %3.1f > with q = %i \t a = %8.6G\t dOmega = %12.10G \t D = %8.6G\n",F_f,(double)MFf2Array[j]/2.0,F_e,(double)MFe2Array[i]/2.0,q-1,aef[i][j][q],omega_ef,Def[i][j][q]);

	if(j < numGStates) {
	  double F_g = ((double)Fg2Array[j]) / 2.0;
	  //***************Trembaly and Jacques Eq 17 ********************
	  aeg[i][j][q] = pow(-1.0, 1.0 + I + J_e + F_e + F_g - MFe);
	  aeg[i][j][q] *= (sqrt( (2*F_g) + 1.0) * sqrt( (2*F_e) + 1.0) * sqrt( (2*J_e) + 1.0));
	  aeg[i][j][q] *= gsl_sf_coupling_3j( Fe2Array[i] , 2 , Fg2Array[j] , -MFe2Array[i], (q-1)*2 , MFg2Array[j]);
	  aeg[i][j][q] *= gsl_sf_coupling_6j( Fe2Array[i] , 2 , Fg2Array[j] , 1            , I2      , Je2         );
	
	  //**************************************************************

	  double omega_eg = fabs(omegaE[i] - omegaG[j]); //ns^-1
	  Deg[i][j][q] = aeg[i][j][q] * sqrt(emitProb/pow(omega_eg,3)); //ns cm A
	  if(!production && debugCoupling && fabs(aeg[i][j][q]) > 0.0) printf(" |%3.1f, %3.1f > ---> |%3.1f, %3.1f > with q = %i \t a = %8.6G\t dOmega = %12.10G \t D = %8.6G\n",F_g,(double)MFg2Array[j]/2.0,F_e,(double)MFe2Array[i]/2.0,q-1,aeg[i][j][q],omega_eg,Deg[i][j][q]);
		  
	}
      }
    }
  }
  

  double startingPop = 0.0;

  for(int i = 0; i < numEStates; i++) {
    startingPop += GSL_REAL(rhoee[i][i]);
	if(i < numFStates) startingPop += GSL_REAL(rhoff[i][i]);
	if(i < numGStates) startingPop += GSL_REAL(rhogg[i][i]);
  }
  cout << "Starting pop = " << startingPop << endl;
      
  bool writeAll = true;

  if(writeAll) { 
    fprintf(eeFile,"%i \t %i \t %6.2G \t %i \t %i \t %i 4.0\n",numEStates+numFStates+numGStates,numEStates,1.0,totalPrint,I2,Je2);
    fprintf(ffFile,"%i \t %i \t %6.2G \t %i \t %i \t %i 4.0\n",numEStates+numFStates+numGStates,numEStates,1.0,totalPrint,I2,Je2);
    fprintf(ggFile,"%i \t %i \t %6.2G \t %i \t %i \t %i 4.0\n",numEStates+numFStates+numGStates,numEStates,1.0,totalPrint,I2,Je2);
    fprintf(efFile,"%i \t %i \t %6.2G \t %i \t %i \t %i 4.0\n",numEStates+numFStates+numGStates,numEStates,1.0,totalPrint,I2,Je2);
    fprintf(egFile,"%i \t %i \t %6.2G \t %i \t %i \t %i 4.0\n",numEStates+numFStates+numGStates,numEStates,1.0,totalPrint,I2,Je2);
    fprintf(fgFile,"%i \t %i \t %6.2G \t %i \t %i \t %i 4.0\n",numEStates+numFStates+numGStates,numEStates,1.0,totalPrint,I2,Je2);
  }
  
  //************************************************************************************************
  //***********START OF TIME LOOP*******************************************************************
  //************************************************************************************************
  double nextPrint = 0.0;
  double update = 0.0;
  for(double time = 0; time < tmax; time += dt) {  

    if(production && fabs(time - update) < pow(10,-2)) {
      cout << time << endl;
      update = update + updateFreq;

    }
    //Print stuff out
    bool printMatrices = false;
    double totalPopulation = 0.0;
    for(int i = 0; i < numEStates; i++) {
      totalPopulation                    += GSL_REAL(rhoee[i][i]);
      if(i < numFStates) totalPopulation += GSL_REAL(rhoff[i][i]);
      if(i < numGStates) totalPopulation += GSL_REAL(rhogg[i][i]);
    }
    //cout << "time = " << time << " nextPrint = " << 
    if(printMatrices  && !production )  {
      fprintf(file,"rhoee : \n");
      printArray2D_complex(&rhoee[0][0],numEStates,numEStates);
      fprintf(file,"\n rhoff: \n");
      printArray2D_complex(&rhoff[0][0],numFStates,numFStates);
      fprintf(file,"\n rhogg : \n");
      printArray2D_complex(&rhogg[0][0],numGStates,numGStates);
      fprintf(file,"\n zrhoef : \n");
      printArray2D_complex(&zrhoef[0][0],numEStates,numFStates);
      fprintf(file,"\n zrhoeg : \n");
      printArray2D_complex(&zrhoeg[0][0],numEStates,numGStates);
      fprintf(file,"\n zrhofg : \n");
      printArray2D_complex(&zrhofg[0][0],numFStates,numGStates);
      fprintf(file,"___________________________________________________________________________________________________________________________________\n");
      fprintf(file,"\n"); 
      
    } else if( fabs(time - nextPrint) < pow(10,-2) ) {
      nextPrint = nextPrint + printFreq;
      //cout << "hit on time = " << time << " now nextPrint = " << nextPrint << endl;

      fprintf(file,"%12.6G ",time);     
      for(int g = 0; g < numGStates; g++) {
	fprintf(file,"%12.6G ",GSL_REAL(rhogg[g][g]));
      }
      for(int f = 0; f < numFStates; f++) {
	fprintf(file,"%12.6G ",GSL_REAL(rhoff[f][f]));
      }
      for(int e = 0; e < numEStates; e++) {
	fprintf(file,"%12.6G ",GSL_REAL(rhoee[e][e]));
      }
      
    
      
      fprintf(file,"%12.8G ",totalPopulation);


      //Get the nuclear polarization
      bool debugPolarization = false;
      double nucPol = 0.0;
      double align = 0.0;
      for(int i = 0; i < (numGStates + numFStates); i++) {	
	if(debugPolarization && !production) printf("i = %i\n",i);
	for( int fzindex = 0; fzindex < (numFStates + numGStates); fzindex++) {
	  if( i < numGStates) {
	    nucPol += Izground[fzindex] * GSL_REAL(rhogg[i][i]) * groundDecomp[fzindex][i] * groundDecomp[fzindex][i];
	    align += Izground[fzindex] * Izground[fzindex] * GSL_REAL(rhogg[i][i]) * groundDecomp[fzindex][i] * groundDecomp[fzindex][i];
	    if(debugPolarization && !production) printf("\tIz = %1.3G, Decomp = %8.6G, Population = %8.6G, Pol = %8.6G, Align = %8.6G \n",Izground[fzindex],groundDecomp[fzindex][i], GSL_REAL(rhogg[i][i]), nucPol,align);
	  } else      {
	    nucPol += Izground[fzindex] * GSL_REAL(rhoff[i-numGStates][i-numGStates]) * groundDecomp[fzindex][i] * groundDecomp[fzindex][i];
	    align += Izground[fzindex]*Izground[fzindex] * GSL_REAL(rhoff[i-numGStates][i-numGStates]) * groundDecomp[fzindex][i] * groundDecomp[fzindex][i];
	    if(debugPolarization && !production) printf("\tIz = %1.3G, Decomp = %8.6G, Population = %8.6G, Pol = %8.6G, Align = %8.6G \n",Izground[fzindex],groundDecomp[fzindex][i], GSL_REAL(rhoff[i][i]), nucPol,align);
	  }

	}
      }
      for(int e = 0; e < numEStates; e++) {
	for( int fzindex = 0; fzindex < numEStates; fzindex++) {
	  nucPol += Izexcited[fzindex] * GSL_REAL(rhoee[e][e]) * excitedDecomp[fzindex][e] * excitedDecomp[fzindex][e];
	  align += Izexcited[fzindex] * Izexcited[fzindex] * GSL_REAL(rhoee[e][e]) * excitedDecomp[fzindex][e] * excitedDecomp[fzindex][e];
	}
      }

      double I = ((double)I2)/2.0;
      nucPol /= I;
      nucPol /= totalPopulation;
      fprintf(file,"%12.6G ",nucPol);
      
      align /= totalPopulation;
      align = (I * (I + 1.0)) - (3 * align);
      align /= I * (2.0*I - 1.0);
      fprintf(file,"%12.6G ",align);
    

      fprintf(file,"\n"); 

      if(writeAll) {
	fprintf(eeFile,"%12.6G\t",time);
	for(int e = 0; e < numEStates; e++) {
	  for(int ep = 0; ep < numEStates; ep++) {
	    fprintf(eeFile,"%12.6G\t %12.6G\t",GSL_REAL(rhoee[e][ep]),GSL_IMAG(rhoee[e][ep]));
	  }
	}
	fprintf(eeFile,"\n");
      
      
      
	fprintf(ffFile,"%12.6G\t",time);
	for(int f = 0; f < numFStates; f++) {
	  for(int fp = 0; fp < numFStates; fp++) {
	    fprintf(ffFile,"%12.6G\t %12.6G\t",GSL_REAL(rhoff[f][fp]),GSL_IMAG(rhoee[f][fp]));
	  }
	}
	fprintf(ffFile,"\n");
      

      	fprintf(ggFile,"%12.6G\t",time);
	for(int g = 0; g < numGStates; g++) {
	  for(int gp = 0; gp < numGStates; gp++) {
	    fprintf(ggFile,"%12.6G\t %12.6G\t",GSL_REAL(rhogg[g][gp]),GSL_IMAG(rhoee[g][gp]));
	  }
	}
	fprintf(ggFile,"\n");
      
      	fprintf(efFile,"%12.6G\t",time);
	for(int e = 0; e < numEStates; e++) {
	  for(int f = 0; f < numFStates; f++) {
	    fprintf(efFile,"%12.6G\t %12.6G\t",GSL_REAL(zrhoef[e][f]),GSL_IMAG(zrhoef[e][f]));
	  }
	}
	fprintf(efFile,"\n");
            
      	fprintf(fgFile,"%12.6G\t",time);
	for(int f = 0; f < numFStates; f++) {
	  for(int g = 0; g < numGStates; g++) {
	    fprintf(fgFile,"%12.6G\t %12.6G\t",GSL_REAL(zrhofg[f][g]),GSL_IMAG(zrhofg[f][g]));
	  }
	}
	fprintf(fgFile,"\n");
      
      	fprintf(egFile,"%12.6G\t",time);
	for(int e = 0; e < numEStates; e++) {
	  for(int g = 0; g < numGStates; g++) {
	    fprintf(egFile,"%12.6G\t %12.6G\t",GSL_REAL(zrhoeg[e][g]),GSL_IMAG(zrhoeg[e][g]));
	  }
	}
	fprintf(egFile,"\n");
      
      }
      

    }

    for(int i = 0; i < numGStates; i ++) {
      for(int f = 0; f < numGStates; f++) {
	//printf("zrhoeg[%i][%i] = %10.6G + %10.6Gi\n",i,f,GSL_REAL(rhogg[i][f]),GSL_IMAG(rhogg[i][f]));
      }
    }


    if(fabs (startingPop - totalPopulation) > pow(10,-8)) {
      cout << "PARTICLE NUMBER NOT CONSERVED.  ABORTING.  time = " << time << endl;
      exit(1);
    }

    if(!isHermitian(&rhoee[0][0],numEStates) || !isHermitian(&rhoff[0][0], numFStates) || !isHermitian(&rhogg[0][0], numGStates)) {
      printf("DENSITY MATRIX NOT HERMITIAN!!!\n");
      printf("time = %.2G\n",time);
      exit(1);
    }

    bool debugGLaser_ee = false;
    bool debugFLaser_ee = false;
    bool debugSpontDecay_ee = false;
    
    //Excited state populations (Eq 32)
    for(int e = 0; e < numEStates; e++)  { 
      for(int ep = 0; ep < numEStates; ep++) {

	//Very important - check to make sure time step is small enough (the delta omega terms are a time scale of the wiggles!!!)
	if( time < 2.0 && (1.0/fabs(omegaE[e] - omegaE[ep]) < dt*10)) {
	  cout << "ERROR:  TIME STEP NOT SMALL ENOUGH.  omegaE[" << e << "]-omegaE[" << ep << "] = " << omegaE[e]-omegaE[ep] << " ns^-1 ---> " << 1.0/(omegaE[e]-omegaE[ep]) << "ns time scale " << endl;
	  exit(1);
	}
	drhoee[e][ep] = gsl_complex_rect(-gamma,-(omegaE[e]-omegaE[ep]));
	//drhoee[e][ep] = gsl_complex_rect(-gamma,0.0);
	drhoee[e][ep] = gsl_complex_mul(drhoee[e][ep],rhoee[e][ep]);
	if(!production && debugSpontDecay_ee && time < 1) printf("e = %i, ep = %i, gamma = %8.6G, deltaOmega = %8.6G\n",e,ep,gamma,omegaE[e]-omegaE[ep]);
	
	
	//Now for the laser interaction!!!
	gsl_complex gLaserTerm_ee = gsl_complex_rect(0.0,0.0);
	gsl_complex fLaserTerm_ee = gsl_complex_rect(0.0,0.0);
	for(int iq = 0; iq < 3; iq ++) { //Dot product term
	  //if(!production && debugGLaser_ee) fprintf(file,"\t iq = %i\n",iq);
	  for (int gpp = 0; gpp < numGStates; gpp++) { //sum over gpp in (T&J 32)
	    
	    //gsl_complex temp = gsl_complex_mul( gsl_complex_rect( Deg[e][gpp][iq], 0.0), gsl_complex_conjugate(zrhoeg[ep][gpp]));
	    gsl_complex temp  = gsl_complex_mul_real( gsl_complex_conjugate(zrhoeg[ep][gpp]), pow(1.0,iq-1)*Deg[e][gpp][iq]);

	    if(!production && debugGLaser_ee && iq == 2) fprintf(file,"\t Deg[%i][%i][%i] = %8.6G  \t ",e,gpp,iq,Deg[e][gpp][iq]);
	    if(!production && debugGLaser_ee && iq == 2) fprintf(file,"\t left-temp = %12.9G \t + %12.9G i ", GSL_REAL(temp), GSL_IMAG(temp));
	    
	    gsl_complex temp2 = gsl_complex_mul( gsl_complex_rect( Deg[ep][gpp][iq], 0.0), zrhoeg[e][gpp]);
	    if(!production && debugGLaser_ee && iq == 2) fprintf(file,"\t right-temp = %12.9G \t + %12.9G i ", GSL_REAL(temp2), GSL_IMAG(temp2));

	    temp = gsl_complex_sub(temp,temp2);
	    temp = gsl_complex_mul_real(temp,laser_pol[iq]);
	    if(!production && debugGLaser_ee && iq == 2) fprintf(file,"\t after polarization = %12.9G \t + %12.9G i \n", GSL_REAL(temp), GSL_IMAG(temp));

	    gLaserTerm_ee = gsl_complex_add(gLaserTerm_ee, temp);
	  }
	} //End dot product for laser one
	gLaserTerm_ee = gsl_complex_mul_real( gLaserTerm_ee, laser_E_ge);

	//Laser 2 term goes here!
	
	for(int iq = 0; iq < 3; iq++) { //Dot prodcut term
	  if(!production && debugFLaser_ee) fprintf(file,"\t iq = %i\n",iq);
	  for(int fpp = 0; fpp < numFStates;fpp++) { //sum over fpp in T&J 32
	    gsl_complex temp  = gsl_complex_mul( gsl_complex_rect( Def[e][fpp][iq],0.0) , gsl_complex_conjugate(zrhoef[ep][fpp]));
	    if(!production && debugFLaser_ee) fprintf(file,"\t left-temp = %12.9G \t + %12.9G i ", GSL_REAL(temp), GSL_IMAG(temp));
	    gsl_complex temp2 = gsl_complex_mul( gsl_complex_rect( Def[ep][fpp][iq],0.0), zrhoef[e][fpp]);
	    if(!production && debugFLaser_ee) fprintf(file,"\t right-temp = %12.9G \t + %12.9G i ", GSL_REAL(temp2), GSL_IMAG(temp2));
	    temp = gsl_complex_sub(temp,temp2);
	    temp = gsl_complex_mul_real(temp,laser_pol[iq]);
	
	    if(!production && debugFLaser_ee) fprintf(file,"\t after polarization = %12.9G \t + %12.9G i \n", GSL_REAL(temp), GSL_IMAG(temp));
	    fLaserTerm_ee = gsl_complex_add(fLaserTerm_ee, temp);

	  }
	} //end dot product loop for laser 2
	fLaserTerm_ee = gsl_complex_mul_real( fLaserTerm_ee, laser_E_fe);
	
	
	gLaserTerm_ee = gsl_complex_add(gLaserTerm_ee, fLaserTerm_ee);
	gLaserTerm_ee = gsl_complex_div_imag( gLaserTerm_ee, -2.0*hbar);
	drhoee[e][ep] = gsl_complex_add(drhoee[e][ep],gLaserTerm_ee);
	
	
	drhoee[e][ep] = gsl_complex_mul_real(drhoee[e][ep],dt);
	
      }
    } //END EXCITED STATE POPULATIONS

    
    //F-Ground state populations (Eq 33)
    bool debugSpontDecay_ff = false;
    bool debugZCoherences_ff = false;
    for(int f = 0; f < numFStates; f++) {  
      for(int fp = 0; fp < numFStates; fp++) {
	
	drhoff[f][fp] = gsl_complex_rect(0.0,0.0);
	
	//**********************
	//SPONTANEOUS DECAY TERM
	//**********************
	if(!production && debugSpontDecay_ff) fprintf(file,"f = %i   fp = %i   m_f = %i/2   m_fp = %i/2 \n",f,fp,MFf2Array[f],MFf2Array[fp]);
	for(int e = 0; e < numEStates; e++) {
	  for(int ep = 0; ep < numEStates; ep++) {
	    
	    if( (MFe2Array[e] - MFe2Array[ep]) == (MFf2Array[f] - MFf2Array[fp])) {
	      
	      int qef   = (MFe2Array[e]  - MFf2Array[f]) / 2;
	      int qepfp = (MFe2Array[ep] - MFf2Array[fp]) / 2;
	      if( abs(qef) <= 1 && abs(qepfp) <=1) {
	
		double temp = aef[e][f][qef+1] * aef[ep][fp][qepfp+1];
		gsl_complex spontDecay = gsl_complex_mul_real(rhoee[e][ep],temp);
		drhoff[f][fp] = gsl_complex_add(drhoff[f][fp],spontDecay);
		if(!production && debugSpontDecay_ff) fprintf(file,"\t\t spontDecay = %6.4f + %14.10G i \t",GSL_REAL(spontDecay),GSL_IMAG(spontDecay));
		if(!production && debugSpontDecay_ff) fprintf(file,"\t\t drhoff[%i][%i] = %6.4f + %6.4G i \n", f,fp,GSL_REAL(drhoff[f][fp]), GSL_IMAG(drhoff[f][fp]));
	      }
	    }
	   
	  } //end ep loop
	} //end e loop
	
	drhoff[f][fp] = gsl_complex_mul(drhoff[f][fp], gsl_complex_rect(gamma,0.0));


	//Very important - check to make sure time step is small enough (the delta omega terms are a time scale of the wiggles!!!)
	if(time < 2.0 && 1.0/fabs(omegaF[f] - omegaF[fp]) < dt*10) {
	  cout << "ERROR:  TIME STEP NOT SMALL ENOUGH.  omegaF[" << f << "]-omegaF[" << fp << "] = " << omegaF[f]-omegaF[fp] << " ns^-1 ---> " << 1.0/(omegaF[f]-omegaF[fp]) << "ns time scale " << endl;
	  exit(1);
	}

	gsl_complex zCoherencesF = gsl_complex_mul_imag(rhoff[f][fp], (omegaF[f] - omegaF[fp]));
	if(!production && debugZCoherences_ff) printf("zCoherence[%i][%i] = %8.6G + %8.6G i \n",f,fp,GSL_REAL(zCoherencesF),GSL_IMAG(zCoherencesF));
	drhoff[f][fp] = gsl_complex_sub( drhoff[f][fp], zCoherencesF);
      
	
	//********END SPONTANEOUS DECAY**************
	
	//EF LASER PUMPING TERM
	gsl_complex laserTerm_ff = gsl_complex_rect(0.0,0.0);
	for(int iq = 0; iq < 3; iq++) {
	  for(int epp = 0; epp < numEStates; epp++) {
	    gsl_complex temp = gsl_complex_mul(  gsl_complex_rect( Def[epp][f][iq] ,0.0), zrhoef[epp][fp]);
	    gsl_complex temp2 = gsl_complex_mul( gsl_complex_rect( Def[epp][fp][iq],0.0), gsl_complex_conjugate(zrhoef[epp][f]));
	    temp = gsl_complex_sub(temp,temp2);
	    temp = gsl_complex_mul_real(temp,laser_pol[iq]);

	    laserTerm_ff = gsl_complex_add(laserTerm_ff,temp);
	  }
	}
	laserTerm_ff = gsl_complex_mul_imag(laserTerm_ff,laser_E_fe/(2.0*hbar));
	drhoff[f][fp] = gsl_complex_add(drhoff[f][fp], laserTerm_ff);
       	//END EF LASER PUMPING TERM
	
	drhoff[f][fp] = gsl_complex_mul(drhoff[f][fp], gsl_complex_rect(dt,0.0));
	
      }//End Fp loop
      
    } //END F-GROUND STATE POPULATIONS


    //G-Ground State populations (Eq 34)
    for(int g = 0; g < numGStates; g++) {
      for(int gp = 0; gp < numGStates; gp++) {

	drhogg[g][gp] = gsl_complex_rect(0.0,0.0);
	
	//*************************
	//***Spontaneous Decay*****
	//*************************
	for(int e = 0; e < numEStates; e++) {
	  for(int ep = 0; ep < numEStates; ep++) {
	    
	    if( (MFe2Array[e] - MFe2Array[ep]) == (MFg2Array[g] - MFg2Array[gp])) {
	      int qeg = (MFe2Array[e] - MFg2Array[g]) / 2;
	      int qepgp = (MFe2Array[ep] - MFg2Array[gp]) / 2;
	      if( abs(qeg) <= 1 && abs(qepgp) <=1) {

		double temp = aeg[e][g][qeg+1] * aeg[ep][gp][qepgp+1];
		gsl_complex spontDecay = gsl_complex_mul_real(rhoee[e][ep],temp);
		drhogg[g][gp] = gsl_complex_add(drhogg[g][gp],spontDecay);
	      }
	    }

	  } //end ep loop
	} //end e loop
	//END SPONTANEOUS EMISSION!!!
	drhogg[g][gp] = gsl_complex_mul(drhogg[g][gp],gsl_complex_rect(gamma,0.0));


	//Very important - check to make sure time step is small enough (the delta omega terms are a time scale of the wiggles!!!)
	if(time < 2.0 && ( 1.0/fabs(omegaG[g] - omegaG[gp]) < dt*10)) {
	  cout << "ERROR:  TIME STEP NOT SMALL ENOUGH.  omegaG[" << g << "]-omegaG[" << gp << "] = " << omegaG[g]-omegaG[gp] << " ns^-1 ---> " << 1.0/(omegaG[g]-omegaG[gp]) << "ns time scale " << endl;
	  exit(1);
	}

	gsl_complex zCoherencesG = gsl_complex_mul_imag(rhogg[g][gp],(omegaG[g] - omegaG[gp]));
	drhogg[g][gp] = gsl_complex_sub( drhogg[g][gp], zCoherencesG);

	
	//EG LASER PUMPING TERM
	gsl_complex laserTerm_gg = gsl_complex_rect(0.0,0.0);
	for(int iq = 0; iq < 3; iq++) {
	  for(int epp = 0; epp < numEStates; epp++) {
	    gsl_complex temp = gsl_complex_mul(  gsl_complex_rect( pow(1.0,iq-1)*Deg[epp][g][iq] ,0.0), zrhoeg[epp][gp]);
	    gsl_complex temp2 = gsl_complex_mul( gsl_complex_rect( Deg[epp][gp][iq],0.0), gsl_complex_conjugate(zrhoeg[epp][g]));
	    temp = gsl_complex_sub(temp,temp2);
	    temp= gsl_complex_mul_real(temp,laser_pol[iq]);
	    laserTerm_gg = gsl_complex_add(laserTerm_gg,temp);
	  }
	}
	laserTerm_gg = gsl_complex_mul_imag(laserTerm_gg,laser_E_ge/(2.0*hbar));
	drhogg[g][gp] = gsl_complex_add(drhogg[g][gp],laserTerm_gg);
	//END EG LASER PUMPING TERM
	
	drhogg[g][gp] = gsl_complex_mul(drhogg[g][gp],gsl_complex_rect(dt   ,0.0));

      } //End Gp loop
    } //END G-GROUND STATE POPULATIONS

    
    //EF-Coherences (Eq 35)
    bool debugDetune_ef = false;
    bool debugLaser2_ef = false;
    for(int e = 0; e < numEStates; e++) {
      for(int f = 0; f < numFStates; f++) {


	double real = (gamma+gamma_2) / 2.0;
	double imag = (omegaE[e] - omegaF[f]) - omega_2;
	if(time < 2.0 && (1.0/fabs(real) < dt*10 || 1.0/fabs(imag) < dt*10)) {
	  cout << "ERROR:  TIME STEP NOT SMALL ENOUGH FOR EF COHERENCES" << endl;
	  cout << "Real = " << real << " ns^-1    Imag = " << imag << " ns^-1" << endl;
	  exit(1);
	}

	if(!production && debugDetune_ef) fprintf(file,"e = %i \t f = %i \t Average linewidth = %4.2G ns^-1 \t Detuning = %12.10G \n",e,f,real,imag);
	
	dzrhoef[e][f] = gsl_complex_mul( gsl_complex_rect(-real, -imag), zrhoef[e][f]);
	if(!production && debugDetune_ef) fprintf(file,"\t dzrhoef[%i][%i] = %8.6G + %8.6G i \n",e,f,GSL_REAL(dzrhoef[e][f]),GSL_IMAG(dzrhoef[e][f]));

       
	//Now to the Laser_2 term
	gsl_complex laser2_total = gsl_complex_rect(0.0,0.0);
	
	if(!production && debugLaser2_ef) fprintf(file,"e = %i   f = %i\n",e,f);

	for(int q = 0; q < 3; q++) {
	  if(!production && debugLaser2_ef) fprintf(file,"\t q = %i \n",q);

	  gsl_complex laser2_ff = gsl_complex_rect(0.0,0.0);
	  gsl_complex laser2_ee = gsl_complex_rect(0.0,0.0);

	  for(int fpp = 0; fpp < numFStates; fpp++) {
	    gsl_complex temp = gsl_complex_mul_real(rhoff[fpp][f],Def[e][fpp][q]);
	    laser2_ff = gsl_complex_add(laser2_ff, temp);
	    if(!production && debugLaser2_ef) fprintf(file,"\t\t fpp = %i \t Def[%i][%i][%i] = %10.6G    and contribution is = %10.6G + %10.6G i \n", fpp, e,fpp,q,Def[e][fpp][q],GSL_REAL(temp), GSL_IMAG(temp));
	  }
	  for(int epp = 0; epp < numEStates; epp++) {
	    gsl_complex temp = gsl_complex_mul_real(rhoee[e][epp],Def[epp][f][q]);
	    laser2_ee = gsl_complex_add(laser2_ee, temp);
	    if(!production && debugLaser2_ef) fprintf(file,"\t\t epp = %i \t and contribution is = %10.6G + %10.6G i \n", epp, GSL_REAL(temp), GSL_IMAG(temp));
	  }
	
	  laser2_total = gsl_complex_sub(laser2_ff,laser2_ee);
	  laser2_total = gsl_complex_mul_imag(laser2_total, laser_E_fe*laser_pol[q]/(2.0*hbar));
	  if(!production && debugLaser2_ef) fprintf(file,"\t laser2_total = %10.6G + %10.6G i\n", GSL_REAL(laser2_total), GSL_IMAG(laser2_total));
	  dzrhoef[e][f] = gsl_complex_add(dzrhoef[e][f],laser2_total);
	}
	//End laser 2 term;
	
	//Now the laser 1 term

	//LASER 1 TERM
	gsl_complex laser1_total = gsl_complex_rect(0.0,0.0);
	for(int q = 0; q < 3; q++) {
	  gsl_complex qsum = gsl_complex_rect(0.0,0.0);
	  for(int gpp = 0; gpp < numGStates; gpp++) {
	    gsl_complex temp = gsl_complex_mul_real( gsl_complex_conjugate(zrhofg[f][gpp]), Deg[e][gpp][q]);
	    qsum = gsl_complex_add(qsum,temp);
	  }
	  qsum = gsl_complex_mul_real(qsum,laser_E_ge*laser_pol[q]);
	  laser1_total = gsl_complex_add(laser1_total,qsum);
	}
	laser1_total = gsl_complex_div_imag(laser1_total,-2.0*hbar);
	dzrhoef[e][f] = gsl_complex_add(dzrhoef[e][f],laser1_total);
	//End laser 1 term 
	
	
	dzrhoef[e][f] = gsl_complex_mul_real(dzrhoef[e][f], dt);
      }
    } //END EF-COHERENCES


    
    //EG-Coherences (Eq 36)
    for(int e = 0; e < numEStates; e++) {
      for(int g = 0; g < numGStates; g++) {
	
	double real = (gamma +gamma_1) / 2.0;
	double imag = (omegaE[e] - omegaG[g]) - omega_1;
	if(time < 2.0 && (1.0/fabs(real) < dt*10 || 1.0/fabs(imag) < dt*10)) {
	  cout << "ERROR:  TIME STEP NOT SMALL ENOUGH FOR EG COHERENCES" << endl;
	  exit(1);
	}


	dzrhoeg[e][g] = gsl_complex_mul( gsl_complex_rect(-real,-imag),zrhoeg[e][g]);

	//Now to the Laser 1 term
	
	gsl_complex laser1_total = gsl_complex_rect(0.0,0.0);
	for(int q = 0; q < 3; q++) {
	  gsl_complex laser1_gg = gsl_complex_rect(0.0,0.0);
	  gsl_complex laser1_ee = gsl_complex_rect(0.0,0.0);
	  for(int gpp = 0; gpp < numGStates; gpp ++) {
	    gsl_complex temp = gsl_complex_mul_real(rhogg[gpp][g],Deg[e][gpp][q]);
	    laser1_gg = gsl_complex_add(laser1_gg,temp);
	  }
	  for(int epp = 0; epp < numEStates; epp++) {
	    gsl_complex temp = gsl_complex_mul_real(rhoee[e][epp],Deg[epp][g][q]);
	    laser1_ee = gsl_complex_add(laser1_ee, temp);
	  }
	  laser1_total = gsl_complex_sub(laser1_gg,laser1_ee);
	  laser1_total = gsl_complex_mul_imag(laser1_total,laser_E_ge*laser_pol[q]/(2.0*hbar));
	  dzrhoeg[e][g] = gsl_complex_add(dzrhoeg[e][g],laser1_total);
	}
	//End laser 1 term

	//Laser 2 term
	
	bool debugEGLaser2 = false;
	if( debugEGLaser2 && !production) printf(" e = %i, g = %i\n",e,g);
	gsl_complex laser2_total = gsl_complex_rect(0.0,0.0);
	
	for(int q = 0; q < 3; q++) {
	  if( debugEGLaser2 && !production) printf("\tq = %i\n",q);
	  gsl_complex qsum = gsl_complex_rect(0.0,0.0);
	  for(int fpp = 0; fpp < numFStates; fpp++) {
	    gsl_complex temp = gsl_complex_mul_real(gsl_complex_conjugate(zrhofg[fpp][g]), Def[e][fpp][q]);
	    qsum = gsl_complex_add(qsum,temp);
	    if( debugEGLaser2 && !production) printf("\t\tzrhofg[%i][%i] = %8.4G, Def[%i][%i][%i] = %8.4G, Running sum = %8.4G\n",fpp,g,GSL_REAL(zrhofg[fpp][g]),e,fpp,q,Def[e][fpp][q],GSL_REAL(qsum));
	  }
	  qsum = gsl_complex_mul_real(qsum,laser_E_fe*laser_pol[q]);
	  laser2_total = gsl_complex_add(laser2_total,qsum);
	}
	//End laser 2 term
	laser2_total = gsl_complex_div_imag(laser2_total,-2.0*hbar);
	dzrhoeg[e][g] = gsl_complex_add(dzrhoeg[e][g],laser2_total);
	
	dzrhoeg[e][g] = gsl_complex_mul_real( dzrhoeg[e][g]           , dt);
      }
    } //END EG-COHERENCES

    
    //FG-Coherences (Eq 37)
    bool debugFG = false;
    bool debugFGLaser = false;
    for(int f = 0; f < numFStates; f++) {
      for(int g = 0; g < numGStates; g++) {

	
	double real = (gamma_1 + gamma_2) / 2.0;
	double imag = ( (omegaF[f] - omegaG[g]) - (omega_1 - omega_2));
	if(time < 2.0 && (1.0/fabs(real) < dt*10 || 1.0/fabs(imag) < dt*10)) {
	  cout << "ERROR:  TIME STEP NOT SMALL ENOUGH FOR FG COHERENCES" << endl;
	  exit(1);
	}


	dzrhofg[f][g] = gsl_complex_rect( - real, -imag);
	dzrhofg[f][g] = gsl_complex_mul(dzrhofg[f][g], zrhofg[f][g]);
		
	
	
	//Now for the two laser terms
	//Laser 1
	gsl_complex laser1_total = gsl_complex_rect(0.0,0.0);
	for(int q = 0; q < 3; q++) {
	  if( debugFGLaser && !production) printf("\t q = %i\n",q);
	  gsl_complex qsum = gsl_complex_rect(0.0,0.0);
	  for(int epp = 0; epp < numEStates; epp++) {
	    gsl_complex temp = gsl_complex_mul_real( gsl_complex_conjugate(zrhoef[epp][f]), Deg[epp][g][q]);
	    qsum = gsl_complex_add(qsum,temp);
	    if( debugFGLaser && !production) printf("\t\t epp = %i, zrhoef* = %8.4G + %8.4Gi, Deg = %8.4G, qsum = %8.4G + %8.4Gi\n",epp,GSL_REAL(zrhoef[epp][f]),-GSL_IMAG(zrhoef[epp][f]),Deg[epp][g][q],GSL_REAL(qsum),GSL_IMAG(qsum));
	  }
	  qsum = gsl_complex_mul_real(qsum, laser_E_ge*laser_pol[q]);
	  laser1_total = gsl_complex_add(laser1_total,qsum);
	}
	

	//Laser 2
	
	gsl_complex laser2_total = gsl_complex_rect(0.0,0.0);
	
	for(int q = 0; q < 3; q++) {
	  gsl_complex qsum = gsl_complex_rect(0.0,0.0);
	  for(int epp = 0; epp < numEStates; epp++) {
	    gsl_complex temp = gsl_complex_mul_real( zrhoeg[epp][g], Def[epp][f][q]);
	    qsum = gsl_complex_add(qsum, temp);
	  }
	  qsum = gsl_complex_mul_real(qsum, laser_E_fe*laser_pol[q]);
	  laser2_total = gsl_complex_add(laser2_total,qsum);
	}
	
	if( debugFGLaser && !production) printf("\t\t laser1-total = %8.4G +%8.4G i \t laser2-total = %8.4G + %8.4G i\n",GSL_REAL(laser1_total),GSL_IMAG(laser1_total),GSL_REAL(laser2_total),GSL_IMAG(laser2_total));
	
	gsl_complex laser_total = gsl_complex_sub(laser2_total,laser1_total);
	laser_total = gsl_complex_div_imag( laser_total, -2.0*hbar);
	dzrhofg[f][g] = gsl_complex_add(dzrhofg[f][g],laser_total);
		
	dzrhofg[f][g] = gsl_complex_mul_real(dzrhofg[f][g], dt);
	
      }
    
    } //END FG COHERENCES
    

    
    //Update with all the new populations
    for(int i = 0; i < numEStates; i++) {
      for(int j = 0; j < numEStates; j++) { 
	if(i < numEStates && j < numEStates)  rhoee[i][j] = gsl_complex_add( rhoee[i][j], drhoee[i][j]);
	if(i < numFStates && j < numFStates)  rhoff[i][j] = gsl_complex_add( rhoff[i][j], drhoff[i][j]);
	if(i < numGStates && j < numGStates)  rhogg[i][j] = gsl_complex_add( rhogg[i][j], drhogg[i][j]);
	if(i < numEStates && j < numFStates) zrhoef[i][j] = gsl_complex_add(zrhoef[i][j],dzrhoef[i][j]);
	if(i < numEStates && j < numGStates) zrhoeg[i][j] = gsl_complex_add(zrhoeg[i][j],dzrhoeg[i][j]);
	if(i < numFStates && j < numGStates) zrhofg[i][j] = gsl_complex_add(zrhofg[i][j],dzrhofg[i][j]);
	
      }
    }

    for(int e = 0; e < numEStates; e++) {
      //printf("drhoee[%i][%i] = %12.6G + %12.6Gi\n",e,e,GSL_REAL(drhoee[e][e]),GSL_IMAG(rhoee[e][e]));
    }
    for(int f= 0; f < numFStates; f++) {
      //printf("drhoff[%i][%i] = %12.6G + %12.6Gi\n",f,f,GSL_REAL(drhogg[f][f]),GSL_IMAG(rhogg[f][f]));
    }

  }  //END TIME LOOP
  
} 

void printArray2D_complex(gsl_complex *array, int row, int col) { 
  

  gsl_complex * value;
  value = new gsl_complex[col*row];
  value = array;
  //printf(" col = %i row = %i\n",col,row);
  for(int i = 0; i < col*row; i++) {
    //printf(" value[%i] = %8.6G + %8.6G i\n",i,GSL_REAL(value[i]),GSL_IMAG(value[i]));
  }
  
  int ir = 0;
  int ii = 0;
  for(int r = 0; r < row; r++) {

    for(int c = 0; c < col; c++) {
      fprintf(stdout,"%10.6G   ",GSL_REAL(value[ir]));
      ir++;
    }
    fprintf(stdout,"\t");
    for(int c = 0; c < col; c++) {
      fprintf(stdout,"%10.6G i   ",GSL_IMAG(value[ii]));
      ii++;
    }
    fprintf(stdout,"\n");

  }
}

void printArray2D_real(gsl_complex *array, int col, int row) {
  gsl_complex value;
  
  for( int c = 0; c < col; c++) {
    for(int r = 0; r < row; r++) {
      value = *(array+c+(col*r));
      //      cout << setw(10) << left << GSL_REAL(value);
      fprintf(stdout,"%12.8G",GSL_REAL(value));
    }
    fprintf(stdout,"\n");
  }
  
}

void printArray2D_imag(gsl_complex *array, int col, int row) {  
  
  gsl_complex value;
  
  for( int c = 0; c < col; c++) {
    for(int r = 0; r < row; r++) {
      value = *(array+c+(col*r));
      //      cout << setw(10) << left << GSL_REAL(value);
      fprintf(stdout,"%12.8G ",GSL_IMAG(value));
    }
    fprintf(stdout,"\n");
  }
  
}

bool isHermitian(gsl_complex *array, int size) {

  bool hermit = true;
  double eps = pow(10,-8);
  gsl_complex * value;
  value = new gsl_complex[size*size];
  value = array;

  //First check that the diagonals are real
  for(int i = 0 ; i < size; i++) {
    if( fabs(GSL_IMAG(value[i + (i*size)])) > eps) hermit = false;
    //printf("Checking for i = %i; Imag = %8.4G\n",i,GSL_IMAG(value[i+ (i*size)]));
  }

  gsl_complex upper;
  gsl_complex lower;
  
  for(int c = 0; c < size; c++) {
    for(int r = 0; r < size; r++) {
      upper = value[c + (size*r)];
      lower = value[r + (size*c)];
      //printf("for c = %i, r = %i; Lower = %8.4G + %8.4G \t Upper = %8.4G + %8.4G \t",c,r,GSL_REAL(lower),GSL_IMAG(lower),GSL_REAL(upper),GSL_IMAG(upper));
      double ImSum = fabs( GSL_IMAG(upper) + GSL_IMAG(lower));
      double ReDif = fabs( GSL_REAL(upper) - GSL_REAL(lower));

      //printf("Dif of reals: %8.4G \t Sum of imaginaries: %8.4G \n", ReDif,ImSum);
      if( ImSum > eps) hermit = false;
      if( ReDif > eps) hermit = false;

    }
  }
  return hermit;
}
