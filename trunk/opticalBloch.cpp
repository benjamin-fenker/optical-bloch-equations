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
  
  const double epsilon_0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY * pow(10,27); // A^2 ns^4 / cm^3 g
  const double speed_of_light = GSL_CONST_CGSM_SPEED_OF_LIGHT * pow(10,-9); // cm/ns
  const double hbar = GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR    * pow(10,-9); // g cm^2 / ns

  DiagMatWithB decomp;

  string isotope = "37K";
  int I2 = 3;
  double Aj[2];
  double g_I = 0.39094;
  int Je2 = 1;
  double tmax = 1.0; //ns
  double dt = 1.0; //ns
  double B_z = 2.0; //G
  double omega_Di = 2*M_PI*(speed_of_light / (769.9*pow(10,-7))) ; //ns^-1 omega = c/(2pi*lambda) and for Potassium the D1 laser has lambda = 769.9 nm.  
  //printf("omega_Di = %12.8G ns^-1",omega_Di);
  //This value is the D1 line in K37



  if(argc > 1) {  
    if(strcmp(argv[1],"-h")==0) {
      cout << "First paramter is isotope (K37 or K38) [K38]" << endl;
      cout << "Second parameter is OP transition: 1 for D1 line 3 for D2 line [1]" << endl;
      cout << "Third parameter is tmax in ns [1]" << endl;
      cout << "Fourth paramter is B_ext in G [2]" << endl;
      cout << "Fifth paramter is dt in ns [1]" << endl;
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
  if( strcmp(isotope.c_str(),"K37")  == 0 || strcmp(isotope.c_str(),"37K")  ==0 ) { I2 = 3; Aj[0] = 0.1201336 ; Aj[1] = 0.01445  ; g_I = .39094 ; if(Je2 == 3) {omega_Di = 2*M_PI*(speed_of_light / (766.7 *pow(10,-7)));}} //see Besch 1968 and Dan's thesis
  if( strcmp(isotope.c_str(),"K38")  == 0 || strcmp(isotope.c_str(),"38K")  ==0 ) { I2 = 6; Aj[0] = 0.70765   ; Aj[1] = 0.0853   ; g_I = 0.2029 ; if(Je2 == 3) {omega_Di = 2*M_PI*(speed_of_light / (766.7 *pow(10,-7)));}} //see Besch 1968 and Dan's thesis
  if( strcmp(isotope.c_str(),"K47")  == 0 || strcmp(isotope.c_str(),"47K")  ==0 ) { I2 = 1; Aj[0] = 0.70765   ; Aj[1] = 0.0853   ; g_I = 0.2029 ; if(Je2 == 3) {omega_Di = 2*M_PI*(speed_of_light / (766.7 *pow(10,-7)));}} //see Besch 1968 and Dan's thesis THIS IS NOT REAL NUMBRES, JUST NEEDED AN ISOTOPE WITH I = 1/2 FOR TESTING!
  
  
  //The idea is to adjust these parameters even though they will later be switched to better units
  //*********Inputs, constants*******
  double laser_I_fe = 0.200; //mW/cm^2
  double laser_I_ge = 0.145; //mW/cm^2
  double laser_pol[3] = {0.0,0.0,1.0};
  double detune = 1.0; //MHz
  double gamma_1 = 0.2; //Mhz = Linewidth of laser on the g->e transition omega_1 in T&J
  double gamma_2 = 0.2; //Mhz = Linewidth of laser on the f->e transition omega_2 in T&J
  //*********************************
  
  
  
  double I = (double)I2/2.0;
  detune *= pow(10,-3); //puts detune in ns^-1
  double omega_2 = omega_Di - (0.5*Aj[0]* ( ((I+0.5)*(I+1.5)) - ( I * (I+1.0)) - 0.75)) - detune;  //Frequency of laser tuned to the f->e transition (ns^-1)
  double omega_1 = omega_Di  - (0.5*Aj[0]* ( ((I-0.5)*(I+0.5)) - ( I * (I+1.0)) - 0.75)) + detune;  //Frequency of laser tuned to the g->e transition (ns^-1)
  
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
  
  
  if(!production) {
    cout << "************** OPTICAL PUMPING WITH OPTICAL-BLOCH EQUATIONS **************" << endl << endl;
    cout << "Pumping " << isotope << "\t As = " << Aj[0] << " ns^1 \t Ap = " << Aj[1] << " ns^-1 \t g-factor: " << g_I << endl;
    if(I2%2 == 0) { cout << endl << "Nuclear spin: " << I2/2 << endl; } 
    else{           cout << endl << "Nuclear spin: " << I2 << "/2" << endl;}
    cout << "tau: " <<  tau << " ns \t gamma = " << gamma << "n s^-1" << endl;
    cout << "Tmax: " << tmax << " ns.  Time step: " << dt << " ns." << endl;
    
    cout << "Magnetic field: " << B_z << " G" << endl;
    
    printf("Laser 1 frequency %12.10G ns^-1 \t linewidth = %4.2G ns^-1 \t Intensity = %G g/ns^3 \t E-field = %10.8G g*cm/(A*ns^3) \n", omega_1, gamma_1, laser_I_ge, laser_E_ge);
    printf("Laser 2 frequency %12.10G ns^-1 \t linewidth = %4.2G ns^-1 \t Intensity = %G g/ns^3 \t E-field = %10.8G g*cm/(A*ns^3) \n", omega_2, gamma_2, laser_I_fe, laser_E_fe);
    //cout << "Epsilon_0 [MKS] = " << GSL_CONST_MKSA_VACUUM_PERMITTIVITY << endl; //A^2 s^4/(m^3 kg)
    //cout << "Speed of light [CGS] = " << GSL_CONST_CGSM_SPEED_OF_LIGHT << endl; //cm/s
    //cout << "Hbar [cgs] = " << GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR << endl; //g cm^2 /s
  }  
  
  
  //e States are any excited state either with J = 1/2, or 3/2 (D1 or D2)
  //f States are ground states with F = I + J
  //g States are ground states with F = I - J
  const int numEStates = (I2+1)*(Je2+1);
  const int numFStates = I2 + 2;
  const int numGStates = I2;

  //Stuff to get the I,J decomposition
  double groundDecomp[numFStates+numGStates][numFStates+numGStates];
  decomp.diagH(I2, 0, 1 , mu_B, B_z, 0.0, g_I, Aj[0], &groundDecomp[0][0]);
  double excitedDecomp[numEStates][numEStates];
  decomp.diagH(I2, 1, Je2,mu_B, B_z, 0.0, g_I, Aj[1], &excitedDecomp[0][0]);

  bool debugDecomp = true;
  if(!production && debugDecomp) {
    printf("\n\n IN OPTICALBLOCH:\n");
    for(int i = 0; i < numFStates+numGStates; i++) {
      for(int j = 0; j < numFStates+numGStates; j++) {
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

  //Getting an ordered list of Iz (same order that the helper class uses)
  double Izground[numFStates+numGStates];
  double Izexcited[numFStates+numGStates];
  for(int i = 0; i < numFStates + numGStates; i++)  Izground[i]  = -((double)I2/2.0) + (double)(i%(I2+1));
  for(int i = 0; i < numEStates             ; i++)  Izexcited[i] = -((double)I2/2.0) + (double)(i%(I2+1)); 



  const int maxTOut = 1000; //Maxiumum number of time steps to write to a file
  int printFreq; //Number of ns between print statements
  int totalPrint; //Total number of print statements
  if(tmax <= maxTOut) {
    printFreq = tmax;
    totalPrint = tmax;
  }
  else {
    printFreq = tmax / maxTOut;
    totalPrint = maxTOut;
  }

  FILE * file;
  if(production) {
    file = fopen("opData.dat","w");
    stdout = file;
    fprintf(file,"%i \t %6.2G \t %i 4.0\n",numEStates+numFStates+numGStates,8.0,totalPrint);
  }
  else {
    file = stdout;
  }
  
  if(!production) {
    cout << "Number of excited states: " << numEStates << endl;
    cout << "f States with 2F = " << I2 + 1 << " -->  " << numFStates << endl;
    cout << "g States with 2F = " << I2 - 1 << " -->  " << numGStates << endl;
  }

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
	if( i == j) {
	  fill = 1.0;
	}
	rhoee[i][j] = gsl_complex_rect(0.0,0.0);
	drhoee[i][j] = gsl_complex_rect(0.0,0.0);

	if(i < numFStates && j < numFStates) {
	  rhoff[i][j] = gsl_complex_rect(fill,0.0);
	  drhoff[i][j] = gsl_complex_rect(0.0,0.0);
	}
	if(j < numFStates) {
	  zrhoef[i][j] = gsl_complex_rect(0.0,0.0);
	  dzrhoef[i][j] = gsl_complex_rect(0.0,0.0);
	}
	if(i < numGStates && j < numGStates) {
	  rhogg[i][j] = gsl_complex_rect(fill,0.0);
	  drhogg[i][j] = gsl_complex_rect(0.0,0.0);
	} 
	if(j < numGStates) {
	  zrhoeg[i][j] = gsl_complex_rect(0.0,0.0);
	  dzrhoeg[i][j] = gsl_complex_rect(0.0,0.0);
	}	
	if(i < numFStates && j < numGStates) {
	  zrhofg[i][j] = gsl_complex_rect(0.0,0.0);
	  dzrhofg[i][j] = gsl_complex_rect(0.0,0.0);
	}


      }
   }
  
  //Array to hold the relative energy of each state E = 0 is at the S1/2 level before the hyperfine splitting!!!
  double omegaE[numEStates]; //Ordering of states: F = I - J --> F = I + J ;;; Mf = -F --> Mf = F
  double omegaF[numFStates]; 
  double omegaG[numGStates];
  
  int Fe2 = abs(I2 - Je2);
  int MFe2 = -Fe2;
  int Ff2 = (I2 + 1);
  int MFf2 = -Ff2;
  int Fg2 = (I2 - 1);
  int MFg2 = -Fg2;
    
  //Get all energies for all states
  bool debugEnergy = false;
  for(int i = 0; i < numEStates; i++) { 
    
    double I = ((double)I2)/2.0;
    double J_excited = ((double)Je2/2.0);
    double F_excited = ((double)Fe2)/2.0;
    
    //Excited states
    omegaE[i] = omega_Di;//Fine splitting [ns^-1]
    omegaE[i] += 0.5*Aj[1]*( (F_excited * (F_excited+1.0)) - (I * (I+1.0)) - (J_excited * (J_excited+1.0))); //Hyperfine [ns^-1]
    double g_Fe = decomp.calc_gf(Fe2, Je2, I2, 2, 1, g_I);
    omegaE[i] += g_Fe * mu_B * B_z * (double)MFe2 / 2.0; //[ns-1]
    if(!production && debugEnergy) printf(" F = %3.0f \t Mf = %d/2 \t omegaE = %f12.8 ns^-1 \n",F_excited,MFe2,omegaE[i]);
    
    //F states
    if( i < numFStates) {

      double F_f = (double)Ff2/2.0;
      
      omegaF[i] = 0.5*Aj[0]*( (F_f * (F_f + 1.0)) - (I * (I+1.0)) - 0.75); //[ns^-1]
      double g_Ff = decomp.calc_gf(Ff2, 1, I2, 0, 1,g_I);
      omegaF[i]+= g_Ff * mu_B * B_z * (double)MFf2 / 2.0; // ns^-1

      
      MFf2 += 2;
    }
    //G states
    if( i < numGStates) {
      double F_g = (double)Fg2/2.0;
      omegaG[i] = 0.5*Aj[0]*( (F_g * (F_g + 1.0)) - (I * (I+1.0)) - 0.75); //ns^-1
      double g_Fg = decomp.calc_gf(Fg2, 1, I2, 0, 1, g_I);
      omegaG[i] += g_Fg * mu_B * B_z * (double)MFg2 /2.0;  //ns^-1

      MFg2 +=2;
    }
    
    MFe2 += 2;
    
    
    if(MFe2 > Fe2) {
      Fe2 += 2;
      MFe2 = -Fe2;
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
  
  Fe2 = abs(I2 - Je2);
  MFe2 = -Fe2;
  Ff2 = (I2 + 1);
  MFf2 = -Ff2;
  Fg2 = (I2 - 1);
  MFg2 = -Fg2;
  
  bool debugCoupling = false;
  for(int i = 0; i < numEStates; i++) {  
    for(int j = 0; j < numFStates; j++) {
      for(int q = 0; q < 3; q++) {
	double I = ((double)I2) / 2.0;
	double J_e = ((double)Je2) / 2.0;
	double F_e = ((double)Fe2) / 2.0;
	double F_f = ((double)Ff2) / 2.0;
	double MFe = ((double)MFe2) / 2.0;
	aef[i][j][q] = pow(-1.0, 1 + I + J_e + F_e + F_f - MFe);
	aef[i][j][q] *= (sqrt( (2*F_f) + 1.0) * sqrt( (2*F_e) + 1.0) * sqrt( (2*J_e) + 1.0));
	if(!production && debugCoupling) cout << "F_f = " << F_f << "\t F_e = " << F_e << "\t J_e = " << J_e << endl;
	aef[i][j][q] *= gsl_sf_coupling_3j( Fe2 , 2.0 , Ff2 , -MFe2, (q-1)*2 , MFf2);
	aef[i][j][q] *= gsl_sf_coupling_6j( Fe2 , 2.0 , Ff2 , 1    , I2      , Je2 );
	
	double omega_ef = omegaE[i] - omegaF[j]; //ns^-1
	if(!production && debugCoupling) printf("omega_ef = %14.10G ns^-1 \n",omega_ef);
	Def[i][j][q] = aef[i][j][q] * sqrt(emitProb/pow(omega_ef,3)); // ns*cm*A
	if(!production && debugCoupling) printf("i = %i \t j = %i \t q = %i \t F_e = %2.1f \t Mfe = %2.1f \t F_f = %2.1f \t MF_f = %2.1f \t q = %2i \t a = %6.4f \t D = %10.8G ns*cm*A\n",i,j,q,F_e,MFe,(double)Ff2/2.0,(double)MFf2/2.0,q-1,aef[i][j][q],Def[i][j][q]);
	//	if(!production && debugCoupling) printf("aef[%i][%i][%i] = %8.6f \n",i,j,q,aef[i][j][q]);
	//if(!production && debugCoupling) printf("Def[%i][%i][%i] = %8.6G \n",i,j,q,Def[i][j][q]);	
	

	if(j < numGStates) {
	  double F_g = ((double)Fg2) / 2.0;
	  double MFg = ((double)MFe2) /2.0;
	  aeg[i][j][q] = pow(-1.0, I + J_e + F_e + F_g - MFg);
	  aeg[i][j][q] *= (sqrt( (2*F_g) + 1.0) * sqrt( (2*F_e) + 1.0) * sqrt( (2*J_e) + 1.0));
	  aeg[i][j][q] *= gsl_sf_coupling_3j( Fe2 , 2.0 , Fg2 , -MFe2, (q-1)*2 , MFg2);
	  aeg[i][j][q] *= gsl_sf_coupling_6j( Fe2 , 2.0 , Fg2 , 1    , I2      , Je2 );
	  
	  double omega_eg = omegaE[i] - omegaG[j]; //ns^-1
	  Deg[i][j][q] = aeg[i][j][q] * sqrt(emitProb/pow(omega_eg,3)); //ns cm A
	  //	    cout << "F_e = " << F_e << "\t MFe = " << MFe << "\t F_g = " << (double)Fg2/2 << "\t MF_g = " << (double)MFg2/2.0 << "\t q = " << q-1 << "\t a = " << aeg[i][j][q] << endl;
	}
      }
      MFf2 += 2;
      MFg2 += 2;
      
    }
    
    MFe2 += 2;
    if(MFe2 > Fe2) {
      Fe2 += 2;
      MFe2 = -Fe2;
    }
    //And reset the ground states
    Ff2 = I2 + 1;
    MFf2 = -Ff2;
    
    Fg2 = I2 - 1;
    MFg2 = -Fg2;
  }
  

  bool zCoherences = true;
  //bool oCoherences = false;
  //bool gsHCoherences = false;
  //bool esHCoherences = false;
  
  double totalPopulation = 0.0;
  for(int i = 0; i < numEStates; i++) {
    totalPopulation                    += GSL_REAL(rhoee[i][i]);
    if(i < numFStates) totalPopulation += GSL_REAL(rhoff[i][i]);
    if(i < numGStates) totalPopulation += GSL_REAL(rhogg[i][i]);
  }
  
  for(double time = 0; time < tmax; time += dt) {  

    //Print stuff out
    bool printMatrices = false;

    if(printMatrices  && !production ) {
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
      fprintf(file,"___________________________________________________________________________________________________________________________________\n");
      fprintf(file,"\n"); 
      
    } else if( (int)time%printFreq == 0) {
      fprintf(file,"%12.6G ",time);     
      for(int g = 0; g < numGStates; g++) {
	fprintf(file,"%10.5G ",GSL_REAL(rhogg[g][g]));
      }
      for(int f = 0; f < numFStates; f++) {
	fprintf(file,"%10.5G ",GSL_REAL(rhoff[f][f]));
      }
      for(int e = 0; e < numEStates; e++) {
	fprintf(file,"%10.5G ",GSL_REAL(rhoee[e][e]));
      }
      
      fprintf(file,"%10.5G ",totalPopulation);

      //Get the nuclear polarization
      bool debugPolarization = true;
      double nucPol = 0.0;
      for(int i = 0; i < (numGStates + numFStates); i++) {	
	if(debugPolarization && !production) printf("i = %i\n",i);
	for( int fzindex = 0; fzindex < (numFStates + numGStates); fzindex++) {
	  if( i < numGStates) {
	    nucPol += Izground[fzindex] * GSL_REAL(rhogg[i][i]) * groundDecomp[fzindex][i] * groundDecomp[fzindex][i];
	    if(debugPolarization && !production) printf("\tIz = %1.3G, Decomp = %8.6G, Population = %8.6G, Running Sum = %8.6G\n",Izground[fzindex],groundDecomp[fzindex][i], GSL_REAL(rhogg[i][i]), nucPol);
	  } else      {
	    nucPol += Izground[fzindex] * GSL_REAL(rhoff[i-numGStates][i-numGStates]) * groundDecomp[fzindex][i] * groundDecomp[fzindex][i];
	    if(debugPolarization && !production) printf("\tIz = %1.3G, Decomp = %8.6G, Population = %8.6G, Running Sum = %8.6G\n",Izground[fzindex],groundDecomp[fzindex][i], GSL_REAL(rhoff[i][i]), nucPol);
	  }

	}
      }
      for(int e = 0; e < numEStates; e++) {
	for( int fzindex = 0; fzindex < numEStates; fzindex++) {
	  nucPol += Izexcited[fzindex] * GSL_REAL(rhoee[e][e]) * excitedDecomp[fzindex][e] * excitedDecomp[fzindex][e];
	}
      }
      nucPol /= (double)I2/2.0;
      nucPol /= totalPopulation;
      fprintf(file,"%10.5G ",nucPol);


      fprintf(file,"\n"); 
    }

    

    bool debugGLaser_ee = false;
    bool debugFLaser_ee = false;
    //Excited state populations (Eq 32)
    for(int e = 0; e < numEStates; e++)  { 
      for(int ep = 0; ep < numEStates; ep++) {

	if( (!production && debugGLaser_ee) || (!production && debugFLaser_ee)) fprintf(file,"e = %i   ep = %i \n",e,ep);
	drhoee[e][ep] = gsl_complex_rect(-gamma, 0.0);
	if(zCoherences) { 
	  drhoee[e][ep] = gsl_complex_add( drhoee[e][ep], gsl_complex_rect( 0.0, -(omegaE[ep]-omegaE[e])));
	}
	drhoee[e][ep] = gsl_complex_mul(drhoee[e][ep], rhoee[e][ep]);


	//Now for the laser interaction!!!
	gsl_complex gLaserTerm_ee = gsl_complex_rect(0.0,0.0);
	gsl_complex fLaserTerm_ee = gsl_complex_rect(0.0,0.0);

	for(int iq = 0; iq < 3; iq ++) { //Dot product term
	  if(!production && debugGLaser_ee) fprintf(file,"\t iq = %i",iq);
	  for (int gpp = 0; gpp < numGStates; gpp++) { //sum over gpp in (T&J 32)
	    
	    gsl_complex temp = gsl_complex_mul( gsl_complex_rect( Deg[e][gpp][iq], 0.0), gsl_complex_conjugate(zrhoeg[ep][gpp]));
	    if(!production && debugGLaser_ee) fprintf(file,"\t left-temp = %12.9G \t + %12.9G i ", GSL_REAL(temp), GSL_IMAG(temp));
	    
	    gsl_complex temp2 = gsl_complex_mul( gsl_complex_rect( Deg[ep][gpp][iq], 0.0), zrhoeg[e][gpp]);
	    if(!production && debugGLaser_ee) fprintf(file,"\t right-temp = %12.9G \t + %12.9G i ", GSL_REAL(temp2), GSL_IMAG(temp2));

	    temp = gsl_complex_sub(temp,temp2);
	    temp = gsl_complex_mul_real(temp,laser_pol[iq]);
	    if(!production && debugGLaser_ee) fprintf(file,"\t after polarization = %12.9G \t + %12.9G i \n", GSL_REAL(temp), GSL_IMAG(temp));

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
	gLaserTerm_ee = gsl_complex_mul_imag( gLaserTerm_ee, 0.5/hbar);

	drhoee[e][ep] = gsl_complex_add(drhoee[e][ep],gLaserTerm_ee);
	drhoee[e][ep] = gsl_complex_mul_real(drhoee[e][ep],dt);
      }
    } //END EXCITED STATE POPULATIONS


    Ff2 = I2 + 1;
    MFf2 = -Ff2;
    //F-Ground state populations (Eq 33)
    bool debugSpontDecay_ff = false;
    for(int f = 0; f < numFStates; f++) {  
      
      int Ffp2 = I2 + 1;
      int MFfp2 = -Ffp2;

      for(int fp = 0; fp < numFStates; fp++) {
	
	drhoff[f][fp] = gsl_complex_rect(0.0,0.0);

	Fe2 = I2-Je2;
	int mFe2 = -Fe2;

	if(!production && debugSpontDecay_ff) fprintf(file,"f = %i   fp = %i   m_f = %i/2   m_fp = %i/2 \n",f,fp,MFf2,MFfp2);
	for(int e = 0; e < numEStates; e++) {

	  int Fep2 = I2-Je2;
	  int mFep2 = -Fep2;

	  for(int ep = 0; ep < numEStates; ep++) {
	    
	    if( (mFe2 - mFep2) == (MFf2 - MFfp2) ) {

	      
	      int qef = (mFe2 - MFf2)/2;
	      int qepfp = (mFep2 - MFfp2)/2;
	      if( abs(qef) <= 1 && abs(qepfp) <=1) { 

	
		if(f == fp) {
		  if(!production && debugSpontDecay_ff) fprintf(file,"\t e = %i   ep = %i     F_e = %i/2    m_Fe = %i/2    qef = %i    F_ep = %i/2     m_Fe = %i/2     qepfp = %i \n",e,ep,Fe2,mFe2,qef,Fep2,mFep2,qepfp);
		  if(!production && debugSpontDecay_ff) fprintf(file,"\t\t aef[%i][%i][%i] = %8.6f \t aef[%i][%i][%i] = %8.6f \n",e,f,qef+1,aef[e][f][qef+1],ep,fp,qepfp+1,aef[ep][fp][qepfp+1]);
		}
		gsl_complex spontDecay = gsl_complex_mul( gsl_complex_rect(aef[e][f][qef+1],0.0) , gsl_complex_rect(aef[ep][fp][qepfp+1],0.0));
		spontDecay             = gsl_complex_mul( spontDecay, rhoee[e][ep]);
		drhoff[f][fp] = gsl_complex_add(drhoff[f][fp],spontDecay);
		if(f == fp) {
		  if(!production && debugSpontDecay_ff) fprintf(file,"\t\t spontDecay = %6.4f + %14.10G i \t",GSL_REAL(spontDecay),GSL_IMAG(spontDecay));
		  if(!production && debugSpontDecay_ff) fprintf(file,"\t\t drhoff[%i][%i] = %6.4f + %6.4G i \n", f,fp,GSL_REAL(drhoff[f][fp]), GSL_IMAG(drhoff[f][fp]));
		}
		
	      }

	      if(zCoherences) {
		//Second half of the sum should take place regardless of whether or not q is allowed
		gsl_complex zCoherencesF = gsl_complex_rect(0.0, -(omegaF[fp]-omegaF[f]));
		zCoherencesF = gsl_complex_mul(zCoherencesF, rhoff[f][fp]);
		drhoff[f][fp] = gsl_complex_add( drhoff[f][fp], zCoherencesF);
	      }


	    }

	    mFep2 += 2;
	    if(mFep2 > Fep2) {
	      Fep2 += 2;
	      mFep2 = -Fep2;
	    }
	  } //end ep loop
	  mFe2 +=2;
	  if(mFe2 > Fe2) {
	    Fe2 += 2;
	    mFe2 = -Fe2;
	  }
	} //end e loop
	
	drhoff[f][fp] = gsl_complex_mul(drhoff[f][fp], gsl_complex_rect(gamma,0.0));
	MFfp2 += 2;
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
      
      MFfp2 = -Ffp2;
      MFf2 += 2;
    } //END F-GROUND STATE POPULATIONS


    Fg2 = I2 - 1;
    MFg2 = -Fg2;
    //G-Ground State populations (Eq 34)
    for(int g = 0; g < numGStates; g++) {
      int Fgp2 = I2 - 1;
      int MFgp2 = -Fgp2;

      for(int gp = 0; gp < numGStates; gp++) {

	drhogg[g][gp] = gsl_complex_rect(0.0,0.0);
	
	Fe2 = I2 - Je2;
	int mFe2 = -Fe2;

	for(int e = 0; e < numEStates; e++) {
	  
	  int Fep2 = I2 - Je2;
	  int mFep2 = -Fep2;
	  
	  for(int ep = 0; ep < numEStates; ep++) {
	    
	    if( (mFe2 - mFep2) == (MFg2 - MFgp2) ) {

	      int qeg = (mFe2 - MFg2)/2;
	      int qepgp = (mFep2 - MFgp2)/2;

	      if( abs(qeg) <= 1 && abs(qepgp) <=1) {
		 
		gsl_complex spontDecay = gsl_complex_mul( gsl_complex_rect(aeg[e][g][qeg+1],0.0), gsl_complex_rect(aeg[ep][gp][qepgp+1],0.0));
		spontDecay             = gsl_complex_mul( spontDecay, rhoee[e][ep]);
		drhogg[g][gp] = gsl_complex_add(drhogg[g][gp],spontDecay);


	      }

	      if(zCoherences) {
		gsl_complex zCoherencesG = gsl_complex_rect(0.0,-(omegaG[gp]-omegaG[g]));
		zCoherencesG = gsl_complex_mul(zCoherencesG, rhogg[g][gp]);
		drhogg[g][gp] = gsl_complex_add( drhogg[g][gp], zCoherencesG);
	      }
	    }

	    mFep2 +=2;
	    if(mFep2 > Fep2) {
	      Fep2 += 2;
	      mFep2 = -Fep2;
	    }
	  } //end ep loop

	  mFe2 +=2;
	  if(mFe2 > Fe2) {
	    Fe2 +=2;
	    mFe2 = -Fe2;
	  }
	} //end e loop

	drhogg[g][gp] = gsl_complex_mul(drhogg[g][gp],gsl_complex_rect(gamma,0.0));
	//END SPONTANEOUS EMISSION!!!

	//EG LASER PUMPING TERM
	gsl_complex laserTerm_gg = gsl_complex_rect(0.0,0.0);
	for(int iq = 0; iq < 3; iq++) {
	  for(int epp = 0; epp < numEStates; epp++) {
	    gsl_complex temp = gsl_complex_mul(  gsl_complex_rect( Deg[epp][g][iq] ,0.0), zrhoeg[epp][gp]);
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
	MFgp2 += 2;
      } //End Gp loop

      MFgp2 = -Fgp2;
      MFg2 +=2;
    } //END G-GROUND STATE POPULATIONS

    //EF-Coherences (Eq 35)
    bool debugDetune_ef = false;
    bool debugLaser2_ef = false;
    for(int e = 0; e < numEStates; e++) {
      for(int f = 0; f < numFStates; f++) {


	double real = (gamma+gamma_2) / 2.0;
	double imag = omegaE[e] - omegaF[f] - omega_2;
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
	  
	dzrhoef[e][f] = gsl_complex_mul_real(dzrhoef[e][f], dt);
      }
    } //END EF-COHERENCES

    //EG-Coherences (Eq 36)
    for(int e = 0; e < numEStates; e++) {
      for(int g = 0; g < numGStates; g++) {
	
	double real = (gamma +gamma_1) / 2.0;
	double imag = omegaE[e] - omegaG[g] - omega_1;
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
	
	dzrhoeg[e][g] = gsl_complex_mul_real( dzrhoeg[e][g]           , dt);
      }
    } //END EG-COHERENCES


    //Update with all the new populations
    for(int i = 0; i < numEStates; i++) {
      for(int j = 0; j < numEStates; j++) { 
	if(i < numEStates && j < numEStates)  rhoee[i][j] = gsl_complex_add( rhoee[i][j], drhoee[i][j]);
	if(i < numFStates && j < numFStates)  rhoff[i][j] = gsl_complex_add( rhoff[i][j], drhoff[i][j]);
	if(i < numGStates && j < numGStates)  rhogg[i][j] = gsl_complex_add( rhogg[i][j], drhogg[i][j]);
	if(i < numEStates && j < numFStates) zrhoef[i][j] = gsl_complex_add(zrhoef[i][j],dzrhoef[i][j]);
	if(i < numEStates && j < numGStates) zrhoeg[i][j] = gsl_complex_add(zrhoeg[i][j],dzrhoeg[i][j]);
	
      }
    }

    for(int i = 0; i < numEStates; i++) {
      
      if(                     GSL_REAL(rhoee[i][i]) < 0.0) GSL_SET_REAL(&rhoee[i][i],0.0);
      if( i < numFStates &&   GSL_REAL(rhoff[i][i]) < 0.0) GSL_SET_REAL(&rhoff[i][i],0.0);
      if( i < numGStates &&   GSL_REAL(rhogg[i][i]) < 0.0) GSL_SET_REAL(&rhogg[i][i],0.0);

    }
    
    if(!isHermitian(&rhoee[0][0],numEStates) || !isHermitian(&rhoff[0][0], numFStates) || !isHermitian(&rhogg[0][0], numGStates)) {
	fprintf(file,"DENSITY MATRIX NOT HERMITIAN!!!\n");
	fprintf(file,"time = %.2G",time);
	exit(1);
    }
      

  }  //END TIME LOOP
  
} 

void printArray2D_complex(gsl_complex *array, int col, int row) { 
  
  gsl_complex value;
  
  for( int r = 0; r < row; r++) {
    for(int c = 0; c < col; c++) {
      value = *(array+c+(col*r));
      //      cout << setw(10) << left << GSL_REAL(value);
      fprintf(stdout,"%12.6G    ",GSL_REAL(value));
    }

    fprintf(stdout,"\t\t");
    for(int c = 0; c < col; c++) {
      value = *(array+c+(col*r));
      fprintf(stdout,"%12.6G i    ", GSL_IMAG(value));
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
  bool eps = pow(10,-6);
  gsl_complex value;
  //First check that the diagonals are real
  for(int i = 0 ; i < size; i++) {
    value = *(array + i + (size*i));
    if( fabs(GSL_IMAG(value)) > eps) hermit = false;
  }

  gsl_complex upper;
  gsl_complex lower;
  
  for(int c = 0; c < size; c++) {
    for(int r = 0; r < size; r++) {
      upper = *(array+c+(size*r));
      lower = *(array+r+(size*c));
      
      if( fabs( GSL_REAL(upper) - GSL_REAL(lower) > eps )) hermit = false;
      if( fabs( GSL_IMAG(upper) + GSL_IMAG(lower) > eps )) hermit = false;

    }
  }
  return hermit;
}
