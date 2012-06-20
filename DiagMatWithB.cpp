#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "math.h"
//#include "genAtomicState.cpp"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;

class DiagMatWithB {

public:

  void diagH(int I2, int L, int J2, double mu_B, double B_z, double B_x, double g_I, double Aj, double *arrayToFill);
  double calc_gf(int F2, int J2, int I2,int L2, int S2,double g_I);
  double calc_gj(int J2, int L2, int S2);

};


void DiagMatWithB::diagH(int I2, int L, int J2, double mu_B, double B_z, double B_x, double g_I, double Aj, double *arrayToFill) {

  //bool debug = true;
  bool debug = false;
  printf( "Decomposing nuclear spin I = %i/2 for the L = %i ; J = %i/2 state.  Aj = %6.4G",I2,L,J2,Aj);
  int numBasisStates = (I2 + 1)*(J2+1);
  //genAtomicState::numBasisStates = numBasisStates;
  //  int Iz2[8] = {-3, -1, 1 , 3 , -3, -1, 1, 3}; these define the ordering of the states
  int Iz2[numBasisStates*2];
  int Jz2[numBasisStates*2];

  for(int i = 0; i < numBasisStates*2; i+=2) {
    Iz2[i] = -I2 + 2*((i/2)%(I2+1));
    Iz2[i+1] = 0.0; //The imaginary part
    
    Jz2[i] = -J2+ 2*(i/2/(I2+1));
    Jz2[i+1] = 0.0;
    
    //genAtomicState::Iz2[i] = -I2 + 2*i;
    //genAtomicState::Iz2[i+(numBasisStates/2)] = -I2+(2*i);
    //genAtomicState::Jz2[i] = -1;
    //genAtomicState::Jz2[i+(numBasisStates/2)] = 1;

  }
  
  cout << endl;

  for(int i = 0; i < numBasisStates*2; i+=2) {
    //cout << Iz2[i] << "/2 + " << Iz2[i+1] << " i , " << Jz2[i] << "/2 +" << Jz2[i+1] << " i " << endl;
  }
  //cout << endl;

  
  double I_z[2*numBasisStates*numBasisStates];
  double J_z[2*numBasisStates*numBasisStates];
  double I_plus[2*numBasisStates*numBasisStates];
  double J_plus[2*numBasisStates*numBasisStates];
  double I_minus[2*numBasisStates*numBasisStates];
  double J_minus[2*numBasisStates*numBasisStates];
  double H[2*numBasisStates*numBasisStates];

  if(debug) {cout << "numBasisStates: " << numBasisStates << endl;}
  for(int i = 0; i < (2*numBasisStates*numBasisStates)-1; i+=2) {
    H[i] = 0.0;

    if( i/2 % (numBasisStates+1) == 0) {
      I_z[i] = (double)Iz2[2*((i/2)%numBasisStates)]/2.0;
      J_z[i] = (double)Jz2[2*((i/2)%numBasisStates)]/2.0;


    } else {
      I_z[i] = 0.0;
      J_z[i] = 0.0;
    }
    I_z[i+1] = 0.0;  //Im
    J_z[i+1] = 0.0;
    //cout << "set I_z[" << i << "] = " << I_z[i] << " , I_z[" << i+1 << "] = " << I_z[i+1] << endl;
    
    
    if( ((i/2)+1) % (numBasisStates+1) == 0) { //One below the diagonal only
      double I = (double)I2 / 2.0;
      double Iz = (double)Iz2[2*((i/2)%(numBasisStates))]/2.0;
      I_plus[i] = sqrt( I*(I+1.0) - Iz*(Iz+1.0));
    } else {
      I_plus[i] = 0.0;
    }
    I_plus[i+1] = 0.0;

    

    if( ((i/2)-1) % (numBasisStates+1) == 0) { //One above the diagonal only
      double I = (double)I2 / 2.0; 
      double Iz = (double)Iz2[2*((i/2)%(numBasisStates))]/2.0;
      I_minus[i] = sqrt( I*(I+1) - Iz*(Iz-1.0));
    } else {
      I_minus[i] = 0.0;
    }
    I_minus[i+1] = 0.0;

    
    if( ((i/2)+(numBasisStates/(J2+1)))%(numBasisStates+1) ==0 ) { //Below the diagonal only
      double J = (double)J2 / 2.0;
      double Jz = (double)Jz2[(2*((i/2)%numBasisStates))]/2.0;
      J_plus[i] = sqrt( J*(J+1.0) - Jz*(Jz+1.0));
    } else {
      J_plus[i] = 0.0;
    }
    J_plus[i+1] = 0.0;
    
    
    if( ((i/2)-(numBasisStates/(J2+1)))%(numBasisStates+1) == 0) {//  && (i/2) < (numBasisStates*numBasisStates/2)) { //Above diagonal only
      double J = (double)J2 / 2.0;
      double Jz = (double)Jz2[2* ( ((((i/2)-(numBasisStates/(J2+1)))%numBasisStates)+(numBasisStates/(J2+1))) %numBasisStates)]/2.0;
      //      printf("Index I'm looking at = %i and Jz = %1.2f\n",((((i/2)-(numBasisStates/(J2+1)))%numBasisStates)+(numBasisStates/(J2+1))),Jz);
      if( sqrt(J*(J+1.0) - Jz*(Jz-1.0)) > 0.0 ) {
	J_minus[i] = sqrt(J*(J+1.0) - Jz*(Jz-1.0));
      } else {
	J_minus[i] = 0.0;
      }
    } else {
      J_minus[i] = 0.0;
    }
    J_minus[i+1] = 0.0;
    
    
  }
  
  if(debug) {
    for(int i = 0; i < (2*numBasisStates*numBasisStates) -1; i+=2) {
      
      cout <<  setw(12) << left << I_z[i] ;//<< " + " << I_minus[i+1] << "i , ";
      if( (i+2) %(2*numBasisStates) == 0) {
	cout << endl;
      }
      
    }
    cout << endl;
  }
  
  
  gsl_matrix_complex_view I_z_view     = gsl_matrix_complex_view_array (I_z    , numBasisStates, numBasisStates);
  gsl_matrix_complex_view J_z_view     = gsl_matrix_complex_view_array (J_z    , numBasisStates, numBasisStates);
  gsl_matrix_complex_view I_plus_view  = gsl_matrix_complex_view_array (I_plus , numBasisStates, numBasisStates);
  gsl_matrix_complex_view J_plus_view  = gsl_matrix_complex_view_array (J_plus , numBasisStates, numBasisStates);
  gsl_matrix_complex_view I_minus_view = gsl_matrix_complex_view_array (I_minus, numBasisStates, numBasisStates);
  gsl_matrix_complex_view J_minus_view = gsl_matrix_complex_view_array (J_minus, numBasisStates, numBasisStates);
  gsl_matrix_complex_view H_view       = gsl_matrix_complex_view_array (H      , numBasisStates, numBasisStates);
  //gsl_matrix_complex_view F_z_view     = gsl_matrix_complex_view_array (F_z    , numBasisStates, numBasisStates);

  

  gsl_matrix_complex *Hbz = gsl_matrix_complex_calloc(numBasisStates,numBasisStates); //Initializes to all zeroes
  gsl_matrix_complex *tempz = gsl_matrix_complex_calloc(numBasisStates,numBasisStates); //Initializes to all zeroes
  gsl_matrix_complex *Jx = gsl_matrix_complex_calloc(numBasisStates,numBasisStates);
  gsl_matrix_complex *Ix = gsl_matrix_complex_calloc(numBasisStates,numBasisStates);
  gsl_matrix_complex *Jy = gsl_matrix_complex_calloc(numBasisStates,numBasisStates);
  gsl_matrix_complex *Iy = gsl_matrix_complex_calloc(numBasisStates,numBasisStates);
  gsl_matrix_complex *Hbt = gsl_matrix_complex_calloc(numBasisStates,numBasisStates);
  gsl_matrix_complex *tempt = gsl_matrix_complex_calloc(numBasisStates,numBasisStates); //Initializes to all zeroes

  double memn = 1.0/1836.152701;
  gsl_complex g_Ip = gsl_complex_rect(-g_I * memn,0.0);
  
  gsl_matrix_complex_add  (Hbz, &I_z_view.matrix);
  cout << "Scaling by g_Ip = " << GSL_REAL(g_Ip) << endl;
  gsl_matrix_complex_scale(Hbz, g_Ip);

  gsl_matrix_complex_add (tempz, &J_z_view.matrix);
  cout << "Scaling by g_J = " << calc_gj(1.0,L*2,1.0) << endl << endl;
  gsl_matrix_complex_scale(tempz, gsl_complex_rect(calc_gj(1.0,L*2,1.0),0.0));
		   
  gsl_matrix_complex_add (Hbz, tempz);
  gsl_matrix_complex_scale (Hbz,gsl_complex_rect((mu_B*B_z),0.0));

  //Now for the transverse field components!
  gsl_matrix_complex_add(Jx, &J_plus_view.matrix);
  gsl_matrix_complex_add(Jx, &J_minus_view.matrix);
  gsl_matrix_complex_scale(Jx, gsl_complex_rect(0.5,0.0));
  
  gsl_matrix_complex_add(Ix, &I_plus_view.matrix);
  gsl_matrix_complex_add(Ix, &I_minus_view.matrix);
  gsl_matrix_complex_scale(Ix, gsl_complex_rect(0.5,0.0));

  gsl_matrix_complex_add(Jy, &J_minus_view.matrix);
  gsl_matrix_complex_sub(Jy, &J_plus_view.matrix);
  gsl_matrix_complex_scale(Jy, gsl_complex_rect(0.0,0.5));

  gsl_matrix_complex_add(Iy, &I_minus_view.matrix);
  gsl_matrix_complex_sub(Iy, &I_plus_view.matrix);
  gsl_matrix_complex_scale(Iy, gsl_complex_rect(0.0,0.5));

  //gsl_matrix_complex_fprintf(stdout, Iy, "%g");
  //gsl_matrix_complex_fprintf(stdout, Hb, "%g");

  //Now that I've got my I,J x,y components I need to build on 3.20 from Dan's thesis
  gsl_matrix_complex_add (Hbt, Ix);
  gsl_matrix_complex_scale(Hbt,g_Ip);
  
  gsl_matrix_complex_add(tempt, Jx);
  gsl_matrix_complex_scale(tempt,gsl_complex_rect(calc_gj(1.0,L*2,1.0),0.0));

  gsl_matrix_complex_add(Hbt,tempt);
  gsl_matrix_complex_scale(Hbt,gsl_complex_rect((mu_B*B_x),0.0));
  //gsl_matrix_complex_fprintf(stdout,Hbt,"%g");
  //*******************************************************************
  

  gsl_blas_zgemm( CblasNoTrans , CblasNoTrans , gsl_complex_rect(1.0,0.0), &I_z_view.matrix,     &J_z_view.matrix    , gsl_complex_rect(0.0,0.0) , &H_view.matrix );
  gsl_blas_zgemm( CblasNoTrans , CblasNoTrans , gsl_complex_rect(0.5,0.0), &I_plus_view.matrix,  &J_minus_view.matrix, gsl_complex_rect(1.0,0.0) , &H_view.matrix);
  gsl_blas_zgemm( CblasNoTrans , CblasNoTrans , gsl_complex_rect(0.5,0.0), &I_minus_view.matrix, &J_plus_view.matrix , gsl_complex_rect(1.0,0.0) , &H_view.matrix);
  gsl_matrix_complex_scale(&H_view.matrix, gsl_complex_rect(Aj,0.0));

  gsl_matrix_complex_add(&H_view.matrix, Hbz);
  gsl_matrix_complex_add(&H_view.matrix, Hbt);


  //  gsl_matrix_complex_fprintf(stdout, &J_z_view.matrix, "%g");
  //  printf("[");

  /*
  cout << endl << "H = " << endl;
  //cout << "[";
  for(int i = 0; i < (2*numBasisStates*numBasisStates)-1; i+=2) {
    cout << setw(16) << left << H[i];// << " + " << H[i+1] << " i , ";
    if( (i+2) % (2*numBasisStates) == 0) {
      cout << endl;
    }
  }
  cout << endl << endl << endl;
  */
  

  
  gsl_vector *eval = gsl_vector_alloc(numBasisStates);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(numBasisStates,numBasisStates);
  
  gsl_eigen_hermv_workspace * w = gsl_eigen_hermv_alloc(numBasisStates);

  //ALL THE MAGIC
  gsl_eigen_hermv(&H_view.matrix,eval,evec,w);
  gsl_eigen_hermv_free(w);
  gsl_eigen_hermv_sort(eval,evec, GSL_EIGEN_SORT_VAL_ASC);
  

  int F2[numBasisStates];
  int Fz2[numBasisStates];

  int F2current = abs(I2 - J2);
  int Fz2current = -F2current;

  for(int i = 0; i < numBasisStates; i++) {

    F2[i] = F2current;
    Fz2current += 2;
    if( Fz2current > F2current) {
      F2current += 2;
      Fz2current = -F2current;
    }

  }
  
  for(int i = 0; i < numBasisStates; i++) {
    
    //Now for the Fz value
    int setFz2;
    bool firstAdMix = true;


    //cout << "Doing i = " << i << endl;
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec,i);
      for(int k = 0; k < numBasisStates; k++) { 
	gsl_complex adMix = gsl_vector_complex_get(&evec_i.vector,k);
	if(fabs(GSL_REAL(adMix)) > pow(10,-10)) {
	  //printf("Mixing partly with state Iz2 = %i and Jz2 = %i\n",Iz2[2*k],Jz2[2*k]);
	  if(firstAdMix) { //This is the first time through
	    setFz2 = Iz2[2*k] + Jz2[2*k];
	    firstAdMix = false;
	  } else { //Check that it is consistent
	    //if(Iz2[k]+Jz2[k] != setFz2) { cout << "Houston, we have a problem." << endl; }
	  }
	}

      
      
      }

    Fz2[i] = setFz2;
    //genAtomicState::Fz2vec[L][i] = setFz2;
  }
  

    /*
  for(int i = 0; i < numBasisStates; i++) {
    double eval_i = gsl_vector_get(eval,i);
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec,i);
    cout << "eigenvalue =  " << eval_i ;//<< "\t State has 2F = " << F2[i] << "\t and 2Fz = " << Fz2[i] << endl;
    cout << "\t\t | " << F2[i] << "/2 , " << Fz2[i] << "/2 > = ";
    //     cout << "eigenvector = " << endl;
    for(int k = 0; k < numBasisStates; k++) {
      gsl_complex adMix = gsl_vector_complex_get(&evec_i.vector,k);
      if( fabs(GSL_REAL(adMix)) > pow(10,-10)) {
	cout << GSL_REAL(adMix) << " | " << Iz2[2*k] << "/2 , " << Jz2[2*k] << "/2 > + ";
      }
    }
    cout << endl << endl;
    //gsl_vector_fprintf(stdout,&evec_i.vector,"%g");
    
  }
  */
  int F2desired = abs(I2 - J2);
  int Fz2desired = -F2desired;
  
  double admixture[numBasisStates][numBasisStates];

  for(int i = 0; i < numBasisStates; i++) {
    //printf("Looking for F2 = %i and Fz2 = %i \n",F2desired, Fz2desired);
    for(int j = 0; j < numBasisStates; j++) {
      if( F2[j] == F2desired && Fz2[j] == Fz2desired) {
	//printf("FOUND IT!\n");
	gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec,j);
	for(int k = 0; k < numBasisStates; k++) {
	  gsl_complex adMix = gsl_vector_complex_get(&evec_i.vector,k);
	  admixture[k][i] = GSL_REAL(adMix);
	}
      }
    }

    Fz2desired += 2;
    if( Fz2desired > F2desired ) {
      F2desired+=2;
      Fz2desired = -F2desired;

    }

  }
  

  for(int i = 0; i < numBasisStates; i++) {
    for(int j = 0; j < numBasisStates; j++) {
      
      *(arrayToFill+j+(i*numBasisStates)) = admixture[i][j];

      //printf("%10.6G\t ",admixture[i][j]);
    }
    //printf("\n");
  }
  
  /*  
  for(int i = 0; i < numBasisStates; i++) {
    double eval_i = gsl_vector_get(eval,i);
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec,i);
    printf("eigenvalue = %g\n", eval_i);
    printf("eigenvector = \n");
    gsl_vector_complex_fprintf(stdout, &evec_i.vector,"%g");
    printf("\n\n");

  }
  */
 
  /*
  for(int i = 0; i < numBasisStates;i++) {

    double eval_i = gsl_vector_get(eval,i);
    cout << setw(14) << left << 1000.0*eval_i;
  }
  
  cout << endl << endl;
  
  for(int i = 0; i < numBasisStates; i++) {

    for(int j = 0; j < numBasisStates; j++) {
      gsl_vector_complex_view evec_j = gsl_matrix_complex_column(evec,j);
      gsl_complex ev = gsl_vector_complex_get(&evec_j.vector,i);
      cout << setw(14) << left << GSL_REAL(ev); //<< " + " << setw(12) << left << GSL_IMAG(ev) << "i";
    }
    cout << endl;
  }
  cout << endl;
  */
  //genAtomicState::eval[L] = eval;
  //genAtomicState::evec[L] = evec;
  //  gsl_vector_free(eval);
  //gsl_matrix_complex_free(evec);
  

  
}

double DiagMatWithB::calc_gj(int J2, int L2, int S2) {

  if(J2 == 0) {
    return 0.0;
  } else {
    double g_L = 1.0; //orbital g-factor
    double g_s = 2.0023; //electron spin g-factor
    
    double J = (double)J2/2.0;
    double L = (double)L2/2.0;
    double S = (double)S2/2.0;
    double num = ((J*(J+1.0)) + (L*(L+1.0)) - (S*(S+1.0)))*g_L;
    double den = 2*J*(J+1.0);
    
    double g_J = num/den;
    
    num = ((J*(J+1.0)) - (L*(L+1.0)) + (S*(S+1.0)))*g_s;
    g_J += (num/den);
    //cout << "For L = " << L << " g_J = " << g_J << endl;
    return g_J;
  }

}
 
double DiagMatWithB::calc_gf(int F2, int J2, int I2,int L2, int S2,double g_I) {

  if(F2 == 0) {
    return 0;
  } else {

    double g_J = calc_gj(J2, L2, S2);
    double memn = 1.0/1836.152701; //ratio of electron and nucleon masses
    double g_Ip = -g_I*memn;
    //cout << "g_Ip = " << g_Ip << endl;
    
    double F = (double)F2/2.0;
    double J = (double)J2/2.0;
    double I = (double)I2/2.0;
    
    double num = ((F*(F+1.0)) + (J*(J+1.0)) - (I*(I+1.0)))*g_J;
    double den = 2*F*(F+1.0);
  
    double g_F = num/den;
    
    num = ((F*(F+1.0)) - (J*(J+1.0)) + (I*(I+1.0)))*g_Ip;
    //g_F += (num/den);

    //cout << "g_F = " << g_F << endl;
    return g_F;
  }
}
