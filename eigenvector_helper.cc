// Authors: Benjamin Fenker 2013
// Copyright 2012 Benjamin Fenker


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <iomanip>
#include <vector>
#include "/usr/include/math.h"
#include "include/eigenvector_helper.h"
#include "include/units.h"

using std::string;
extern bool op_verbose;
extern bool op_batch;

Eigenvector_Helper::Eigenvector_Helper(atom_data set_atom,
                                       magnetic_field_data set_field)
    :atom(set_atom), field(set_field),
     nuclear_spin_ground(2*atom.numBasisStates_ground, 0),
     nuclear_spin_excited(2*atom.numBasisStates_excited, 0),
     total_atomic_spin_ground(2*atom.numBasisStates_ground, 0),
     total_atomic_spin_excited(2*atom.numBasisStates_excited, 0),
     IzJz_decomp_ground(atom.numBasisStates_ground,
                        vector<double>(atom.numBasisStates_ground, 0.0)),
     IzJz_decomp_excited(atom.numBasisStates_excited,
                         vector<double>(atom.numBasisStates_excited, 0.0)),
     eShift_excited(atom.numEStates, 0.0),
     eShift_ground(atom.numFStates+atom.numGStates, 0.0) {
  if (op_verbose) printf("Eigenvector_Helper::Eigenvector_Helper(...)\n");
  nuclear_spin_ground = get_nuclear_spin(atom.numBasisStates_ground);
  nuclear_spin_excited = get_nuclear_spin(atom.numBasisStates_excited);
  total_atomic_spin_ground = get_total_atomic_spin(atom.numBasisStates_ground,
                                                   1);
  total_atomic_spin_excited = get_total_atomic_spin(atom.numBasisStates_excited,
                                                    atom.Je2);
  IzJz_decomp_ground = diagH(0);
  IzJz_decomp_excited = diagH(1);
  /*
  printf("Excited decomp:\n");
  for (int i = 0; i < atom.numBasisStates_excited; i++) {
    for (int j = 0; j < atom.numBasisStates_excited; j++) {
      printf("%8.6G   ", IzJz_decomp_excited[i][j]);
    }
    printf("\n");
  }
  */
}

vector<vector<double> > Eigenvector_Helper::diagH(int L) {
  bool debug = false;
  // Figure out some parameters to use based on J
  int J2;
  double Aj;
  vector<int> Iz2, Jz2;
  if (L == 0) {
    J2 = 1;
    Aj = atom.Aj_g;
    Iz2 = nuclear_spin_ground;
    Jz2 = total_atomic_spin_ground;
  } else if (L == 1) {
    J2  = atom.Je2;
    Aj = atom.Aj_e;
    Iz2 = nuclear_spin_excited;
    Jz2 = total_atomic_spin_excited;
  } else {
    printf("FATAL ERROR.  L MUST EQUAL 0 OR 1\n L = %d\n", L);
    exit(1);
  }
  if (op_verbose) {
    printf("Decomposing nuclear spin I = %i/2 for the L = %i ; ", atom.I2, L);
    printf("J = %i/2 state.\n\tAj = %6.4G MHz\n", J2, Aj/_MHz);
    printf("\tB = %6.4G z + %6.4G x G\n\n", field.B_z/_G, field.B_x/_G);
  }
  int numBasisStates = (atom.I2 + 1)*(J2+1);
  // genAtomicState::numBasisStates = numBasisStates;

  if (op_verbose) printf("\n");

  double *I_z = new double[2*numBasisStates*numBasisStates];
  double *J_z = new double[2*numBasisStates*numBasisStates];
  double *I_plus = new double[2*numBasisStates*numBasisStates];
  double *J_plus = new double[2*numBasisStates*numBasisStates];
  double *I_minus = new double[2*numBasisStates*numBasisStates];
  double *J_minus = new double[2*numBasisStates*numBasisStates];
  double *H = new double[2*numBasisStates*numBasisStates];

  if (debug) printf("numBasisStates: %d\n", numBasisStates);
  for (int i = 0; i < (2*numBasisStates*numBasisStates)-1; i+=2) {
    H[i] = 0.0;

    if (i/2 % (numBasisStates+1) == 0) {
      I_z[i] = static_cast<double>(Iz2[
                                       2*((i/2)%numBasisStates)])/2.0;
      J_z[i] = static_cast<double>(Jz2[
                                       2*((i/2)%numBasisStates)])/2.0;
    } else {
      I_z[i] = 0.0;
      J_z[i] = 0.0;
    }
    I_z[i+1] = 0.0;  // Im
    J_z[i+1] = 0.0;

    if (((i/2)+1) % (numBasisStates+1) == 0) {  // One below the diagonal only
      double I = static_cast<double>(atom.I2) / 2.0;
      double Iz = static_cast<double>(Iz2[
                                          2*((i/2)%(numBasisStates))])/2.0;
      I_plus[i] = sqrt(I*(I+1.0) - Iz*(Iz+1.0));
    } else {
      I_plus[i] = 0.0;
    }
    I_plus[i+1] = 0.0;

    if (((i/2)-1) % (numBasisStates+1) == 0) {  // One above the diagonal only
      double I = static_cast<double>(atom.I2) / 2.0;
      double Iz = static_cast<double>(Iz2[
                                          2*((i/2)%(numBasisStates))])/2.0;
      I_minus[i] = sqrt(I*(I+1) - Iz*(Iz-1.0));
    } else {
      I_minus[i] = 0.0;
    }
    I_minus[i+1] = 0.0;

    if (((i/2)+(numBasisStates/(J2+1)))%(numBasisStates+1) == 0) {
      //  Below the diagonal only
      double J = static_cast<double>(J2) / 2.0;
      double Jz = static_cast<double>(Jz2[
                                          (2*((i/2)%numBasisStates))])/2.0;
      J_plus[i] = sqrt(J*(J+1.0) - Jz*(Jz+1.0));
    } else {
      J_plus[i] = 0.0;
    }
    J_plus[i+1] = 0.0;

    if (((i/2)-(numBasisStates/(J2+1)))%(numBasisStates+1) == 0) {
      //  && (i/2) < (numBasisStates*numBasisStates/2)) { //Above diagonal only
      double J = static_cast<double>(J2) / 2.0;
      int index = 2* ( ((((i/2)-(numBasisStates/(J2+1)))%numBasisStates)+
                        (numBasisStates/(J2+1))) %numBasisStates);
      double Jz = static_cast<double>(Jz2[index])/2.0;
      if (sqrt(J*(J+1.0) - Jz*(Jz-1.0)) > 0.0) {
        J_minus[i] = sqrt(J*(J+1.0) - Jz*(Jz-1.0));
      } else {
        J_minus[i] = 0.0;
      }
    } else {
      J_minus[i] = 0.0;
    }
    J_minus[i+1] = 0.0;
  }
  if (debug) {
    for (int i = 0; i < (2*numBasisStates*numBasisStates) -1; i+=2) {
      printf("%6.4G", J_plus[i]);
      if ((i+2) %(2*numBasisStates) == 0) {
        printf("\n");
      }
    }
    printf("\n");
  }

  gsl_matrix_complex_view I_z_view = gsl_matrix_complex_view_array(I_z,
                                         numBasisStates, numBasisStates);
  gsl_matrix_complex_view J_z_view = gsl_matrix_complex_view_array(J_z,
                                         numBasisStates, numBasisStates);
  gsl_matrix_complex_view I_plus_view = gsl_matrix_complex_view_array(I_plus,
                                            numBasisStates, numBasisStates);
  gsl_matrix_complex_view J_plus_view = gsl_matrix_complex_view_array(J_plus,
                                            numBasisStates, numBasisStates);
  gsl_matrix_complex_view I_minus_view = gsl_matrix_complex_view_array(I_minus,
                                             numBasisStates, numBasisStates);
  gsl_matrix_complex_view J_minus_view = gsl_matrix_complex_view_array(J_minus,
                                             numBasisStates, numBasisStates);
  gsl_matrix_complex_view H_view = gsl_matrix_complex_view_array(H,
                                       numBasisStates, numBasisStates);

  // Initializes to all zeroes
  gsl_matrix_complex *Hbz = gsl_matrix_complex_calloc(numBasisStates,
                                                      numBasisStates);
  gsl_matrix_complex *tempz = gsl_matrix_complex_calloc(numBasisStates,
                                                        numBasisStates);
  gsl_matrix_complex *Jx = gsl_matrix_complex_calloc(numBasisStates,
                                                     numBasisStates);
  gsl_matrix_complex *Ix = gsl_matrix_complex_calloc(numBasisStates,
                                                     numBasisStates);
  // gsl_matrix_complex *Jy = gsl_matrix_complex_calloc(numBasisStates,
  //                                                    numBasisStates);
  // gsl_matrix_complex *Iy = gsl_matrix_complex_calloc(numBasisStates,
  //                                                    numBasisStates);
  gsl_matrix_complex *Hbt = gsl_matrix_complex_calloc(numBasisStates,
                                                      numBasisStates);
  gsl_matrix_complex *tempt = gsl_matrix_complex_calloc(numBasisStates,
                                                        numBasisStates);
  const double memn = 1.0/1836.152701;
  gsl_complex g_Ip = gsl_complex_rect(-atom.g_I * memn, 0.0);

  gsl_matrix_complex_add(Hbz, &I_z_view.matrix);
  if (debug) printf("Scaling Iz by g_Ip = %8.6G\n", GSL_REAL(g_Ip));
  gsl_matrix_complex_scale(Hbz, g_Ip);

  gsl_matrix_complex_add(tempz, &J_z_view.matrix);
  if (debug) printf("Scaling Jz by g_J = %8.6G\n\n", calc_gj(1, L*2, 1));
  gsl_matrix_complex_scale(tempz,
                           gsl_complex_rect(calc_gj(1, L*2, 1), 0.0));

  gsl_matrix_complex_add(Hbz, tempz);
  gsl_matrix_complex_scale(Hbz,
                   gsl_complex_rect((_bohr_magneton/_planck_h*field.B_z), 0.0));

  // Now for the transverse field components!
  gsl_matrix_complex_add(Jx, &J_plus_view.matrix);
  gsl_matrix_complex_add(Jx, &J_minus_view.matrix);
  gsl_matrix_complex_scale(Jx, gsl_complex_rect(0.5, 0.0));

  gsl_matrix_complex_add(Ix, &I_plus_view.matrix);
  gsl_matrix_complex_add(Ix, &I_minus_view.matrix);
  gsl_matrix_complex_scale(Ix, gsl_complex_rect(0.5, 0.0));

  // gsl_matrix_complex_add(Jy, &J_minus_view.matrix);
  // gsl_matrix_complex_sub(Jy, &J_plus_view.matrix);
  // gsl_matrix_complex_scale(Jy, gsl_complex_rect(0.0, 0.5));

  // gsl_matrix_complex_add(Iy, &I_minus_view.matrix);
  // gsl_matrix_complex_sub(Iy, &I_plus_view.matrix);
  // gsl_matrix_complex_scale(Iy, gsl_complex_rect(0.0, 0.5));

  // gsl_matrix_complex_fprintf(stdout, Iy, "%g");
  // gsl_matrix_complex_fprintf(stdout, Hb, "%g");

  // Now that I've got my I,J x,y components
  // I need to build on 3.20 from Dan's thesis
  gsl_matrix_complex_add(Hbt, Ix);
  gsl_matrix_complex_scale(Hbt, g_Ip);

  gsl_matrix_complex_add(tempt, Jx);
  gsl_matrix_complex_scale(tempt,
                           gsl_complex_rect(calc_gj(1.0, L*2, 1.0), 0.0));

  gsl_matrix_complex_add(Hbt, tempt);
  // Calculating in units of MHz
  gsl_matrix_complex_scale(Hbt,
                       gsl_complex_rect((_bohr_magneton*field.B_x/_planck_h),
                                        0.0));


  // *****************************************************************

  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1.0, 0.0),
                 &I_z_view.matrix, &J_z_view.matrix, gsl_complex_rect(0.0, 0.0),
                 &H_view.matrix);

  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(0.5, 0.0),
                 &I_plus_view.matrix, &J_minus_view.matrix,
                 gsl_complex_rect(1.0, 0.0), &H_view.matrix);

  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(0.5, 0.0),
                 &I_minus_view.matrix, &J_plus_view.matrix,
                 gsl_complex_rect(1.0, 0.0), &H_view.matrix);

  gsl_matrix_complex_scale(&H_view.matrix, gsl_complex_rect(Aj, 0.0));


  // gsl_matrix_complex_fprintf(stdout, Hbt, "%g");
  // printf("*******************\n");
  // gsl_matrix_complex_fprintf(stdout, Hbz, "%g");

  gsl_matrix_complex_add(&H_view.matrix, Hbz);
  gsl_matrix_complex_add(&H_view.matrix, Hbt);


  // printf("\nH=\n[");
  // for(int i = 0; i < (2*numBasisStates*numBasisStates)-1; i+=2) {
  //   printf("%8.6G\t", H[i]);
  //   if( (i+2) % (2*numBasisStates) == 0) {
  //     printf("\n");
  //   }
  // }
  // printf("\n\n\n");


  gsl_vector *eval = gsl_vector_alloc(numBasisStates);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(numBasisStates,
                                                      numBasisStates);
  gsl_eigen_hermv_workspace * w = gsl_eigen_hermv_alloc(numBasisStates);

  // ALL THE MAGIC
  gsl_eigen_hermv(&H_view.matrix, eval, evec, w);
  gsl_eigen_hermv_free(w);
  gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

  int *F2 = new int[numBasisStates];
  int *Fz2 = new int[numBasisStates];

  int F2current = abs(atom.I2 - J2);
  int Fz2current = -F2current;

  // Assumes weak field.
  // Assumes that each F-manifold is well seperated and that the energies
  // increase with increasing F
  for (int i = 0; i < numBasisStates; i++) {
    F2[i] = F2current;
    Fz2current += 2;
    if (Fz2current > F2current) {
      F2current += 2;
      Fz2current = -F2current;
    }
  }

  for (int i = 0; i < numBasisStates; i++) {
    // printf("i = %d", i);
    // Now for the Fz value
    int kWithMaxAdmix = -1;
    double max_real_admix = -1.0;
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);
      for (int k = 0; k < numBasisStates; k++) {
        gsl_complex adMix = gsl_vector_complex_get(&evec_i.vector, k);
        // printf("\tk = %d, trying adMix = %g\n", k, GSL_REAL(adMix));
        if (fabs(GSL_REAL(adMix)) > max_real_admix) {
          // printf("\t\t...its the new biggest\n");
          max_real_admix = fabs(GSL_REAL(adMix));
          kWithMaxAdmix = k;
        }
      }
      Fz2[i] = Iz2[2*(kWithMaxAdmix)] + Jz2[2*(kWithMaxAdmix)];
    // genAtomicState::Fz2vec[L][i] = setFz2;
  }

  // // here
  // for(int i = 0; i < numBasisStates; i++) {
  //   double eval_i = gsl_vector_get(eval,i);
  //   gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec,i);
  //   printf("eigenvalue = %8.6G ", eval_i);
  //   printf("\t\t | %d/2, %d/2 > = eigenvector = \n", F2[i], Fz2[i]);
  //   for(int k = 0; k < numBasisStates; k++) {
  //     gsl_complex adMix = gsl_vector_complex_get(&evec_i.vector,k);
  //     if( fabs(GSL_REAL(adMix)) > pow(10,-10)) {
  //       printf("%8.6G | %d", GSL_REAL(adMix), Iz2[2*k]);
  //       printf("/2, %d/2 > + ", Jz2[2*k]);
  //     }
  //   }
  //   printf("\n\n");
  //   //gsl_vector_fprintf(stdout,&evec_i.vector,"%g");
  // }
  // // here

  int F2desired = abs(atom.I2 - J2);
  int Fz2desired = -F2desired;
  /*
  double **admixture = new double*[numBasisStates];
  for (int i = 0; i < numBasisStates; i++) {
    admixture[i] = new double[numBasisStates];
  }
  */
  vector<vector<double> > admixture(numBasisStates,
                                    vector<double>(numBasisStates, 0.0));

  for (int i = 0; i < numBasisStates; i++) {
    //  printf("Looking for F2 = %i and Fz2 = %i \n",F2desired, Fz2desired);
    for (int j = 0; j < numBasisStates; j++) {
      if (F2[j] == F2desired && Fz2[j] == Fz2desired) {
        //  printf("FOUND IT!\n");
        gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, j);
        for (int k = 0; k < numBasisStates; k++) {
          gsl_complex adMix = gsl_vector_complex_get(&evec_i.vector, k);
          admixture[k][i] = GSL_REAL(adMix);
        }
      }
    }
    Fz2desired += 2;
    if (Fz2desired > F2desired) {
      F2desired+=2;
      Fz2desired = -F2desired;
    }
  }
  for (int i = 0; i < numBasisStates; i++) {
    for (int j = 0; j < numBasisStates; j++) {
      //      *(decomp_array+j+(i*numBasisStates)) = admixture[i][j];
      //  printf("%10.6G\t ",admixture[i][j]);
    }
    //  printf("\n");
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
  vector<double> eShift(numBasisStates, 0.0);
  for (int i = 0; i < numBasisStates; i++) {
    double eval_i = gsl_vector_get(eval, i);
    eShift[i] = eval_i;
    if (debug) printf("%12.10G\t", eval_i/_MHz);
  }
  if (debug) printf("\n\n");
  for (int i = 0; i < numBasisStates; i++) {
    for (int j = 0; j < numBasisStates; j++) {
      gsl_vector_complex_view evec_j = gsl_matrix_complex_column(evec, j);
      gsl_complex ev = gsl_vector_complex_get(&evec_j.vector, i);
      if (debug) printf("%12.10G\t", GSL_REAL(ev));
    }
    if (debug) printf("\n");
  }
  if (debug) printf("\n");

  if (L == 0) {
    eShift_ground = eShift;
  } else if (L == 1) {
    eShift_excited = eShift;
  } else {
    printf("Eigenvector_Helper::diagH --> L must be 0 or 1\n\n");
    exit(1);
  }
  // printf("for l = 0\n");
  // for (int i = 0; i < numBasisStates; i++) {
  //   printf("%g\n", eShift_ground[i]);
  // }

  //  delete[] Iz2;
  //  delete[] Jz2;
  delete[] I_z;
  delete[] J_z;
  delete[] I_plus;
  delete[] J_plus;
  delete[] I_minus;
  delete[] J_minus;
  delete[] H;
  delete[] F2;
  delete[] Fz2;
  return admixture;
  //  for (int i = 0; i < numBasisStates; i++)  delete[] admixture[i];
  // genAtomicState::eval[L] = eval;
  // genAtomicState::evec[L] = evec;
  // gsl_vector_free(eval);
  // gsl_matrix_complex_free(evec);
}

double Eigenvector_Helper::calc_gj(int J2, int L2, int S2) {
  double g_J = 0.0;
  if (J2 != 0) {
    double g_L = 1.0;  // orbital g-factor
    double g_s = 2.0023;  // electron spin g-factor
    double J = static_cast<double>(J2)/2.0;
    double L = static_cast<double>(L2)/2.0;
    double S = static_cast<double>(S2)/2.0;
    double num = ((J*(J+1.0)) + (L*(L+1.0)) - (S*(S+1.0)))*g_L;
    double den = 2*J*(J+1.0);
    g_J = num/den;
    num = ((J*(J+1.0)) - (L*(L+1.0)) + (S*(S+1.0)))*g_s;
    g_J += (num/den);
  }
  // printf("g_J = %6.4G\t", g_J);
  return g_J;
}

double Eigenvector_Helper::calc_gf(int F2, int J2, int I2, int L2,
                                   int S2, double g_I) {
  double g_F = 0.0;
  double g_J = calc_gj(J2, L2, S2);
  if (F2 != 0) {
    double memn = 1.0/1836.152701;  // ratio of electron and nucleon masses
    double g_Ip = -g_I*memn;
    double F = static_cast<double>(F2)/2.0;
    double J = static_cast<double>(J2)/2.0;
    double I = static_cast<double>(I2)/2.0;

    double num = ((F*(F+1.0)) + (J*(J+1.0)) - (I*(I+1.0)))*g_J;
    double den = 2*F*(F+1.0);
    g_F = num/den;
    num = ((F*(F+1.0)) - (J*(J+1.0)) + (I*(I+1.0)))*g_Ip;
    g_F += (num/den);
  }
  // printf("g_F = %6.4G\t", g_F);
  return g_F;
}

vector<int> Eigenvector_Helper::get_nuclear_spin(int numBasisStates) {
  /*
  int *Iz2_local = new int[numBasisStates*2]; // Times two to leave room for
  // imaginary component
  for (int i = 0; i < numBasisStates*2; i+=2) {
    Iz2_local[i] = -I2 + 2*((i/2)%(I2+1));
    Iz2_local[i+1] = 0.0;       // The imaginary part

    *(Iz2+i) = Iz2_local[i];
    *(Iz2+i+1) = Iz2_local[i+1];
  }
  */
  vector<int> Iz2_local(numBasisStates*2, 0);
  for (int i = 0; i < numBasisStates*2; i+=2) {
    Iz2_local[i] = -atom.I2 + 2*((i/2)%(atom.I2+1));
    Iz2_local[i+1] = 0;
  }
  return (Iz2_local);
}

vector<int> Eigenvector_Helper::get_total_atomic_spin(int numBasisStates,
                                                          int J2) {
  /*
  int *Jz2_local = new int[numBasisStates*2]; // Times two to leave room for
  // imaginary component
  for (int i = 0; i < numBasisStates*2; i+=2) {
    Jz2_local[i] = -J2 + 2*(i/2/(I2+1));
    Jz2_local[i+1] = 0.0;       // The imaginary part

    *(Jz2+i) = Jz2_local[i];
    *(Jz2+i+1) = Jz2_local[i+1];
  }
  */
  vector<int> Jz2_local(numBasisStates*2, 0);
  for (int i = 0; i < numBasisStates*2; i+=2) {
    Jz2_local[i] = -J2 + 2*(i/2/(atom.I2+1));
    Jz2_local[i+1] = 0;
  }
  return Jz2_local;
}

