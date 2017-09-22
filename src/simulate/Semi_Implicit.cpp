// ---------------------------------------------------------------------
//
// (c) author: Razvan Ababei
// date: 14/09/2017
//---------------------------------------------------------------------
//
// This code contains Semi-implicit integration method
//
// References:
// [1] J. H. Mentink et al., 2010, JOP
// [2] M. Ellis et al., 2012, PRB
// -------------------------------------------------------------------------------
// The method is described by Matt Ellis et al. in the reference above mentioned
// The method respects several steps:
// 1) defining the A function as: A(S)= -prefac*Dt*[Heff(S) + lambda*S x Heff)]
// 2) calculate the Cayley transformation: Cay(A)S=S + (A x S + 0.5 A x A x S)/(1+0.25*|A|^2)
// 3) calculate the predictor step S_prime=Cay(A)S
// 4) repeat the procedure and calculate the new spin 
//
#include <iostream>
#include <cstdlib>
#include <cmath>

// include headers
#include "atoms.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "material.hpp"
#include "sim.hpp"
// functions
int calculate_spin_fields(const int, const int);
int calculate_external_fields(const int, const int);

namespace sim{

int Simplicit(){
   // check calling of routine if error checking is activated
   if(err::check==true){std::cout<< "sim::Semi-Implicit has been called"<<std::endl;}
   using namespace LLG_arrays;

   if(LLG_set==false) sim::LLGinit();

   // Store initial spin positions
   for(int atom=0; atom<atoms::num_atoms; atom++){
      x_initial_spin_array[atom] = atoms::x_spin_array[atom];
      y_initial_spin_array[atom] = atoms::y_spin_array[atom];
      z_initial_spin_array[atom] = atoms::z_spin_array[atom];
   }

   // Calculate fields
   calculate_spin_fields(0,atoms::num_atoms);
   calculate_external_fields(0, atoms::num_atoms);

   // Declare A-function
   //        -gamma Dt
   // A(S) = -------------- [H + lambda S x H ]
   //        (1+lambda^2)mu   

   double A[3]={};

   // Declare Cayley transformation
   //                 A x S + A x A x S 
   // Cay(A)S = S + ---------------------
   //                     1 + |A|^2

   double Cay[3]={};

   for (int atom=0; atom <atoms::num_atoms; atom++){
      // import the LLG pre-factors 
      const int imaterial=atoms::type_array[atom];
      const double one_oneplusalpha_sq = mp::dt/(1.0 + mp::material[imaterial].alpha*mp::material[imaterial].alpha);
      const double alpha_oneplusalpha_sq = one_oneplusalpha_sq * mp::material[imaterial].alpha;
      
      // Store local spin in sand local field H (similar as previous method)
      const double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};
      const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                           atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                           atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]};
      // Calculate A
      A[0] = one_oneplusalpha_sq*H[0] + alpha_oneplusalpha_sq*(S[1]*H[2]-S[2]*H[1]);
      A[1] = one_oneplusalpha_sq*H[1] + alpha_oneplusalpha_sq*(S[2]*H[0]-S[0]*H[2]);
      A[2] = one_oneplusalpha_sq*H[2] + alpha_oneplusalpha_sq*(S[0]*H[1]-S[1]*H[0]);
     
      // Declare a variable to store the square module of A
      double A_sq=A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
      double one_Asq = 1.0/(1.0+0.25*A_sq);

      double AxS[3]={A[1]*S[2]-A[2]*S[1], A[2]*S[0]-A[0]*S[2], A[0]*S[1]-A[1]*S[0]};

      double AxAxS[3]={A[1]*AxS[2]-A[2]*AxS[1],
                       A[2]*AxS[0]-A[0]*AxS[2],
                       A[0]*AxS[1]-A[1]*AxS[0]
                        };

      // Calculate Cayley transformation
      Cay[0]=S[0] + one_Asq*(AxS[0] + 0.5*AxAxS[0]);
      Cay[1]=S[1] + one_Asq*(AxS[1] + 0.5*AxAxS[1]);
      Cay[2]=S[2] + one_Asq*(AxS[2] + 0.5*AxAxS[2]);
     
      // Calculate the average of predictor and initial spin value
      x_spin_storage_array[atom] = 0.5*(Cay[0]+S[0]);
      y_spin_storage_array[atom] = 0.5*(Cay[1]+S[1]);
      z_spin_storage_array[atom] = 0.5*(Cay[2]+S[2]);

  }// end of for 

  // Copy new spins to spin array
  for(int atom=0;atom<atoms::num_atoms;atom++){
     atoms::x_spin_array[atom]=x_spin_storage_array[atom];
     atoms::y_spin_array[atom]=y_spin_storage_array[atom];
     atoms::z_spin_array[atom]=z_spin_storage_array[atom];
   }                          
   // Recalculate spin dependent fields
   calculate_spin_fields(0,atoms::num_atoms);
   calculate_external_fields(0, atoms::num_atoms);
  
  // Calculate new spin
  for(int atom=0; atom<atoms::num_atoms;atom++){
      // import the LLG pre-factors 
      const int imaterial=atoms::type_array[atom];
      const double one_oneplusalpha_sq = mp::dt/(1.0 + mp::material[imaterial].alpha*mp::material[imaterial].alpha);
      const double alpha_oneplusalpha_sq = one_oneplusalpha_sq * mp::material[imaterial].alpha;

      // store initial spins 
      double Si[3]={x_initial_spin_array[atom], y_initial_spin_array[atom],  z_initial_spin_array[atom]};

      // Store local spin in sand local field H (similar as previous method)
      const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
      const double H[3] = {atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom],
                           atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom],
                           atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom]};
      // Calculate A
      A[0] = one_oneplusalpha_sq*H[0] + alpha_oneplusalpha_sq*(S[1]*H[2]-S[2]*H[1]);
      A[1] = one_oneplusalpha_sq*H[1] + alpha_oneplusalpha_sq*(S[2]*H[0]-S[0]*H[2]);
      A[2] = one_oneplusalpha_sq*H[2] + alpha_oneplusalpha_sq*(S[0]*H[1]-S[1]*H[0]);
     
      // Declare a variable to store the square module of A
      double A_sq=A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
      double one_Asq = 1.0/(1.0+0.25*A_sq);

      double AxSi[3]={A[1]*Si[2]-A[2]*Si[1], A[2]*Si[0]-A[0]*Si[2], A[0]*Si[1]-A[1]*Si[0]};

      double AxAxSi[3]={A[1]*AxSi[2]-A[2]*AxSi[1],
                        A[2]*AxSi[0]-A[0]*AxSi[2],
                        A[0]*AxSi[1]-A[1]*AxSi[0]
                        };

      // Calculate Cayley transformation
      Cay[0]=Si[0] + one_Asq*(AxSi[0] + 0.5*AxAxSi[0]);
      Cay[1]=Si[1] + one_Asq*(AxSi[1] + 0.5*AxAxSi[1]);
      Cay[2]=Si[2] + one_Asq*(AxSi[2] + 0.5*AxAxSi[2]);
     
      // Calculate the predictor step
      x_spin_storage_array[atom] = Cay[0];
      y_spin_storage_array[atom] = Cay[1];
      z_spin_storage_array[atom] = Cay[2];
   }
  // Copy new spins to spin array
  for(int atom=0;atom<atoms::num_atoms;atom++){
     atoms::x_spin_array[atom]=x_spin_storage_array[atom];
     atoms::y_spin_array[atom]=y_spin_storage_array[atom];
     atoms::z_spin_array[atom]=z_spin_storage_array[atom];
  }// end of for 
}//end of Simplicit

}// end of namespace
