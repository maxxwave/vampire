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

#include <iostream>
#include <cstdlib>
#include <cmath>

// include headers
#include "atoms.h"
#include "errors.h"
#include "LLG.hpp"
#include "material.hpp"

// functions
int calculate_spin_fields(const int, const int);
int calculate_external_fields(const int, const int);
namespace LLG_arrays{
   std::vector<double> x_euler_array;
   std::vector<double> y_euler_array;
   std::vector<double> z_euler_array;

   std::vector<double> x_heun_array;
   std::vector<double> y_heun_array;
   std::vector<double> z_heun_array;

	std::vector <double> x_spin_storage_array;
	std::vector <double> y_spin_storage_array;
	std::vector <double> z_spin_storage_array;

   std::vector<double> x_initial_spin_array;
   std::vector<double> y_initial_spin_array;
   std::vector<double> z_initial_spin_array;
   
   bool simplicit_set=true;
   return EXIT_SUCCESS;
}

namespace sim{

   // check calling of routine if error checking is activated
   if(err::check==true){std::cout<< "sim::Semi-Implicit has been called"<<std::endl;}
   using namespace LLG_arrays;

   if(LLG_set==false) sim::LLG_init();

   // Store initial spin positions
   for(int atom=0; atom<num_atoms;atom++){
      x_initial_spin_array[atom] = atoms::x_spin_array[atom];
      y_initial_spin_array[atom] = atoms::y_spin_array[atom];
      z_initial_spin_array[atom] = atoms::z_spin_array[atom];
   }

   // Calculate fields
   calculate_spin_fields(0,num_atoms);
   calculate_external_fields(0, num_atoms);

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

   // Declare prime spin
   double S_prime[3]={};

   for (int atom=0; atom <num_atoms; atom++){

      cons int imaterial=atoms::type_array[atom];
		const double alpha_oneplusalpha_sq = mp::material[imaterial].alpha_oneplusalpha_sq;
		const double one_oneplusalpha_sq = mp::material[imaterial].one_oneplusalpha_sq; // material specific alpha and gamma
      
      // Store local spin in sand local field H (similar as previous method)
		const double S[3] = {atoms::x_spin_array[atom],atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
		const double H[3] = {atoms::x_total_spin_field_array[atom]+atoms::x_total_external_field_array[atom],
                           atoms::y_total_spin_field_array[atom]+atoms::y_total_external_field_array[atom],
                           atoms::z_total_spin_field_array[atom]+atoms::z_total_external_field_array[atom]};
      // Calculate A
      A[0] = one_oneplusalpha_sq*H[0] + alpha_oneplusalpha_sq*(S[1]H[2]-S[2]H[1]);
      A[1] = one_oneplusalpha_sq*H[1] + alpha_oneplusalpha_sq*(S[2]H[0]-S[0]H[2]);
      A[2] = one_oneplusalpha_sq*H[2] + alpha_oneplusalpha_sq*(S[0]H[1]-S[1]H[0]);
     
      // declare a variable to store the square module of A
      double A_sq=A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
      double one_Asq = 1.0/(1+A_sq);

      double AxS[3]={A[1]*S[2]-A[2]*S[1], A[2]*S[0]-A[0]*S[2], A[0]*S[1]-A[1]*S[0]};

      double AxAxS[3]={
                        };

      // Calculate Cayley transformation
      Cay[0]=S[0] + one_Asq*(AxS[0] + 0.25*AxAxS[0]);
      Cay[1]=S[1] + one_Asq*(AxS[1] + 0.25*AxAxS[1]);
      Cay[2]=S[2] + one_Asq*(AxS[2] + 0.25*AxAxS[2]);
     
      // Calculate the predictor step
      S_new[0] = S[0] + Cay[0];
      S_new[1] = S[1] + Cay[1];
      S_new[2] = S[2] + Cay[2];

    }// end of for 


}// end of namespace 
