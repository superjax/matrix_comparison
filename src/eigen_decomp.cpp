#include <stdio.h>
#include <math.h>
#include <iostream>
#include <chrono>
#include <cstdlib>
#include <unistd.h>

#include "TNT/tnt.h"
#include "JAMA/jama_eig.h"

#include "Eigen/Eigenvalues"

#include "../Eigen/Eigen/Dense" 
typedef Eigen::Matrix<float, 10, 10, Eigen::RowMajor> MagCalM10x10; 
typedef Eigen::Matrix<float, 6, 4, Eigen::RowMajor> MagCalM6x4; 
typedef Eigen::Matrix<float, 4, 4, Eigen::RowMajor> MagCalM4x4; 
typedef Eigen::Matrix<float, 4, 6, Eigen::RowMajor> MagCalM4x6; 
typedef Eigen::Matrix<float, 3, 3, Eigen::RowMajor> MagCalM3x3; 
typedef Eigen::VectorXcf MagCalVc; 
typedef Eigen::Vector3f  MagCalV3;

typedef Eigen::Matrix<float, 6, 6, Eigen::RowMajor> MagCalM6x6; 
typedef Eigen::VectorXf  MagCalV; 
typedef Eigen::ArrayXf   MagCalA; 
typedef Eigen::MatrixXf  MagCalM; 
typedef Eigen::MatrixXcf MagCalMc; 

MagCalM6x6 AE; 
MagCalMc VE(6, 6); 
MagCalA lambda6E(6); 
Eigen::EigenSolver<MagCalM> ES; 

typedef float f_t;

namespace TNT
{

Array1D<f_t> get_col(int col, Array2D<f_t>& array)
{
    int size = array.dim1();
    Array1D<f_t> out(size);
    for (int i = 0; i < size; i++)
    {
        out[i] = array[i][col];
    }
    return out;
}

float get_max_coeff(int* index, Array1D<f_t>& array)
{
    *index = 0;
    double max_e_value = array[0];
    for (int i = 1; i < array.dim1(); i++)
    {
        if (array[i] > max_e_value)
        {
            *index = i;
            max_e_value = array[i];
        }
    }
    return max_e_value;
}

void conv_to_c_arr(f_t* c_arr, const Array1D<f_t>& tnt_arr, int len)
{
    int i = 0;
    while (i < len)
    {
        *(c_arr+i) = tnt_arr[i++];
    }
}

void conv_to_2DTNT(const f_t* c_arr, Array2D<f_t>& tnt_arr, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            tnt_arr[i][j] = c_arr[i*rows + j];
        }
    }
}

void find_largest_eigen_vector(float* S11, float *v, float* max_eval)
{
    Array1D<f_t> evalues(6);
    Array2D<f_t> evectors(6, 6);
    Array2D<f_t> mat(6, 6);
    conv_to_2DTNT(S11, mat, 6, 6);
    JAMA::Eigenvalue<f_t> es_tnt(mat);
    es_tnt.getRealEigenvalues(evalues);
    es_tnt.getV(evectors);

    int ii;
    get_max_coeff(&ii, evalues);
    *max_eval = evalues[ii];

    Array1D<f_t> max_tnt_evec = get_col(ii, evectors);
    conv_to_c_arr(v, max_tnt_evec, 6);
}
}

int compare_eigen_TNT_largest_eigenvalue_calc()
{    
  float S11[36];
  float v_eigen[6];
  float v_TNT[6];
  
  // Initialize Matrix
  std::cout << "[";
  for (int i = 0; i < 6; i++)
  {
    std::cout << "\n [";
    for (int j = 0; j < 6; j++)
    {
//      S11[i+6*j] = (double)(rand())/double(RAND_MAX)*2.-1.; // this one almost always gives bogus results
//      S11[i+6*j] = (double)(rand())/double(RAND_MAX); // this one sometimes flips the e-vector
      S11[i+6*j] = (double)(rand())/double(RAND_MAX)-0.5; // This one sometimes gives bogus results
      if (i == j)
        S11[i+6*j] += 10.;
      std::cout << S11[i+6*j] << ", ";
    }
    std::cout << "],";
  }
  std::cout << "\n]" << std::endl;
  
  // Eigen Code
  AE = Eigen::Map<MagCalM6x6>(S11, 6, 6); 
  ES.compute(AE); 
  lambda6E = ES.eigenvalues().real(); 
  VE = ES.eigenvectors(); 
  
  int ii;
  lambda6E.maxCoeff(&ii); 
  float e_max_eval = lambda6E[ii];
  MagCalV::Map(v_eigen, 6) = VE.col(ii).real(); 
  
  // TNT Version
  float TNT_max_eval;
  TNT::find_largest_eigen_vector(S11, v_TNT, &TNT_max_eval);
  
  // Check equivalence
  std::cout << "eigen\tTNT\n";
  bool success = true;
  
  float v_TNT_flipped[6];
  for (int i = 0; i < 6; i++)
  {
    v_TNT_flipped[i] = v_TNT[i] * -1.0;
  }
  
  for (int i = 0; i < 6; i++)
  {
    std::cout << v_eigen[i] << "\t" << v_TNT[i]  << "\n";
    if (std::fabs(v_eigen[i] - v_TNT[i]) > 1e-4 
        && std::fabs(v_eigen[i] - v_TNT_flipped[i]) > 1e-4 )
    {
        success = false;
    }
  }
  std::cout << std::endl;
  if (!success)
    std::cout << "FAILED\n\n";
  else
    std::cout << "PASSED\n\n";
}

int main()
{
  for (int i = 0; i < 10; i++)
  {
    compare_eigen_TNT_largest_eigenvalue_calc();
  }
}

