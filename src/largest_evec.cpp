#include <stdio.h>
#include <chrono>
#include <cstdlib>
#include <iostream>

#include "TNT/tnt.h"
#include "JAMA/jama_eig.h"

#include "Eigen/Eigenvalues"

template<typename T>
TNT::Array1D<T> get_col(int col, TNT::Array2D<T>& array)
{
  int size = array.dim1();
  TNT::Array1D<T> out(size);
  for (int i = 0; i < size; i++)
  {
    out[i] = array[i][col];
  }
  return out;
}

bool check_close(double a, double b)
{
  if (fabs(fabs(a) - fabs(b)) > 0.0000001)
  {
    printf("problem\n");
  }
}


int main()
{
  int size = 6;
  TNT::Array1D<double> evalues(size, size);
  TNT::Array2D<double> evectors(size, size);

  Eigen::MatrixXd mat(size, size);
  TNT::Array2D<double> pei(size, size);
  printf("iter: %d/25", 0);
  for (int iter = 0; iter < 2500; iter++)
  {
    printf("\riter: %d/2500", iter+1);


    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        pei[i][j] = (double)(rand())/double(RAND_MAX);
        mat(i,j) = pei[i][j];
      }
    }

    for (int i = 0; i < size; i++)
    {
      pei[i][i] += 1.0;
      mat(i,i) += 1.0;
    }

    JAMA::Eigenvalue<double> es_tnt(pei);
    es_tnt.getRealEigenvalues(evalues);
    es_tnt.getV(evectors);

    int max_index = 0;
    double max_e_value = evalues[0];
    for (int i = 1; i < evalues.dim1(); i++)
    {
      if (evalues[i] > max_e_value)
      {
        max_index = i;
        max_e_value = evalues[i];
      }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> es;
    es.compute(mat, true);
    Eigen::MatrixXcd eig_evectors = es.eigenvectors();
    Eigen::VectorXd eig_evalues = es.eigenvalues().real();
    double max_tnt_eval = eig_evalues[max_index];
    TNT::Array1D<double> max_tnt_evec = get_col(max_index, evectors);

    //  std::cout << "TNT mat\n" << pei << "\nTNT real evals\n" << evalues << "\nTNT evecs\n" << evectors << "\nmax_eval_index\n" << max_index << "\nmax_tnt_evec\n" << max_tnt_evec << "\n";

    int max_index_eig = 0;
    eig_evalues.maxCoeff(&max_index_eig);
    double max_eig_eval = eig_evalues(max_index_eig);
    Eigen::VectorXd max_eig_evec = eig_evectors.col(max_index_eig).real();

    check_close(max_tnt_eval, max_eig_eval);
    for (int i = 0; i < size; i++)
    {
      check_close(max_tnt_evec[i], max_eig_evec(i));
    }
  }

  //  std::cout << "Eig mat\n" << mat << "\nreal Eig evals\n" << eig_evalues << "\nEig evecs\n" << eig_evectors << "\nmax_eval_index\n" << max_index_eig << "\nmax_eig_evec\n" << max_eig_evec << "\n";
}

