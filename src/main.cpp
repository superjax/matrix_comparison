#include <stdio.h>
#include <chrono>
#include <cstdlib>

#include "TNT/tnt.h"
#include "JAMA/jama_eig.h"

#include "Eigen/Eigenvalues"

#include "gsl/gsl_eigen.h"
#include "gsl/gsl_math.h"


double run_TNT_symm(int size, double alpha)
{
  TNT::Array2D<double> pei(size, size, (double)1.0);
  TNT::Array2D<double> id(size, size, (double)0);
  for (int i = 0; i < size; i++)
    id[i][i] = 1.;
  pei += id;

  TNT::Array2D<double> evalues(size, size);
  TNT::Array2D<double> evectors(size, size);

  double avg_sum = 0.0;
  for (int i = 0; i < 25; i++)
  {
    auto start = std::chrono::steady_clock::now();
    JAMA::Eigenvalue<double> es(pei);
    es.getD(evalues);
    es.getV(evectors);
    avg_sum += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
  }
  return avg_sum/25.0;
}

double run_Eigen_symm(int size, double alpha)
{
  Eigen::MatrixXd eigen_pei(size, size);
  eigen_pei.setOnes();
  eigen_pei += Eigen::MatrixXd::Identity(size, size);
  Eigen::EigenSolver<Eigen::MatrixXd> es;
  double avg_sum = 0.0;
  for (int j = 0; j < 25; j++)
  {
    auto start = std::chrono::steady_clock::now();
    es.compute(eigen_pei, true);
    Eigen::MatrixXcd evectors = es.eigenvectors();
    Eigen::MatrixXcd evalues = es.eigenvalues();
    avg_sum += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
  }
  return avg_sum/25.0;
}

double run_GSL_symm(int size, double alpha)
{
  double data[size*size];
  for (int i = 0; i < size*size; i++)
  {
    data[i] = 1.0;
  }
  for (int i = 0; i < size; i++)
  {
    data[i + i*size] += alpha;
  }

  double avg_sum = 0.0;
  gsl_matrix_view gsl_pei = gsl_matrix_view_array(data, size, size);

  gsl_vector *eval = gsl_vector_alloc(size);
  gsl_matrix *evec = gsl_matrix_alloc(size, size);
  gsl_eigen_symmv_workspace *wspace = gsl_eigen_symmv_alloc(size);
  for (int j = 0; j < 25; j++)
  {
    auto start = std::chrono::steady_clock::now();
    gsl_eigen_symmv (&gsl_pei.matrix, eval, evec, wspace);
    avg_sum += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
  }
  gsl_eigen_symmv_free (wspace);
  return avg_sum/25.0;
}

double run_TNT_asymm(int size)
{
  TNT::Array2D<double> evalues(size, size);
  TNT::Array2D<double> evectors(size, size);

  double avg_sum = 0.0;
  for (int i = 0; i < 25; i++)
  {
    TNT::Array2D<double> pei(size, size);
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        pei[i][j] = (double)(rand())/double(RAND_MAX);
      }
    }
    TNT::Array2D<double> id(size, size, (double)0);
    for (int i = 0; i < size; i++)
      id[i][i] = 1.;
    pei += id;

    auto start = std::chrono::steady_clock::now();
    JAMA::Eigenvalue<double> es(pei);
    es.getD(evalues);
    es.getV(evectors);
    avg_sum += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
  }
  return avg_sum/25.0;
}

double run_Eigen_asymm(int size)
{
  Eigen::EigenSolver<Eigen::MatrixXd> es;
  double avg_sum = 0.0;
  for (int j = 0; j < 25; j++)
  {
    Eigen::MatrixXd mat(size, size);
    mat.setRandom();
    mat += Eigen::MatrixXd::Identity(size, size);
    auto start = std::chrono::steady_clock::now();
    es.compute(mat, true);
    volatile Eigen::MatrixXcd evectors = es.eigenvectors();
    volatile Eigen::MatrixXcd evalues = es.eigenvalues();
    avg_sum += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
  }
  return avg_sum/25.0;
}

double run_GSL_asymm(int size)
{
  double data[size*size];

  double avg_sum = 0.0;

  gsl_vector_complex *eval = gsl_vector_complex_alloc(size);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(size, size);
  gsl_eigen_nonsymmv_workspace *wspace = gsl_eigen_nonsymmv_alloc(size);

  for (int j = 0; j < 25; j++)
  {
    for (int i = 0; i < size*size; i++)
    {
      data[i] = (double)(rand()/(double)(RAND_MAX));
    }
    for (int i = 0; i < size; i++)
    {
      data[i + i*size] += 1.0;
    }
    gsl_matrix_view gsl_pei = gsl_matrix_view_array(data, size, size);

    auto start = std::chrono::steady_clock::now();
    gsl_eigen_nonsymmv (&gsl_pei.matrix, eval, evec, wspace);
    avg_sum += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
  }
  gsl_eigen_nonsymmv_free (wspace);
  return avg_sum/25.0;
}


int main()
{
    FILE *symm_fp = fopen("csv_symm_results.csv", "w");
    printf("size,\tEigen,\tTNT\n");
    fprintf(symm_fp, "size,\tEigen,\tTNT\n");
    for (int size = 5; size < 200; size += 1)
    {
      double tnt_time = run_TNT_symm(size, 0.4);
      double eigen_time = run_Eigen_symm(size, 0.4);
      double gsl_time = run_GSL_symm(size, 0.4);
      printf("%*d,%*.2f,%*.2f,%*.2f\n",7, size, 12, eigen_time, 12, tnt_time, 12, gsl_time);
      fprintf(symm_fp, "%*d,%*.2f,%*.2f,%*.2f\n",7, size, 12, eigen_time, 12, tnt_time, 12, gsl_time);
    }
    fclose(symm_fp);

    FILE* asymm_fp = fopen("csv_asymm_results.csv", "w");
    printf("size,\tEigen,\tTNT\n");
    fprintf(asymm_fp, "size,\tEigen,\tTNT\n");
    for (int size = 5; size < 200; size += 1)
    {
      double tnt_time = run_TNT_asymm(size);
      double eigen_time = run_Eigen_asymm(size);
      double gsl_time = run_GSL_asymm(size);
      printf("%*d,%*.2f,%*.2f,%*.2f\n",7, size, 12, eigen_time, 12, tnt_time, 12, gsl_time);
      fprintf(asymm_fp, "%*d,%*.2f,%*.2f,%*.2f\n",7, size, 12, eigen_time, 12, tnt_time, 12, gsl_time);
    }
    fclose(asymm_fp);
}

