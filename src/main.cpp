#include <iostream>
#include <chrono>

#include "Eigen/Core"
#include "Eigen/Eigenvalues"

#

#define EIGEN_TEST_X10(start, alpha) \
  EIGEN_TEST(start, alpha) \
  EIGEN_TEST(start + 5, alpha) \
  EIGEN_TEST(start + 10, alpha) \
  EIGEN_TEST(start + 15, alpha) \
  EIGEN_TEST(start + 20, alpha) \
  EIGEN_TEST(start + 25, alpha) \
  EIGEN_TEST(start + 30, alpha) \
  EIGEN_TEST(start + 35, alpha) \
  EIGEN_TEST(start + 40, alpha) \
  EIGEN_TEST(start + 45, alpha)

#define EIGEN_TEST_X50(start, alpha) \
  EIGEN_TEST_X10(start, alpha) \
  EIGEN_TEST_X10(start+50, alpha) \
  EIGEN_TEST_X10(start+100, alpha)

//#define EIGEN_TEST(size, alpha) \
//{ \
//  Eigen::Matrix<double, size, size> eigen_pei = alpha * Eigen::Matrix<double, size, size>::Identity()\
//  + Eigen::Matrix<double, size, size>::Constant(1.0);\
//  Eigen::EigenSolver<Eigen::Matrix<double, size, size>> es;\
//  double avg_sum = 0.0; \
//  int avg_count = 0; \
//  for (int j = 0; j < 10; j++)\
//{\
//  auto start = std::chrono::steady_clock::now();\
//  es.compute(eigen_pei, true);\
//  avg_sum += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();\
//  avg_count++;\
//  }\
//  std::cout << "EIGEN: size: " << size << " took " << avg_sum/avg_count << " ms" << std::endl;\
//  }



int main()
{
  for (int size = 5; size < 100; size++)
  {
    Eigen::MatrixXd eigen_pei(size, size);
    eigen_pei.setOnes();
    eigen_pei += Eigen::MatrixXd::Identity(size, size);
    Eigen::EigenSolver<Eigen::MatrixXd> es;
    double avg_sum = 0.0;
    for (int j = 0; j < 10; j++)
    {
      auto start = std::chrono::steady_clock::now();
      es.compute(eigen_pei, true);
      avg_sum += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
    }
    std::cout << "EIGEN: size: " << size << " took " << avg_sum/10 << " ms\n";
  }
}
