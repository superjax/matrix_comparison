# Matrix Comparison Code


A short C++ program to compare the performance of GSL, TNT and Eigen, and a python script to plot it.


## Installation

TNT and Eigen are header-only, so no installation is required.  GSL, on the other hand, is not, and therefore requires installation.  I'm using the [AMPL CMake-enabled version](https://github.com/ampl/gsl/tree/644e768630841bd085cb7121085a688c4ff424d0) of GSL, which is included as a submodule.


To install, perform the following operations

```bash
git clone https://github.com/superjax/matrix_comparison.git
cd matrix_comparison
git submodule update --init --recursive
cd lib/GSL
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4 -l4
sudo make install
```

To run the examples
```bash
cd matrix_comparison
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4 -l4
./matrix_test
```

To plot the comparison
``` bash
cd python
python plot_results.py
```
