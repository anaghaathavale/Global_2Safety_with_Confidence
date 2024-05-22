# Installation
## Build dependencies
1. The build process uses [CMake](https://cmake.org/download/) version 3.16 (or later). 

2. Depends on the Boost and the OpenBLAS libraries, which are automatically downloaded and built when you run make. Library CXXTEST comes included in the repository.

Instructions for Linux machines using CMake:
```
git clone https://github.com/anaghaathavale/Global_2Safety_with_Confidence.git
cd Marabou/
mkdir build 
cd build
cmake ../
cmake --build ./
```
## To run experiments with confidence-based global robustness:
1. For DNNs with 2 output classes:
   To run Query 1:

```
  cp src/engine/Marabou_robustness_Query1.cpp src/engine/Marabou
  cd build 
  make
  cd ..
  ./build/Marabou resources/nnet/'network_name.nnet' 'input_distance(epsilon)_value'
  
```
2. To run Query 2:
```
  cp src/engine/Marabou_robustness_Query1.cpp src/engine/Marabou
  cd build 
  make
  cd ..
  ./build/Marabou resources/nnet/'network_name.nnet' 'input_distance(epsilon)_value'
  
```
Note that in our experiments, the DNNs German Credit, Adult and Law have 2 output classes.

