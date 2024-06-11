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
  ./build/Marabou resources/nnet/'network_name.nnet' resources/properties/smallNnet.txt
  
```
where smallNnet.txt is the property file with confidence and epsilon values set by the user

To run Query 2:
```
  cp src/engine/Marabou_robustness_Query2.cpp src/engine/Marabou
  cd build 
  make
  cd ..
  ./build/Marabou resources/nnet/'network_name.nnet' resources/properties/smallNnet.txt
  
```
Note that in our experiments, the DNNs German Credit, Adult and Law have 2 output classes.

2. For DNNs with 3 output classes:
Run Query X = 1, 2 and 3:

```
  cp src/engine/Marabou_robustness_QueryX.cpp src/engine/Marabou
  cd build 
  make
  cd ..
  ./build/Marabou resources/nnet/'network_name.nnet' resources/properties/smallNnet.txt
  
```
Note that in our experiments, the COMPAS DNN has 3 output classes.
