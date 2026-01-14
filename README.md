# Morphogenesis
Cellular Potts Model for the evolution of morphogenesis.

The source code uses the "Tissue Simulation Toolkit" as a base and requires QT libraries. 

To compile:
Qmake -> Make -> executable
The target file can be changed in "CellularPotts2.pro." The target files are:
sorting.cpp - visualising a single development
evolution.cpp - evolutionary simulation described in the main text.
multisort.cpp - simulate the development of the same organism N times. Used for reproducibility scores.
potency.cpp - simulate the development of the same organisms N times. Used for making cell state spaces.
 
All relevant parameters, including the GRN weights, numbers of genes and types of genes are all changed in "parameter.cpp".





