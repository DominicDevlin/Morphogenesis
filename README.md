# Morphogenesis
Cellular Potts Model for the evolution of morphogenesis. Model in src. Most other folders are for visualsation, smaller models and creating graphs.

Requires Qt libraries. The source code uses the "Tissue Simulation Toolkit" as a base. 

To compile:
Qmake -> Make -> executable
Target file in "CellularPotts2.pro" must be either sorting.cpp, evolution.cpp, multisort.cpp, transition.cpp, potency.cpp etc. 

"Sorting" runs the a cellular potts model of the development of single organism. The gene regulatory network used
can be randomised or manually put in.

"Evolution" runs a multithreaded simulation for the evolution of morphogenesis.
Many organisms are simulated simultaneously. 

Parameters are changed in "parameter.cpp".





