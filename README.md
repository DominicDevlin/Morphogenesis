# Morphogenesis
Cellular Potts Model for the evolution of morphogenesis

Requires Qt libraries.

To compile:
Qmake -> Make -> executable
Target file in "CellularPotts2.pro" must be either sorting.cpp or evolution.cpp. 

"Sorting" runs the a cellular potts model of the development of single organism. The gene regulatory network used
can be randomised or manually put in.

"Evolution" runs a multithreaded simulation for the evolution of morphogenesis.
Many organisms are simulated simultaneously. 

Parameters are changed in "parameter.cpp".

**Beware!** 
There are about 100 parameters that covary. 



