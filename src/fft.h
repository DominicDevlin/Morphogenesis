#ifndef _FFT_H_
#define _FFT_H_


#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>


using namespace std;




class fft
{
public:
	void AllocateGrid(int sizex, int sizey);

	void ImportGrid(int **sigma);

	void PolarTransform();

	void PolarToOutput();





private:
	int **grid;
	int **polar;

	int rho = 180;

	int sizex;
	int sizey;

};





#endif