#ifndef _FFT_H_
#define _FFT_H_


#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include "ca.h"


using namespace std;




class fft
{
public:
	void AllocateGrid(int sizex, int sizey);

	void ImportGrid(int **sigma);

	void ImportGrid(int **sigma, CellularPotts *cpm);

	void PolarTransform();

	void PolarToOutput(string name="polar.png");

	double PolarComparison(int** polar2);

	void FFTransform();

	void ShiftGrid(int **toshift, int n=1);

	void ReflectGrid(int **toshift);

	inline int** GetPolar()
	{
		return polar;
	}

	~fft(void);



private:
	int **grid;
	int **polar;
	int **tmp_polar;

	int rho = 360;
	// THIS IS HARD CODED, BE CAREFUl!!!
	int sizer=125;

	int sizex;
	int sizey;

};





#endif