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

	void ShowOptimal(string name = "shift.png");

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
	// THIS IS HARD CODED, BE CAREFUl!!! needs to be 2*a^2, where a^2 is half the size of x or y on the grid. Longer is better than shorter.
	int sizer=180;

	int sizex;
	int sizey;

	int optimal;

};





#endif