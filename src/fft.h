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

	fft(void);

	fft(int sx, int sy);

	void AllocateGrid(int sizex, int sizey);

	void ImportGrid(int **sigma);

	void ImportCPM(CellularPotts *cpm);

	void PolarTransform(int xcen=0, int ycen=0);

	void PolarToOutput(string name="polar.png");

	void cpmOutput(string name="cpm.png");

	void GridToOutput(string name="grid.png");

	void ShowOptimal(string name = "shift.png");

	double PolarComparison(int** polar2, bool record=false);

	void FFTransform();

	void ShiftGrid(int **toshift, int n=1);

	void ReflectGrid(int **toshift);

	double StickSymmetry();

	inline int** GetPolar()
	{
		return polar;
	}

	void ClearPolar();

	void OutputLoss();

	~fft(void);



private:
	int **grid;
	int **polar;
	int **tmp_polar;

	int rho = 360;
	// THIS IS HARD CODED, BE CAREFUl!!! needs to be at least 2*a^2, where a^2 is half the size of x or y on the grid. Longer is better than shorter.
	int sizer=180;

	int sizex;
	int sizey;

	int optimal;

	CellularPotts *m_CPM;

	vector<double> loss;

};





#endif