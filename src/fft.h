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

	void PolarToOutput(string name="polar.png");

	void PolarComparison(int** polar2);

	void FFTransform();

	void ShiftGrid(int **toshift, int n=1);

	inline int** GetPolar()
	{
		return polar;
	}



private:
	int **grid;
	int **polar;

	int rho = 360;

	int sizex;
	int sizey;

};





#endif