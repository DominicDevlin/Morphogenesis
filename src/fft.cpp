#include "fft.h"
#include "parameter.h"
#include <math.h>
#include <stdio.h>
#include "crash.h"
#include "string.h"

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif

#include <fftw3.h>


extern Parameter par;

using namespace std;





void fft::AllocateGrid(int sx, int sy)
{
	sizex = sx;
	sizey = sy;
  grid=(int **)malloc(sizex*sizeof(int *));
  if (grid==NULL)
    MemoryWarning();
  
  grid[0]=(int *)malloc(sizex*sizey*sizeof(int));
  if (grid[0]==NULL)  
    MemoryWarning();
  
  
  {for (int i=1;i<sizex;i++) 
    grid[i]=grid[i-1]+sizey;}
  
  /* Clear grid */
   {for (int i=0;i<sizex*sizey;i++) 
     grid[0][i]=0; }



  polar=(int **)malloc(rho*sizeof(int *));
  if (polar==NULL)
    MemoryWarning();
  
  polar[0]=(int *)malloc(rho*sizex*sizeof(int));
  if (polar[0]==NULL)  
    MemoryWarning();
  
  
  {for (int i=1;i<rho;i++) 
    polar[i]=polar[i-1]+sizex;}
  
  /* Clear grid */
   {for (int i=0;i<rho*sizex;i++) 
     polar[0][i]=0; }

}


void fft::ImportGrid(int **sigma)
{
	for (int x=0; x < sizex;++x)
	{
		for (int y = 0; y < sizey;++y)
		{
			grid[x][y] = sigma[x][y];
		}
	}

}




void fft::PolarTransform()
{
	// find c.o.m
	int center[] = {0,0};
  int xtotal{};
  int ytotal{};
  int mass{};

  for (int x=1;x<sizex;++x)
    for (int y=1;y<sizey;++y)
    {
      if (grid[x][y] > 0)
      {
        xtotal += x;
        ytotal += y;
        ++mass;
      }
    }
  
  center[0] = (double)xtotal / mass;
  center[1] = (double)ytotal / mass;
	

	// need to iterate through POLAR co-ordinate system and inverse to find what that point is on grid

	// going to make the grid 360 x 250, although we dont need 250 pixels. 


	for (int q = 0; q < rho; ++q)
	{
		for (int r=0; r < sizex;++r)
		{
			// need to turn q into radians
			double qn = q * M_PI / 180.;

			int x = center[0] + round(r*cos(qn));
			int y = center[1] + round(r*sin(qn));


			if (x >= sizex -1 || x <= 0 || y >= sizex-1 || y <= 0)
			{
				polar[q][r] = 0;
			}
			else
				polar[q][r] = grid[x][y];
		}
	}

}



void fft::PolarToOutput(string name)
{
	// for (int q=0;q<rho;++q)
	// {
	// 	for (int r=0;r<sizex;++r)
	// 	{
	// 		cout << polar[q][r];
	// 	}
	// 	cout << endl;
	// }
	cout << name << endl;

	char* nname = new char[name.size() + 1];
	strcpy(nname, name.c_str());


	#ifdef QTGRAPHICS

			
	QtGraphics g(rho,150);
 


	char fname[200];
	sprintf(fname, nname,par.datadir);

	g.BeginScene();
	g.ClearImage();    


  for (int i = 0; i < rho-1; i++ )
    for (int j = 0; j < 150; j++ ) 
		{
      

      int colour;
      if (polar[i][j]<=0) 
			{
				colour=0;
      } 
			else 
			{
				colour = 10;
      }
      
      if (polar[i][j]>0)
			{
				/* if draw */ 
				colour = polar[i][j] % 255;
				g.Point( colour, i, j);
			}  
        
		}

	g.EndScene();

	g.Write(fname);

	#endif

	delete[] nname;
}


void fft::ShiftGrid(int **toshift, int n)
{
	//store 0th array
	for (int i=0;i<n;++i)
	{
		cout << n << endl;
		int tmp[sizex];
		for (int l=0;l<sizex;++l)
		{
			tmp[l] = toshift[0][l];
		}		

		//shift grid across
		for (int l=0;l<sizex*rho-sizex;++l)
		{
			toshift[0][l]=toshift[0][l+sizex];
		}

		for (int l=0;l<sizex;++l)
		{
			toshift[rho-1][l] = tmp[l];
		}	
	}

}



// we dont need to be efficient so can just compare grids with changing rho
void fft::PolarComparison(int** polar2)
{

	vector<int> overlaps;
	double max_overlap=0.;
	
	//we will rotate polar, and keep polar2 constant
	int tmp_polar[rho][sizex];
	

	for (int i = 0; i < rho; ++i)
	{

		int overlap{};
		int outer{};
		double proportion{};

		int tmp_polar[rho][sizex];
		for (int x = 0;x<rho*sizex;++x)
		{
			tmp_polar[0][x] = polar[0][x];
		}


		for (int x = 0; x<rho*sizex;++x)
		{
			if (!tmp_polar[0][x] && !polar2[0][x])
			{
				continue;
			}
			else if (tmp_polar[0][x] > 0 && polar2[0][x] > 0)
			{
				++overlap;
			}
			else if (tmp_polar[0][x] != -1)
			{
				++outer;
			}
		}


		proportion = (double)overlap / (double)(overlap + outer);
		overlaps.push_back(proportion);
		// cout << "Overlap is: " << overlap << endl;
		if (proportion > max_overlap)
			max_overlap = proportion;
	}
	cout << "max overlap:  " << max_overlap << endl;


}




void fft::FFTransform()
{
	// double *in;
	// fftw_complex *out;



	double inputData[rho][sizex-1];
	for (int i=0;i<rho;++i)
		for (int j=0;j<sizex-1;++j)
		{
			if (polar[i][j])
				inputData[i][j]=1;
			else
				inputData[i][j]=0;
		}

 
    double realOut[rho][sizex-1];
    double imagOut[rho][sizex-1];
    double amplitudeOut[rho][sizex-1];
 
    int height = rho;
    int width = sizex-1;
 
    // Two outer loops iterate on output data.
    for (int yWave = 0; yWave < height; yWave++)
    {
        for (int xWave = 0; xWave < width; xWave++)
        {
            // Two inner loops iterate on input data.
            for (int ySpace = 0; ySpace < height; ySpace++)
            {
                for (int xSpace = 0; xSpace < width; xSpace++)
                {
                    // Compute real, imag, and ampltude.
                    realOut[yWave][xWave] += (inputData[ySpace][xSpace] * cos(
                            2 * M_PI * ((1.0 * xWave * xSpace / width) + (1.0
                                    * yWave * ySpace / height)))) / sqrt(
                            width * height);
                    imagOut[yWave][xWave] -= (inputData[ySpace][xSpace] * sin(
                            2 * M_PI * ((1.0 * xWave * xSpace / width) + (1.0
                                    * yWave * ySpace / height)))) / sqrt(
                            width * height);
                    amplitudeOut[yWave][xWave] = sqrt(
                            realOut[yWave][xWave] * realOut[yWave][xWave]
                                    + imagOut[yWave][xWave]
                                            * imagOut[yWave][xWave]);
                }
                cout << realOut[yWave][xWave] << " + " << imagOut[yWave][xWave]
                        << " i (" << amplitudeOut[yWave][xWave] << ")\n";
            }
        }
    }
}















