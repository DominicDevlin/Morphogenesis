#include "fft.h"
#include "parameter.h"
#include <math.h>
#include <stdio.h>
#include "crash.h"
#include "string.h"
#include "ca.h"
#include <fstream>

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif
#include <pthread.h>



// #include <fftw3.h>


extern Parameter par;

using namespace std;

fft::fft(void)
{
	sizex=250;
	sizey=250;
}

fft::fft(int sx, int sy)
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
  
  polar[0]=(int *)malloc(rho*sizer*sizeof(int));
  if (polar[0]==NULL)  
    MemoryWarning();
  
  
  {for (int i=1;i<rho;i++) 
    polar[i]=polar[i-1]+sizer;}
  
  /* Clear grid */
   {for (int i=0;i<rho*sizer;i++) 
     polar[0][i]=0; }


  tmp_polar=(int **)malloc(rho*sizeof(int *));
  if (tmp_polar==NULL)
    MemoryWarning();
  
  tmp_polar[0]=(int *)malloc(rho*sizer*sizeof(int));
  if (tmp_polar[0]==NULL)  
    MemoryWarning();
  
  
  {for (int i=1;i<rho;i++) 
    tmp_polar[i]=tmp_polar[i-1]+sizer;}
  
  /* Clear grid */
   {for (int i=0;i<rho*sizer;i++) 
     tmp_polar[0][i]=0; }
}




// destructor (virtual)
fft::~fft(void) {
  if (grid) {
    free(grid[0]);
    free(grid);
    grid=0;
  }

  if (polar)
  {
    free(polar[0]);
    free(polar);
    polar=0;
  }

	if (tmp_polar)
	{
		free(tmp_polar[0]);
		free(tmp_polar);
		tmp_polar=0;

	}
}


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
  
  polar[0]=(int *)malloc(rho*sizer*sizeof(int));
  if (polar[0]==NULL)  
    MemoryWarning();
  
  
  {for (int i=1;i<rho;i++) 
    polar[i]=polar[i-1]+sizer;}
  
  /* Clear grid */
   {for (int i=0;i<rho*sizer;i++) 
     polar[0][i]=0; }


  tmp_polar=(int **)malloc(rho*sizeof(int *));
  if (tmp_polar==NULL)
    MemoryWarning();
  
  tmp_polar[0]=(int *)malloc(rho*sizer*sizeof(int));
  if (tmp_polar[0]==NULL)  
    MemoryWarning();
  
  
  {for (int i=1;i<rho;i++) 
    tmp_polar[i]=tmp_polar[i-1]+sizer;}
  
  /* Clear grid */
   {for (int i=0;i<rho*sizer;i++) 
     tmp_polar[0][i]=0; }


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

// value of grid corresponds to cell type
void fft::ImportCPM(CellularPotts *cpm)
{
	m_CPM = cpm;
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


	// WE ARE GOING TO USE THE ORIGIN INSTEAD OF CENTER OF MASS FOR NOW!!! Replace par.sizex with center[0] and center[1] for com.


	for (int q = 0; q < rho; ++q)
	{
		for (int r=0; r < sizer;++r)
		{
			// need to turn q into radians
			double qn = q * M_PI / 180.;

			int x = round(par.sizex/2) + round(r*cos(qn));
			int y = round(par.sizey/2) + round(r*sin(qn));


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
	char* nname = new char[name.size() + 1];
	strcpy(nname, name.c_str());


	#ifdef QTGRAPHICS

			
	QtGraphics g(rho*2,sizer*2);
 


	char fname[200];
	sprintf(fname, nname,par.datadir);

	g.BeginScene();
	g.ClearImage();    


  for (int i = 0; i < rho-1; i++ )
    for (int j = 0; j < sizer; j++ ) 
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
				g.Point( colour, 2*i, 2*j);
				g.Point( colour, 2*i+1, 2*j+1);
				g.Point( colour, 2*i, 2*j+1);
				g.Point( colour, 2*i+1, 2*j);
			}  
        
		}

	g.EndScene();

	g.Write(fname);

	#endif

	delete[] nname;
}

void fft::cpmOutput(string name)
{
	char* nname = new char[name.size() + 1];
	strcpy(nname, name.c_str());

	#ifdef QTGRAPHICS

	QtGraphics g(sizex*2,sizey*2);

	char fname[200];
	sprintf(fname, "%s/%07.png",par.data_file.c_str(), nname);

	g.BeginScene();
	g.ClearImage();    

	m_CPM->SearchNandPlot(&g, false);

	g.EndScene();
	g.Write(fname);

	#endif

	delete[] nname;
}



void fft::GridToOutput(string name)
{
	char* nname = new char[name.size() + 1];
	strcpy(nname, name.c_str());


	#ifdef QTGRAPHICS

			
	QtGraphics qtgrid(sizex*2,sizey*2);
 


	char fname[200];
	sprintf(fname, nname,par.datadir);

	qtgrid.BeginScene();
	qtgrid.ClearImage();    


  for (int i = 0; i < sizex; i++ )
    for (int j = 0; j < sizey; j++ ) 
		{
      

      int colour{};
      if (grid[i][j]<=0) 
			{
				colour=0;
      } 
      
      if (grid[i][j]>=0)
			{
				/* if draw */ 
				colour = grid[i][j] % 255;
				qtgrid.Point( colour, 2*i, 2*j);
				qtgrid.Point( colour, 2*i+1, 2*j+1);
				qtgrid.Point( colour, 2*i, 2*j+1);
				qtgrid.Point( colour, 2*i+1, 2*j);
			} 
			// else if (!grid[i][j])
			// {

			// } 

        
		}

	qtgrid.EndScene();

	qtgrid.Write(fname);

	#endif

	delete[] nname;
}







void fft::ShowOptimal(string name)
{
	char* nname = new char[name.size() + 1];
	strcpy(nname, name.c_str());

	for (int i=0;i<rho*sizer;++i)
		tmp_polar[0][i]=polar[0][i];

	if (optimal>rho)
	{
		ReflectGrid(tmp_polar);
		int shifts = optimal % rho;
		ShiftGrid(tmp_polar, shifts);
	}
	else
	{
		ShiftGrid(tmp_polar, optimal);
	}

	#ifdef QTGRAPHICS

			
	QtGraphics s(rho*2,sizer*2);

	char fname[200];
	sprintf(fname, nname,par.data_file);

	s.BeginScene();
	s.ClearImage();    

  for (int i = 0; i < rho-1; i++ )
    for (int j = 0; j < sizer; j++ ) 
		{
      

      int colour;
      if (tmp_polar[i][j]<=0) 
			{
				colour=0;
      } 
			else 
			{
				colour = 10;
      }
      
      if (tmp_polar[i][j]>0)
			{
				// draw colour, double size for more pixels
				colour = tmp_polar[i][j] % 255;
				s.Point( colour, 2*i, 2*j);
				s.Point( colour, 2*i+1, 2*j+1);
				s.Point( colour, 2*i, 2*j+1);
				s.Point( colour, 2*i+1, 2*j);
			}  
        
		}

	s.EndScene();

	s.Write(fname);

	#endif

	delete[] nname;	
}







void fft::ShiftGrid(int **toshift, int n)
{
	//store 0th array
	for (int i=0;i<n;++i)
	{
		int tmp[sizer];
		for (int l=0;l<sizer;++l)
		{
			tmp[l] = toshift[0][l];
		}		

		//shift grid across
		for (int l=0;l<sizer*rho-sizer;++l)
		{
			toshift[0][l]=toshift[0][l+sizer];
		}
		
		//replace final line
		for (int l=0;l<sizer;++l)
		{
			toshift[rho-1][l] = tmp[l];
		}
	}

}

// reflect across midline
void fft::ReflectGrid(int **toshift)
{

	int tmp[rho][sizer];
	for (int i=0;i<rho*sizer;++i)
	{
		tmp[0][i] = toshift[0][i];
	}

	
	for (int i=0;i<rho;++i)
	{
		int swap = rho - i - 1;
		
		
		for (int j=0;j<sizer;++j)
		{
			toshift[i][j]=tmp[swap][j];
		}
	}
}


// we dont need to be efficient so can just compare grids with changing rho using dice coefficient
double fft::PolarComparison(int** polar2, bool record)
{

	vector<int> overlaps;
	double max_overlap=0.;
	double min_overlap=1.;


	int opt_shift{};
	int counter{};

	//initialise tmp polar to polar
	for (int i=0;i<rho*sizer;++i)
		tmp_polar[0][i]=polar[0][i];
	
	//we will rotate and reflect polar, and keep polar2 constant

	for (int i = 0; i < rho; ++i)
	{
		
		int overlap{};
		int outer{};
		double proportion{};

		// shift grid one degree over
		ShiftGrid(tmp_polar);

		// we have to scale radius
		for (int x = 0; x<rho;++x)
		{
			for (int r=0;r<sizer;++r)
			{
				if (!tmp_polar[x][r] && !polar2[x][r])
				{
					continue;
				}
				else if (tmp_polar[x][r] > 0 && polar2[x][r] > 0)
				{
					overlap += (r+1);
				}
				else if (tmp_polar[x][r] != -1)
				{
					outer += r+1;
				}
			}

		}

		proportion = (double)overlap / (double)(overlap + outer);
		overlaps.push_back(proportion);

		++counter;

		if (par.print_fitness)
			cout << "Count is... " << counter << "  Overlap is: " << overlap << "  Outer is: " << outer << "  Proportion is: " << proportion << endl;
			
		if (proportion > max_overlap)
		{
			max_overlap = proportion;
			opt_shift = counter;
		}

		if (proportion < min_overlap)
			min_overlap = proportion;

		if (record)
			loss.push_back(proportion);

		
	}

	// we need to now flip the grid to make sure it is reflection invariant as well. 
	ReflectGrid(tmp_polar);

	// repeat process
	for (int i = 0; i < rho; ++i)
	{

		int overlap{};
		int outer{};
		double proportion{};

		// shift grid one degree over
		ShiftGrid(tmp_polar);

		// polar method over compensates for things close to center, so we are going to remove that by starting
		// at a higher radius, especially because we want to capture morphology
		for (int x = 0; x<rho;++x)
		{
			for (int r=0;r<sizer;++r)
			{
				if (!tmp_polar[x][r] && !polar2[x][r])
				{
					continue;
				}
				else if (tmp_polar[x][r] > 0 && polar2[x][r] > 0)
				{
					overlap += r+1;
				}
				else if (tmp_polar[x][r] != -1)
				{
					outer += r+1;
				}
			}

		}

		proportion = (double)overlap / (double)(overlap + outer);
		overlaps.push_back(proportion);

		// cout << "Post-reflect, i is: " << i << " Overlap is: " << overlap << "  Outer is: " << outer << endl;
		++counter;
		if (proportion > max_overlap)
		{
			max_overlap = proportion;
			opt_shift = counter;
		}
		
		if (proportion < min_overlap)
			min_overlap = proportion;
		if (par.print_fitness)
			cout << "Count is... " << counter << "  Overlap is: " << overlap << "  Outer is: " << outer << "  Proportion is: " << proportion << endl;

		if (record)
			loss.push_back(proportion);		
		
	}
	// cout << "max overlap:  " << max_overlap << endl;
	// cout << "min overlap:  " << min_overlap << endl;

	if (par.print_fitness)
		cout << "optimal arrangement is: " << opt_shift << endl;

	optimal = opt_shift;
	return max_overlap;


}


void fft::OutputLoss()
{
	ofstream outfile;
  string var_name = "loss.txt";
  outfile.open(var_name, ios::app);

	for (unsigned int i=0;i<rho;++i)
	{
		outfile << i << '\t' << loss[i] << endl;
		
	}
	outfile.close();
}




double fft::StickSymmetry()
{
	// we want to check where the longest stick in the organism is, and see if it is also long on the other side. 
	// Anywhere else detracts. 
	double max_rho{};
	double max_r{};
	for (int i=0;i<rho;++i)
		for (int j=0;j<sizer-1;++j)
		{
			if (polar[i][j])
			{
				if (j > max_r)
				{
					max_r = j;
					max_rho = i;
				}
			}
		}


	double same_side_max = (max_rho + 40 < rho) ? max_rho + 40 : (int(max_rho) + 40) % rho;
	double same_side_min = (max_rho - 40 > 0) ? max_rho - 40 : rho + max_rho - 40;

	double opp = (max_rho < rho/2) ? max_rho + rho/2 : max_rho - rho/2;
	double opp_min = (opp - 40 > 0) ? opp - 40 : rho + opp - 40;
	double opp_max = (opp + 40 < rho) ? opp + 40 : (int(opp) + 40) % rho;

	double opp_r{};
	double opp_rho{};
	
	double bad_length{};

	// sub + add loop
	for (int i=0;i<rho;++i)
		for (int j=0;j<sizer-1;++j)
		{
			if (polar[i][j])
			{

				if ((same_side_max < same_side_min && (i > same_side_min || i < same_side_max)) || (i > same_side_min && i <  same_side_max))
				{
					continue;
				}
				else
				{
					if (opp_max < opp_min)
					{
						if (i > opp_min || i < opp_max)
						{
							// larger is better
							if (j > opp_r)
							{
								opp_r = j;
								opp_rho = i;
							}
						}
						else
						{
							// smaller is better
							bad_length += j*j;
						}
					}
					else
					{
						if (i > opp_min && i < opp_max)
						{
							// larger is better
							if (j > opp_r)
							{
								opp_r = j;
								opp_rho = i;
							}
						}
						else
						{
							// smaller is better
							bad_length += j*j;
						}
					}
				}



				
			}
		}



	bad_length = bad_length / 10000000;

	double fitness = (max_r + opp_r - 80) / bad_length;
	fitness = fitness;


	if (par.print_fitness)
	{
		cout << "Fitness: " << fitness << "  length: " << max_r << "  opp length: " << opp_r << endl;
		cout << "max rho: " << max_rho << "  opp rho: " << opp_rho << "  detraction: " << bad_length << endl;
	}


	return fitness;
}


void fft::FFTransform()
{
	// double *in;
	// fftw_complex *out;



	double inputData[rho][sizer-1];
	for (int i=0;i<rho;++i)
		for (int j=0;j<sizer-1;++j)
		{
			if (polar[i][j])
				inputData[i][j]=1;
			else
				inputData[i][j]=0;
		}

 
    double realOut[rho][sizer-1];
    double imagOut[rho][sizer-1];
    double amplitudeOut[rho][sizer-1];
 
    int height = rho;
    int width = sizer-1;
 
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
















