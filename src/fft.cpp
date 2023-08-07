#include "fft.h"
#include "parameter.h"
#include <math.h>
#include <stdio.h>
#include "crash.h"
#include "x11graph.h"


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
	cout << "here1" << endl;
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
			int x = center[0] + round(r*cos(q));
			int y = center[1] + round(r*sin(q));
			if (x >= sizex -1 || x <= 0 || y >= sizex-1 || y <= 0)
			{
				polar[q][r] = 0;
			}
			else
				polar[q][r] = grid[x][y];
		}
	}
	cout << "here2" << endl;

}



void fft::PolarToOutput()
{
	for (int q=0;q<rho;++q)
	{
		for (int r=0;r<sizex;++r)
		{
			cout << polar[q][r];
		}
		cout << endl;
	}


  X11Graphics g(par.sizex*2,par.sizey*2);


	char fname[200];
	sprintf(fname,"polar.png",par.datadir);

	g.BeginScene();
	g.ClearImage();    


  for (int i = 0; i < rho-1; i++ )
    for (int j = 0; j < sizex-1; j++ ) 
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
      
      if (polar[i][j]>0)  /* if draw */ 
        g.Point( colour, 2*i, 2*j);
		}

	g.EndScene();

	g.Write(fname);



}


















