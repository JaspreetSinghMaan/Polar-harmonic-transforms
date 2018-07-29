/******************************************************************************************************************
*																												  *
*	Developed By: Jaspreet Singh,                                                                                 *
*				  Research scholar (Ph.D.),                                                                       *
*				  Department of Computer Science,                                                                 *
*				  Punjabi University, Patiala, Punjab, India-151103.                                              *
*				  Mob.: +919803335164.                                                                            *
*				  Email: maan.jaspreetpitho@gmail.com                                                             *
*				  Date: 29/7/2018.                                                                                *
*  Instructions: 1) This program is implemented in Microsoft Visual Studio 2013. But it is also succesfully	      *
*					tested in Visual Studion 2008 and 2010 with minor changes.                                    *
*				 2) Basic purpose of this program is to read color images from the database text file, compute the*
*				    PCTs, QPCTs, and RQPCTs and write the features into a file.                                   *
*				 3) After succesfull build, allocate the memory to the program as follows: 1) click the "project" *
*				    button, click on "PHTs properties", click on "Linker", click on "System", allocate 100000000  *
*					in "Stack reserve size".                                                                      *
*				 4) Sample database file is provided "Colorimagereadwrite.txt" which contains two images.         *
*                5) If you have any issue regarding program execution feel free to communicate on the above       *
*					contact information.                                                                          *
*                                                                                                                 *
*******************************************************************************************************************
*/
#include "stdafx.h"
#include <iostream>
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <Windows.h>
#define MAXFEATURES 500
#define MAXIMAGES 1000
#define IMSIZE 1500
#define NMAX 15
#define EPS 1.0e-06
#define PI 3.1415926535897
using namespace std;

//Global variable section
const int R = 0, G = 1, B = 2;
int M, N, num, classno, imageno, total_classes, im_per_class, total_im, step = 2;
int label[MAXIMAGES], w_flag = 0, gradimage[IMSIZE][IMSIZE] = { 0 };
double dist[IMSIZE][IMSIZE], xi[IMSIZE], yk[IMSIZE], fact[100];
FILE *fout, *fe_file1, *fe_file2, *fe_file3;

/*****************************************************************************************************************/
/* P.T. Yap, X. Jiang, A. Chichung Kot, Two-dimensional polar harmonic transforms for invariant image            */
/* representation, IEEE Trans. Pattern Anal. Mach. Intell. 32 (2010) 1259–1270. doi:10.1109/TPAMI.2009.119.      */
/*****************************************************************************************************************/
int PCT(int imagenum, int nmax, int fxyr[IMSIZE][IMSIZE], int fxyg[IMSIZE][IMSIZE], int fxyb[IMSIZE][IMSIZE], 
	double Anm[3][2 * NMAX + 1][2 * NMAX + 1], double AnmReal[3][2 * NMAX + 1][2 * NMAX + 1], 
	double AnmImg[3][2 * NMAX + 1][2 * NMAX + 1], int OC)
{
	int i, k, icount, p, q, x, y;
	double rpq = 0, vpq = 0, theta, D, norm_fact, vv1[256][256] = { 0 }, 
		vv2[260][260] = { 0 };

	//Intialization section
	icount = 1; rpq = 0; vpq = 0;

	//For outer rectangle
	if (M != N)
		D = sqrt((double)M*M + (double)N*N);
	//For outer circle
	else if ((M == N) && OC == 1)
	{
		D = (N*sqrt(2.0));
	}
	//For Inner circle
	else
		D = N;

	for (i = 0; i<M; i++)
	{
		xi[i] = ((2.0 * i) + 1.0 - M) / D;
		for (k = 0; k<N; k++)
		{
			yk[k] = ((2.0 * k) + 1.0 - N) / D;
			dist[i][k] = sqrt((xi[i] * xi[i]) + (yk[k] * yk[k]));
		}
	}

	for (p = 0; p <= nmax; p++)
	{
		if (p == 0)
			norm_fact = (1.0 / PI)*(4.0 / D*D);
		else
			norm_fact = (2.0 / PI)*(4.0 / D*D);
		x = p;
		for (q = 0; q <= nmax; q++)
		{
			//y=p+q;
			y = q;
			AnmReal[R][x][y] = 0.0, AnmImg[R][x][y] = 0.0, Anm[R][x][y] = 0.0;
			AnmReal[G][x][y] = 0.0, AnmImg[G][x][y] = 0.0, Anm[G][x][y] = 0.0;
			AnmReal[B][x][y] = 0.0, AnmImg[B][x][y] = 0.0, Anm[B][x][y] = 0.0;
			for (i = 0; i < M; i++)
			{
				for (k = 0; k < N; k++)
				{
					if (dist[i][k] > 1.0)
						continue;
					rpq = cos(PI * p * dist[i][k] * dist[i][k]);

					theta = atan2(yk[k], xi[i]);
					if (theta < 0.0)
						theta = (2.0*3.14159265) + theta;
					vv1[x][y] = rpq*cos(q*theta);
					vv2[x][y] = rpq*sin(q*theta);
					AnmReal[R][x][y] = AnmReal[R][x][y] + (4.0 * norm_fact * (vv1[x][y] * fxyr[i][k])) / (D*D);
					AnmImg[R][x][y] = AnmImg[R][x][y] - (4.0 * norm_fact * (vv2[x][y] * fxyr[i][k])) / (D*D);
					AnmReal[G][x][y] = AnmReal[G][x][y] + (4.0 * norm_fact * (vv1[x][y] * fxyg[i][k])) / (D*D);
					AnmImg[G][x][y] = AnmImg[G][x][y] - (4.0 * norm_fact * (vv2[x][y] * fxyg[i][k])) / (D*D);
					AnmReal[B][x][y] = AnmReal[B][x][y] + (4.0 * norm_fact * (vv1[x][y] * fxyb[i][k])) / (D*D);
					AnmImg[B][x][y] = AnmImg[B][x][y] - (4.0 * norm_fact * (vv2[x][y] * fxyb[i][k])) / (D*D);
				}//closing of k loop
			}//closing of i loop
			//Calculating the magnitude of transform
			Anm[R][x][y] = sqrt(AnmReal[R][x][y] * AnmReal[R][x][y] + AnmImg[R][x][y] * AnmImg[R][x][y]);
			Anm[G][x][y] = sqrt(AnmReal[G][x][y] * AnmReal[G][x][y] + AnmImg[G][x][y] * AnmImg[G][x][y]);
			Anm[B][x][y] = sqrt(AnmReal[B][x][y] * AnmReal[B][x][y] + AnmImg[B][x][y] * AnmImg[B][x][y]);
			
			if (Anm[R][x][y] <= 0 || Anm[G][x][y] <= 0 || Anm[B][x][y] <= 0)
			{
				cout << "magnitude is 0 or negative:";
				_getch();	exit(1);
			}//closing of if block
			//Writing features to the file
			fprintf(fe_file1, "%d:%e\t", icount, Anm[R][x][y]);	icount++;
			fprintf(fe_file1, "%d:%e\t", icount, Anm[G][x][y]); icount++;
			fprintf(fe_file1, "%d:%e\t", icount, Anm[B][x][y]); icount++;
		}//closing of q loop
	}//closing of p loop
	return icount;
}

/*****************************************************************************************************************/
/* Y.N. Li, Quaternion polar harmonic transforms for color images, IEEE Signal Process. Lett. 20 (2013) 803–806. */
/* doi:10.1109/LSP.2013.2267775.                                                                                 */
/*****************************************************************************************************************/
int RQPCT(int imagenum, int nmax, double AnmReal[3][2 * NMAX + 1][2 * NMAX + 1], double AnmImg[3][2 * NMAX + 1][2 * NMAX + 1])
{
	int icount, n, k, m, x, y;
	double magn, Arnm[NMAX + 1][NMAX + 1], Brnm[NMAX + 1][NMAX + 1], Crnm[NMAX + 1][NMAX + 1], Drnm[NMAX + 1][NMAX + 1];
	double A1, B1, C1, D1, sqrt3;

	//Initialization section 
	icount = 1;
	sqrt3 = 1.0 / sqrt((double)3.0);

	// To compute the QPHT of a given color image
	for (n = 0; n <= nmax; n++)
	{
		x = n;
		for (m = 0; m <= nmax; m++)
		{
			y = m;
			Arnm[x][y] = (AnmImg[R][x][y] + AnmImg[G][x][y] + AnmImg[B][x][y])*(-sqrt3);
			Brnm[x][y] = ((AnmImg[G][x][y] - AnmImg[B][x][y])*sqrt3) + AnmReal[R][x][y];
			Crnm[x][y] = ((AnmImg[B][x][y] - AnmImg[R][x][y])*sqrt3) + AnmReal[G][x][y];
			Drnm[x][y] = ((AnmImg[R][x][y] - AnmImg[G][x][y])*sqrt3) + AnmReal[B][x][y];
			magn = sqrt((Arnm[x][y] * Arnm[x][y]) + (Brnm[x][y] * Brnm[x][y]) + (Crnm[x][y] * Crnm[x][y]) + (Drnm[x][y] * Drnm[x][y]));
			fprintf(fe_file2, "%d:%e\t", icount, magn); //Writing magnitude features to file
			icount++;
		}// Closing of m loop
	}// Closing of n loop
	icount = 1;
	// To compute the rotation invariants of QPHT
	for (n = 0; n <= nmax; n++)
	{
		for (k = 0; k <= n; k++)
		{
			for (m = 0; m <= nmax; m++)
			{
				A1 = -1 * (Arnm[n][m] * Arnm[k][m] + Brnm[n][m] * Brnm[k][m] + Crnm[n][m] * Crnm[k][m] + Drnm[n][m] * Drnm[k][m]);
				B1 = -1 * (-Arnm[n][m] * Brnm[k][m] + Brnm[n][m] * Arnm[k][m] - Crnm[n][m] * Drnm[k][m] + Drnm[n][m] * Crnm[k][m]);
				C1 = -1 * (-Arnm[n][m] * Crnm[k][m] + Brnm[n][m] * Drnm[k][m] + Crnm[n][m] * Arnm[k][m] - Drnm[n][m] * Brnm[k][m]);
				D1 = -1 * (-Arnm[n][m] * Drnm[k][m] - Brnm[n][m] * Crnm[k][m] + Crnm[n][m] * Brnm[k][m] + Drnm[n][m] * Arnm[k][m]);
				//Normalizing the feature by maintaining their signs
				if (A1 < 0)
					A1 = sqrt(fabs(A1))*-1.0;
				else
					A1 = sqrt(fabs(A1));
				if (B1 < 0)
					B1 = sqrt(fabs(B1))*-1.0;
				else
					B1 = sqrt(fabs(B1));
				if (C1 < 0)
					C1 = sqrt(fabs(C1))*-1.0;
				else
					C1 = sqrt(fabs(C1));
				if (D1 < 0)
					D1 = sqrt(fabs(D1))*-1.0;
				else
					D1 = sqrt(fabs(D1));
				//Writing features to the file
				fprintf(fe_file3, "%d:%e\t", icount, A1);	icount++;
				fprintf(fe_file3, "%d:%e\t", icount, B1);	icount++;
				fprintf(fe_file3, "%d:%e\t", icount, C1);	icount++;
				fprintf(fe_file3, "%d:%e\t", icount, D1);	icount++;
			}// Closing of m loop
		}// Closing of k loop
	}// // Closing of n loop
	return icount;
}
//Entery point
void main(int argc, char *argv[])
{
	if (argc != 4)
	{
		printf("Wrong number of arguments.\n First argument must be the name of .exe file"
			"\n Second argument must be the name of database file without .txt"
			"\n Third agrument must be 1 for outer circle and 0 for inner circle"
			"\n Fourth argument must be the order of transform, i.e., 7. ");
		_getch(); exit(1);
	}
	FILE *infile;
	char fn[100], fn1[100], buffer[50], filename[100], pht_f[100], qpht_f[100], rqpht_f[100];
	int db, inc = 0, fxyr[IMSIZE][IMSIZE] = { 0 }, fxyg[IMSIZE][IMSIZE] = { 0 }, fxyb[IMSIZE][IMSIZE] = { 0 },
		OC, no_of_images;
	int pmax, totalimages, totalclasses, imagesperclass = 0, i, j, k;
	double Anm[3][2 * NMAX + 1][2 * NMAX + 1] = { 0 }, AnmReal[3][2 * NMAX + 1][2 * NMAX + 1] = { 0 }, AnmImg[3][2 * NMAX + 1][2 * NMAX + 1] = { 0 },
		f_val = 0, sum = 0;
	int total_M1 = 0, total_M2 = 0, total_M3 = 0, total_M4 = 0, total_M5 = 0, flag = 1;

	OC = atoi(argv[2]);
	pmax = atoi(argv[3]);

	printf("\nComputing the moments/transforms at the order: %d\n", pmax);

	strcpy(fn1, argv[1]);
	strcpy(fn, fn1);
	strcpy(pht_f, fn);
	strcpy(qpht_f, fn);
	strcpy(rqpht_f, fn);
	strcat(fn, ".txt");
	if (OC == 1)
	{
		strcat(pht_f, "_MPHT_OC_");
		strcat(qpht_f, "_QPHT_OC_");
		strcat(rqpht_f, "_RQPHT_OC_");
	}
	else
	{
		strcat(pht_f, "_MPHT_IC_");
		strcat(qpht_f, "_QPHT_IC_");
		strcat(rqpht_f, "_RQPHT_IC_");
	}
	sprintf(buffer, "%d", pmax);
	strcat(pht_f, buffer);
	strcat(qpht_f, buffer);
	strcat(rqpht_f, buffer);

	//Opening database file
	if ((infile = fopen(fn, "r")) == NULL)
	{
		printf("can't open file - %s - either the filename is wrong or the"
			" file does not exist\n",fn);
		_getch();
		exit(1);
	}

	fscanf(infile, "%d%d%d", &totalimages, &totalclasses, &imagesperclass);
	total_im = totalimages; total_classes = totalclasses; im_per_class = imagesperclass;
	printf("Total no. of images=%d, No. of classes=%d, No. of images in each class=%d\n", total_im, total_classes, im_per_class);  //_getch();

	no_of_images = total_im;
	//MPHT feature file
	if ((fe_file1 = fopen(pht_f, "w")) == NULL)
	{
		printf("Fail to create PHT feature file.\n");
		_getch();
		exit(0);
	}
	//QPHT feature file
	if ((fe_file2 = fopen(qpht_f, "w")) == NULL)
	{
		printf("Fail to create QPHT feature file.\n");
		_getch();
		exit(0);
	}
	//RQPHTI feature file
	if ((fe_file3 = fopen(rqpht_f, "w")) == NULL)
	{
		printf("Fail to create RQPHT feature file.\n");
		_getch();
		exit(0);
	}

	//**************** Loop for computing features for all images ****************//

	for (j = 0; j < no_of_images; j++)
	{
		fscanf(infile, "%d%s%d%d", &num, filename, &classno, &imageno);
		printf("image no.=%d, filename=%s,  class no.=%d, imageno=%d\n", num, filename, classno, imageno);

		fscanf(infile, "%d%d", &M, &N);
		printf("Processing image no=%d, M=%d, N=%d\n", num, M, N);
		label[j] = classno;

		for (k = 0; k<M; k++)
		{
			for (i = 0; i<N; i++)
			{
				fscanf(infile, "%d", &db);
				fxyr[k][i] = db;
				fscanf(infile, "%d", &db);
				fxyg[k][i] = db;
				fscanf(infile, "%d", &db);
				fxyb[k][i] = db;
			}
		}
		fprintf(fe_file1, "%d\t", classno);
		fprintf(fe_file2, "%d\t", classno);
		fprintf(fe_file3, "%d\t", classno);

		total_M1 = PCT(j, pmax, fxyr, fxyg, fxyb, Anm, AnmReal, AnmImg, OC);
		total_M2 = RQPCT(j, pmax, AnmReal, AnmImg);
		// To insert a new line for next image features
		fprintf(fe_file1, "\n");
		fprintf(fe_file2, "\n");
		fprintf(fe_file3, "\n");
	}
	//After the complete excution of program, it beeps to get the attention
	for (i = 400; i<800; i = i + 50) { Beep(i, 1000); }
	fclose(infile);
	//_getch();
}