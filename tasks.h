/***************************************************************************
 *
 *   File        : tasks.h
 *   Student Id  : <INSERT STUDENT ID HERE>
 *   Name        : <INSERT STUDENT NAME HERE>
 *
 ***************************************************************************/
#include <stdio.h>
#ifndef TASKS_H

typedef struct{
	float strong;
	float weak;
	float M;
	float T;
}shockSolution;

typedef struct{
	shockSolution* array;
	int used;
	int size;
}shockArray;

typedef struct{
	float* array;
	int used;
	int size;
}floatArray;

void shockwave(const char* q2_file);

void linalgbsys(const char* q4_file);

void interp(const char* q5_file, float xo);

void waveeqn(const char* q6_file);

float f(float M, float B, float T);

float df (float M, float B, float T, float G);

float DegtoRad(float degrees);

float RadtoDeg(float rad);

float NewtonRaphson(float M, float B, float G, float T);

void initArray_shock(shockArray *a);

void insertArray_shock(shockArray *a, shockSolution element);

void freeArray_float(floatArray *a);

void freeArray_shock(shockArray *a);

void initArray_float(floatArray *a);

void insertArray_float(floatArray *a, float element);

void thetaSweep(shockArray *results, float M, float T, float B_u, float B_l, float G);

void thomasAlgorithm(floatArray *A,floatArray *B,floatArray *C,floatArray *Q,floatArray *X);

float lagrangeInterp(floatArray x, floatArray y, float target);

float lagrangeMultiplier(floatArray x_points,floatArray y_points,int i, float x);

float cubicSplineInterp(floatArray x, floatArray y, float target);

void tridiagonalCubicSplineGen(int n, float h[], float triMatrix[][n], float y[n+1]);

void intervalValues(int n, float h[], floatArray x);

void printMatrix(int m, int n, float matrix[m][n]);

void thomasWrapper(int n, float triMatrix[][n], float sigTemp[]);

void cSCoeffCalc(int n, float h[], float sig[], float y[], float a[], float b[], float c[], float d[]);

float interpolate(int n, float a[], float b[], float c[], float d[], float h[], float x[], float target);

void initaliseFn(float fn[], float dx);

void exportArray(FILE *data_out, float fn[], int Nx);

void eulerUpwinding(float fn[], float fn1[], float c, float dt, float dx, int Nx);

void updateArray(float fn[], float fn1[], int Nx);

void laxWendroff(float fn[], float fn1[], float c, float dt, float dx, int Nx);

#endif
