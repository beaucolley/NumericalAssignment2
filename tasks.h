/***************************************************************************
 *
 *   File        : tasks.h
 *   Student Id  : <INSERT STUDENT ID HERE>
 *   Name        : <INSERT STUDENT NAME HERE>
 *
 ***************************************************************************/

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
}Array;

void shockwave(const char* q2_file);

void linalgbsys(const char* q4_file);

void interp(const char* q5_file, const double xo);

void waveeqn(const char* q6_file);

float f(float M, float B, float T);

float df (float M, float B, float T, float G);

float DegtoRad(float degrees);

float RadtoDeg(float rad);

float NewtonRaphson(float M, float B, float G, float T);

void initArray(Array *a);

void insertArray(Array *a, shockSolution element);

void freeArray(Array *a);

//void thetaSweep(Array *results, float M, float T, float B_u, float B_l, float G);
Array thetaSweep(float M, float T, float B_u, float B_l, float G);


#endif
