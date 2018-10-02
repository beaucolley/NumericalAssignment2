/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : <INSERT STUDENT ID HERE>
 *   Name        : <INSERT STUDENT NAME HERE>
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"
#include <assert.h>

//Constants
#define BUFFER 50

//DEBUG SWITCHES
#define DEBUG_NEWTON_RAPHSON 1
#define DEBUG_SHOCKWAVE 1

void shockwave(const char* q2_file)
{
    printf("shockwave - called\n");

    float B_l, B_u, M, T, G;

    //Open File
    FILE *data_in;
	data_in = fopen(q2_file, "r");
	assert(data_in!=NULL);

	//Part 2.3

	//Part A
	//Check Headings
	char headings[BUFFER];
	fscanf(data_in, "%s",headings);
	if(strcmp(headings,"M,theta,beta_l,beta_u,gamma")!=0){
		printf("Incorrect data format - expecting <M,theta,beta_l,beta_u,gamma>\n");
		printf("Data is in %s\n",headings);
		return;
	}

	//Retrieve Parameters
	fscanf(data_in,"%f,%f,%f,%f,%f",&M,&T,&B_l,&B_u,&G);

	//Solve
	shockSolution partA;
	partA.strong = NewtonRaphson(M,B_u,G,T);
	partA.weak = NewtonRaphson(M,B_l,G,T);

	//Part B
	//Check Headings
	fscanf(data_in, "%s",headings);
	if(strcmp(headings,"M")!=0){
		printf("Incorrect data format - expecting <M>\n");
		printf("Data is in %s\n",headings);
		return;
	}

	//Retrieve Parameters
	fscanf(data_in,"%f",&M);

	//Set Parameters
	Array partB;
	initArray(&partB);
	T = 0; //Starting Theta

	thetaSweep(&partB,M,T,B_u,B_l,G);

	if(DEBUG_SHOCKWAVE){
		for(int i=0;i<partB.used;i++){
			printf("%f,%f,%f,%f\n",partB.array[i].M,partB.array[i].T,
					partB.array[i].weak,partB.array[i].strong);
		}
	}

	freeArray(&partB);

	//PartC
	//Check Headings
	fscanf(data_in, "%s",headings);
	if(strcmp(headings,"M")!=0){
		printf("Incorrect data format - expecting <M>\n");
		printf("Data is in %s\n",headings);
		return;
	}

	float mVals[100];
	int i = 0;
	float tempM;
	while(fscanf(data_in,"%f",&tempM) > 0){
		mVals[i] = tempM;
		i++;
		printf("%f\n",tempM);
	}

	Array partC;
	initArray(&partC);

	for(int j=0;j<i;j++){
		thetaSweep(&partC,mVals[j],T,B_u,B_l,G);
	}

	if(DEBUG_SHOCKWAVE){
		for(int i=0;i<partC.used;i++){
			printf("%f,%f,%f,%f\n",partC.array[i].M,partC.array[i].T,
					partC.array[i].weak,partC.array[i].strong);
		}
	}


	fclose(data_in);

	int noSolution = 0;


}

void linalgbsys(const char* q4_file)
{
    printf("linalgbsys() - IMPLEMENT ME!\n");
    exit(EXIT_FAILURE);
}

void interp(const char* q5_file, const double xo)
{
    printf("interp() - IMPLEMENT ME!\n");
    exit(EXIT_FAILURE);
}

void waveeqn(const char* q6_file)
{
    printf("heateqn() - IMPLEMENT ME!\n");
    exit(EXIT_FAILURE);
}

float f(float M, float B, float T)
{
	float G = 1.4;
	float X = pow(M,2)*pow(sin(B),2)-1;
	float Y = pow(M,2)*(G+cos(2*B))+2;

    return (2/tan(B))*(X/Y)-tan(T);
}

float df (float M, float B, float T, float G)
{
	return (2*(2*pow(M,2)*pow(cos(B),2)*(2+pow(M,2)*(G+cos(2*B)))-
			(2+G*pow(M,2)+pow(M,2)*cos(2*B))*(pow(M,2)-pow(1/sin(B),2))+
			2*pow(M,2)*(cos(B)/sin(B))*(-1+pow(M,2)*pow(sin(B),2))*sin(2*B)))/
			pow((2+pow(M,2)*(G+cos(2*B))),2);
}
float DegtoRad(float degrees){
	return degrees*(M_PI/180);
}
float RadtoDeg(float rad){
	return rad*(180/M_PI);
}
float NewtonRaphson(float M, float B, float G, float T){

	int itr;
	int maxmitr = 100;
	float h, x1;
	float allerr = 0.01;

	T = DegtoRad(T);
	B = DegtoRad(B);

	for (itr=1; itr<=maxmitr; itr++)
	{
		h=f(M,B,T)/df(M,B,T,G);

		x1=B-h;
		if(DEBUG_NEWTON_RAPHSON){
			printf(" At Iteration no. %3d, x = %9.6f h = %f\n", itr, x1, h);
		}

		if (fabs(h) < allerr)
		{
			if(DEBUG_NEWTON_RAPHSON){
				printf("After %3d iterations, root = %8.6f\n", itr, RadtoDeg(x1));
			}
			return RadtoDeg(x1);
		}
		B=x1;
	}
	if(DEBUG_NEWTON_RAPHSON){
		printf(" The required solution does not converge or iterations are insufficient\n");
	}
	return -1;
}

void initArray(Array *a)
{
    // Allocate initial space
    a->array = (shockSolution *)malloc(sizeof(shockSolution));

    a->used = 0; 	// no elements used
    a->size = 1; 	// set initial size

}

// Add element to array
void insertArray(Array *a, shockSolution element)
{
    if (a->used == a->size)
    {
        a->size *= 2;
        a->array = (shockSolution *)realloc(a->array, a->size * sizeof(shockSolution));
    }

    // Add element
    a->array[a->used] = element;

    a->used++;
}

void freeArray(Array *a)
{

	// Now free the array
    free(a->array);
    a->array = NULL;

    a->used = 0;
    a->size = 0;
}

void thetaSweep(Array* results, float M, float T, float B_u, float B_l, float G){
	int solution = 1;
	while(solution){
		shockSolution temp;
		temp.strong = NewtonRaphson(M,B_u,G,T);
		temp.weak = NewtonRaphson(M,B_l,G,T);
		temp.M = M;
		temp.T = T;

		if(temp.strong != -1 && temp.weak != -1){
			insertArray(results,temp);
		}else{
			solution = 0;
		}
		T++;
	}
	return;
}

float* readMvals(FILE* data_in){

}

