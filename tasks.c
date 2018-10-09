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
#define PI 3.1415

//DEBUG SWITCHES
#define DEBUG_NEWTON_RAPHSON 0
#define DEBUG_SHOCKWAVE 0
#define DEBUG_LINALSYS 0
#define DEBUG_INTERP 1

void shockwave(const char* q2_file)
{
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

	printf("%f,%f,%f,%f\n",partA.M,partA.T,partA.weak,partA.strong);

	if(DEBUG_SHOCKWAVE){
		printf("Part A Complete\n");
	}

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
	shockArray partB;
	initArray_shock(&partB);
	T = 0; //Starting Theta

	//Solve varying Theta Values
	thetaSweep(&partB,M,T,B_u,B_l,G);

	if(DEBUG_SHOCKWAVE){
		for(int i=0;i<partB.used;i++){
			printf("%f,%f,%f,%f\n",partB.array[i].M,partB.array[i].T,
					partB.array[i].weak,partB.array[i].strong);
		}
	}


	freeArray_shock(&partB);

	if(DEBUG_SHOCKWAVE){
		printf("Part B Complete\n");
	}

	//PartC
	//Check Headings
	fscanf(data_in, "%s",headings);
	if(strcmp(headings,"M")!=0){
		printf("Incorrect data format - expecting <M>\n");
		printf("Data is in %s\n",headings);
		return;
	}

	//Create Dynamic Array to store M values
	floatArray mVals;
	initArray_float(&mVals);
	float temp;
	while(fscanf(data_in,"%f",&temp) > 0){
		insertArray_float(&mVals,temp);
	}

	fclose(data_in);

	//Create Array to store shock solutions
	shockArray partC;
	initArray_shock(&partC);

	for(int i=0;i<mVals.used;i++){
		thetaSweep(&partC,mVals.array[i],T,B_u,B_l,G);
	}

	if(DEBUG_SHOCKWAVE){
		for(int i=0;i<partC.used;i++){
			printf("%f,%f,%f,%f\n",partC.array[i].M,partC.array[i].T,
					partC.array[i].weak,partC.array[i].strong);
		}
	}

	//Open output File
	FILE* data_out;
	data_out = fopen("out_shock.csv","w+");

	//Print to file
	for(int i=0;i<partC.used;i++){
		fprintf(data_out,"%f,%f,%f,%f\n",partC.array[i].M,partC.array[i].T,
			partC.array[i].weak,partC.array[i].strong);
	}

	freeArray_shock(&partC);
	freeArray_float(&mVals);

	fclose(data_out);


	if(DEBUG_SHOCKWAVE){
		printf("Part C Complete\n");
	}


}

void linalgbsys(const char* q4_file)
{

    //Open File
	FILE *data_in;
	data_in = fopen(q4_file, "r");
	assert(data_in!=NULL);
	char headings[BUFFER];

	//Initialise Arrays
	floatArray A,B,C,Q,X;
	initArray_float(&A);
	initArray_float(&B);
	initArray_float(&C);
	initArray_float(&X);
	initArray_float(&Q);

	fscanf(data_in, "%s",headings);
	if(strcmp(headings,"a,b,c,q")!=0){
		printf("Incorrect data format - expecting <a,b,c,q>\n");
		printf("Data is in %s\n",headings);
		return;
	}
	float aTemp, bTemp, cTemp, qTemp;

	while(fscanf(data_in,"%f,%f,%f,%f",&aTemp,&bTemp,&cTemp,&qTemp) > 0){
		insertArray_float(&A, aTemp);
		insertArray_float(&B, bTemp);
		insertArray_float(&C, cTemp);
		insertArray_float(&Q, qTemp);
	}

	if(DEBUG_LINALSYS){
		for(int i = 0;i<A.used;i++){
			printf("%f,%f,%f,%f\n",A.array[i],B.array[i],C.array[i],Q.array[i]);
		}
	}

	thomasAlgorithm(&A,&B,&C,&Q,&X);

	if(DEBUG_LINALSYS){
		for(int i = 0;i<A.used;i++){
			printf("%f\n",X.array[i]);
		}
	}

	fclose(data_in);

	//Free Memory
	freeArray_float(&A);
	freeArray_float(&B);
	freeArray_float(&C);
	freeArray_float(&Q);
	freeArray_float(&X);

}

void interp(const char* q5_file, const double xo)
{
    //Open File
	FILE *data_in;
	data_in = fopen(q5_file, "r");
	assert(data_in!=NULL);
	char headings[BUFFER];

	//Read Headings
	fscanf(data_in, "%s",headings);
	if(strcmp(headings,"x,f(x)")!=0){
		printf("Incorrect data format - expecting <x,f(x)\n");
		printf("Data is in %s\n",headings);
		return;
	}

	//Define Variables
	floatArray x, y;
	float xTemp, yTemp;
	float target = 5;

	//Initialise Arrays
	initArray_float(&x);
	initArray_float(&y);

	//Scan in Data
	while(fscanf(data_in,"%f,%f",&xTemp,&yTemp) > 0){
		insertArray_float(&x,xTemp);
		insertArray_float(&y,yTemp);
	}

	//DEBUG Print Data
	if (DEBUG_INTERP) {
		for(int i = 0; i<x.used; i++){
			printf("%f, %f\n",x.array[i],y.array[i]);
		}
	}

	float lagrangeEstimate = lagrangeInterp(x,y,target);
	printf("\nlagrange\n%f\n",lagrangeEstimate);

	float cubicSlineEstimate = cubicSplineInterp(x,y,target);
	printf("\ncubic\n%f\n",cubicSlineEstimate);



	//Free Arrarys
	freeArray_float(&x);
	freeArray_float(&y);

}

void waveeqn(const char* q6_file)
{
	//Open File
	FILE *data_in;
	data_in = fopen(q6_file, "r");
	assert(data_in!=NULL);
	char headings[BUFFER];

	//Read Headings
	fscanf(data_in, "%s",headings);
	if(strcmp(headings,"c,Nx,CFL,out_iter")!=0){
		printf("Incorrect data format - expecting <c.Nx,CFL,out_iter>\n");
		printf("Data is in %s\n",headings);
		return;
	}

	//Define Variables for input parameters
	float c, CFL, out_iter;
	int Nx;
	fscanf(data_in,"%f,%d,%f,%f",&c,&Nx,&CFL,&out_iter);

	fclose(data_in);

	//Open output File
	FILE *data_out;
	data_out = fopen("out_waveeqn.csv","w+");

	//Define internal parameters
	float dx = 1/(float)Nx;
	float dt = 0.01;
	float T = 10;
	float fn[Nx+1];
	float fn1[Nx+1];

	//Initialise Array
	initaliseFn(fn,dx);

	//Save Data
	exportArray(data_out,fn, Nx);

	//Perform Upwinding
	float t = 0;

	while(t<T){
		eulerUpwinding(fn,fn1,c,dt,dx,Nx);
		exportArray(data_out,fn1,Nx);
		updateArray(fn,fn1,Nx);
		t += dt;
	}


	fclose(data_out);

}

void updateArray(float fn[], float fn1[], int Nx){
	for (int i = 0; i <= Nx; ++i) {
		fn[i] = fn1[i];
	}
}

void eulerUpwinding(float fn[], float fn1[], float c, float dt, float dx, int Nx){
	//Solve with boundry condition
	fn1[0] = fn[0] -c*(dt/dx)*(fn[0]- fn[Nx]);

	//Solve for remaining points
	for(int i = 1; i<=Nx; i++){
		fn1[i] = fn[i] -c*(dt/dx)*(fn[i]- fn[i-1]);
	}

}

void exportArray(FILE* data_out, float fn[], int Nx){
	fprintf(data_out,"%f",fn[0]);
	for (int i = 1; i <= Nx; ++i) {
		fprintf(data_out,",%f",fn[i]);
	}
	fprintf(data_out,"\n");
}

void initaliseFn(float fn[], float dx){
	float x = 0;
	int i = 0;
	while(x<=1){
		if(x<0.125)
			fn[i] = 0;
		else if(x<=0.375)
			fn[i] = 0.5*(1-cos(8*PI*(x-0.125)));
		else if(x<=1)
			fn[i] = 0;
		else
			printf("Initialise Out of Bounds\n");
		i++;
		x += dx;
	}
}

void tridiagonalCubicSplineGen(int n, float h[], float triMatrix[][n], float y[]){
    int i;
    for(i=0;i<n-1;i++){
        triMatrix[i][i]=2*(h[i]+h[i+1]);
    }
    for(i=0;i<n-2;i++){
        triMatrix[i][i+1]=h[i+1];
        triMatrix[i+1][i]=h[i+1];
    }
    for(i=1;i<n;i++){
        triMatrix[i-1][n-1]=(y[i+1]-y[i])*6/(float)h[i]-(y[i]-y[i-1])*6/(float)h[i-1];
    }
}

float cubicSplineInterp(floatArray x, floatArray y, float target){
	int n = x.used - 1;
	float a[n]; //array to store the ai's
	float b[n]; //array to store the bi's
	float c[n]; //array to store the ci's
	float d[n]; //array to store the di's
	float h[n];
	float sig[n+1]; //array to store Si's
	float sigTemp[n-1]; //array to store the Si's except S0 and Sn

	//Nature Spline
	sig[0]=0;
	sig[n]=0;

	//Calculate Interval Values
	intervalValues(n,h,x);

	//Matrix to store the tridiagonal system of equations that will solve for Si's
	float tri[n-1][n];
	tridiagonalCubicSplineGen(n,h,tri,y.array); //to initialize tri[n-1][n]

	//DEBUG Print Matrix
	if (DEBUG_INTERP) {
		printMatrix(n-1,n,tri);
	}

	//Wrapper function to apply Thomas Algorithm
	thomasWrapper(n,tri,sigTemp);

	//Merge Si's
	for(int i=1;i<n;i++){
		sig[i]=sigTemp[i-1];
	}

	//DEBUG Print Si's
	if (DEBUG_INTERP) {
		for(int i=0;i<n+1;i++){
			printf("\nSig[%d] = %lf\n",i,sig[i]);
		}
	}

	//Calculate Spline Coefficients
	cSCoeffCalc(n,h,sig,y.array,a,b,c,d);

	//DEBUG Print Spline Equations
	if (DEBUG_INTERP) {
		printf("\nThe equations of cubic interpolation polynomials between the successive intervals are:\n\n");
		for(int i=0;i<n;i++){
			printf("P%d(x) b/w [%f,%f] = %f*(x-%f)^3+%f*(x-%f)^2+%f*(x-%f)+%f\n",i,x.array[i],x.array[i+1],a[i],x.array[i],b[i],x.array[i],c[i],x.array[i],d[i]);
		}
	}

	return interpolate(n,a,b,c,d,h,x.array,target);

}

float interpolate(int n, float a[], float b[], float c[], float d[], float h[], float x[], float target){

	//Determine which spline region target resides in
	int i = 0;

	//Check Target is within interpolation range
	if(target < x[0] || target > x[n]){
		printf("Point is outside of spline interpolation range\n");
		return 0;
	}

	while(target > h[i]){
		i++;
	}

	return (a[i] + b[i]*(target-x[i]) + c[i]*pow((target-x[i]),2) + d[i]*pow((target - x[i]),3));

}

void cSCoeffCalc(int n, float h[], float sig[], float y[], float a[], float b[], float c[], float d[]){

	int i;
    for(i=0;i<n;i++){
        d[i]=y[i];
        b[i]=sig[i]/2.0;
        a[i]=(sig[i+1]-sig[i])/(h[i]*6.0);
        c[i]=(y[i+1]-y[i])/h[i]-h[i]*(2*sig[i]+sig[i+1])/6.0;
//    	a[i] = y[i];
//    	c[i] = sig[i];


    }

    if (DEBUG_INTERP) {
		printf("a[]\t\tb[]\t\tc[]\t\td[]\t\ty[]\t\th[]\n");
		for (int i = 0; i < n; ++i) {
			printf("%f\t%f\t%f\t%f\t%f\t%f\n",a[i],b[i],c[i],d[i],y[i],h[i]);
		}
	}
}

void thomasWrapper(int n, float triMatrix[][n], float sigTemp[]){
	//Initialise Arrays
	floatArray A,B,C,Q,X;
	initArray_float(&A);
	initArray_float(&B);
	initArray_float(&C);
	initArray_float(&X);
	initArray_float(&Q);

	for(int i=0; i<n-1; i++){
		insertArray_float(&A,triMatrix[i][i]);
	}


	for(int i=0; i<n-2; i++){
		insertArray_float(&B,triMatrix[i][i+1]);
	}

	insertArray_float(&B,0);
	insertArray_float(&C,0);

	for(int i=0; i<n-2; i++){
		insertArray_float(&C,triMatrix[i+1][i]);
	}

	for(int i=0; i<n-1; i++){
		insertArray_float(&Q,triMatrix[i][n-1]);
	}

	thomasAlgorithm(&A,&B,&C,&Q,&X);

	for(int i=0; i<n-1; i++){
		sigTemp[i] = X.array[i];
	}

}

void thomasAlgorithm(floatArray *A,floatArray *B,floatArray *C,floatArray *Q,floatArray *X){

	//Initialise X
	for(int i=0;i<A->used;i++){
		insertArray_float(X,0.0);
	}

	//Forward Sweep
	for(int i=1;i<A->used;i++){
		A->array[i] = A->array[i] - (C->array[i]*B->array[i-1])/A->array[i-1];
		Q->array[i] = Q->array[i] - (C->array[i]*Q->array[i-1])/A->array[i-1];
	}


	if (DEBUG_INTERP) {
		printf("[Q,A]\n");
		for(int i=0;i<A->used;i++){
			printf("%f,%f\n",A->array[i],Q->array[i]);
		}
	}

	//Solve Xn
	X->array[X->used-1] = Q->array[Q->used-1]/A->array[A->used-1];

	//Reverse Sweep
	for(int i=X->used-2;i>=0;i--){
		X->array[i] = (Q->array[i]-B->array[i]*X->array[i+1])/A->array[i];
	}


}

void intervalValues(int n, float h[], floatArray x){

	for(int i=0;i<n;i++){
		h[i]=x.array[i+1]-x.array[i];
	}

	if (DEBUG_INTERP) {
		printf("\nh[]\n");
		for(int i = 0; i < n; i++) {
			printf("%f\n",h[i]);
		}
	}
}

void printMatrix(int m, int n, float matrix[m][n]){
    int i,j;
    printf("\nTriBand Matrix\n");
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            printf("%f\t",matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

float lagrangeInterp(floatArray x, floatArray y, float target){

	float estimate = 0;

	for(int i=0; i<x.used; i++)
	{
		estimate += lagrangeMultiplier(x,y,i,target)*y.array[i];
	}
	return estimate;
}

float lagrangeMultiplier(floatArray x_points,floatArray y_points,int i, float x){

	float output =1;

	for(int j=0; j<x_points.used; j++){
		if(i!=j){
			output *= (x-x_points.array[j])/(x_points.array[i]-x_points.array[j]);
		}
	}

	return output;
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

void initArray_shock(shockArray *a)
{
    // Allocate initial space
    a->array = (shockSolution *)malloc(sizeof(shockSolution));

    a->used = 0; 	// no elements used
    a->size = 1; 	// set initial size

}

// Add element to array
void insertArray_shock(shockArray *a, shockSolution element)
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

void initArray_float(floatArray *a)
{
    // Allocate initial space
    a->array = (float*)malloc(sizeof(float));

    a->used = 0; 	// no elements used
    a->size = 1; 	// set initial size

}

// Add element to array
void insertArray_float(floatArray *a, float element)
{
    if (a->used == a->size)
    {
        a->size *= 2;
        a->array = (float *)realloc(a->array, a->size * sizeof(float));
    }

    // Add element
    a->array[a->used] = element;

    a->used++;
}

void freeArray_shock(shockArray *a)
{

	// Now free the array
    free(a->array);
    a->array = NULL;

    a->used = 0;
    a->size = 0;
}

void freeArray_float(floatArray *a)
{

	// Now free the array
    free(a->array);
    a->array = NULL;

    a->used = 0;
    a->size = 0;
}

void thetaSweep(shockArray* results, float M, float T, float B_u, float B_l, float G){
	int solution = 1;
	while(solution){
		shockSolution temp;
		temp.strong = NewtonRaphson(M,B_u,G,T);
		temp.weak = NewtonRaphson(M,B_l,G,T);
		temp.M = M;
		temp.T = T;

		if(temp.strong != -1 && temp.weak != -1){
			insertArray_shock(results,temp);
		}else{
			solution = 0;
		}
		T++;
	}
	return;
}



