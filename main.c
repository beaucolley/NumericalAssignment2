/***************************************************************************
 *
 *   File        : main.c
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

//TASK Switches
#define TASK2 1
#define TASK3 0
#define TASK4 0
#define TASK5 0

int main(int argc, char *argv[]) {
	

	char* q2_file = argv[1];
	char* q4_file = argv[2];
	char* q5_file = argv[3];
	char* q6_file = argv[5];

	float xo = atof(argv[4]);


	struct timeval start;
	struct timeval stop;

    
	/* Question 2 */
	if (TASK2) {
		gettimeofday(&start, NULL);
		shockwave(q2_file);
		gettimeofday(&stop, NULL);
	}
	double elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	if(TASK2){
			elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
			printf("Root Finding:  %.2f milliseconds\n", elapsed_ms);
			printf("\n");
	}
	
	/* Question 4 */
	if (TASK3) {
		gettimeofday(&start, NULL);
		linalgbsys(q4_file);
		gettimeofday(&stop, NULL);
		elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
		elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
		printf("Linear Algebraic Systems:  %.2f milliseconds\n", elapsed_ms);
		printf("\n");
	}
	
	/* Question 5 */
	if (TASK4) {
		gettimeofday(&start, NULL);
		interp(q5_file,xo);
		gettimeofday(&stop, NULL);
		elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
		elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
		printf("Interpolation:  %.2f milliseconds\n", elapsed_ms);
		printf("\n");
	}
	
	/* Question 6 */
	if (TASK5) {
		gettimeofday(&start, NULL);
		waveeqn(q6_file);
		gettimeofday(&stop, NULL);
		elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
		elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
		printf("Wave Equations:  %.2f milliseconds\n", elapsed_ms);
		printf("\n");
	}

	printf("All tasks Completed\n");
    
	return (EXIT_SUCCESS);
}
