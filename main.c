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
#define TASK3 1
#define TASK4 1
#define TASK5 1

int main(int argc, char *argv[]) {
	
	/* TODO: Parse Command Line Arguments
	DONOT explicitly set arguments to filenames */
	char* q2_file = argv[1];
	char* q4_file = argv[2];
	char* q5_file = argv[3];
	char* q6_file = argv[4];

	double xo;

	/* TODO: Add timing for each task and output running time in ms */
    
	/* Question 2 */
	if (TASK2) {
		printf("Shockwave - Started\n\n");
		shockwave(q2_file);
		printf("Shockwave - Finished\n\n");
	}
	
	/* Question 4 */
	if (TASK3) {
		printf("Linear Systems - Started\n\n");
		linalgbsys(q4_file);
		printf("Linear Systems - Finished\n\n");
	}
	
	/* Question 5 */
	if (TASK4) {
		printf("Interpolation - Started\n\n");
		interp(q5_file,xo);
		printf("Interpolation - Finished\n\n");
	}
	
	/* Question 6 */
	if (TASK5) {
		printf("WaveEquation - Started\n\n");
		waveeqn(q6_file);
		printf("WaveEquation - Finished\n\n");
	}

	printf("All tasks Completed\n");
    
	return (EXIT_SUCCESS);
}
