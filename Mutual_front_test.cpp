#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <algorithm>
#include "Mutual_front_v3.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;				
				
int main(int argc, char *argv[])
{
	int dimension=atoi(argv[1]);
	int linear_size=atoi(argv[2]);
	int neighbor_width=atoi(argv[3]);
	int num_states=atoi(argv[4]);
	double mutualism_strength=atof(argv[5]);
	double p_swap=atof(argv[6]);
	int tfin=atoi(argv[7]);
	int trials=atoi(argv[8]);
	unsigned int refresh_time=(atoi(argv[9]))*60; //min
	int thread=atoi(argv[10]);
	
	char FILE_NAME [100];
	sprintf(FILE_NAME,"Mutual_front_%d_%d_%d_%g_%g_%d_%d.txt",linear_size,neighbor_width,num_states,mutualism_strength,p_swap,tfin,thread);

	double **aij;
	
	aij =  new double *[num_states];
	for (int i = 0; i < num_states; ++i) aij[i] = new double[num_states];
	
	for(int i=0;i<num_states;i++){
		for(int j=0;j<num_states;j++){
			if (i!=j) aij[i][j]=mutualism_strength;
			else aij[i][j]=0;
		}
	}

	System system(dimension,linear_size,neighbor_width,num_states,aij,p_swap,1);
	Visualize_evolution visualize(system,tfin,trials,refresh_time,thread,FILE_NAME);

	visualize.simulate();

	return EXIT_SUCCESS;

}
