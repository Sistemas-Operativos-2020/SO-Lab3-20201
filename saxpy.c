/**
 * @defgroup   SAXPY saxpy
 *
 * @brief      This file implements an iterative saxpy operation
 * 
 * @param[in] <-p> {vector size} 
 * @param[in] <-s> {seed}
 * @param[in] <-n> {number of threads to create} 
 * @param[in] <-i> {maximum itertions} 
 *
 * @author     Danny Munera
 * @date       2020
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <pthread.h>

//  {}  \n  []    ||

// Variables to obtain command line parameters
unsigned int seed = 1;
int p = 10000000;
int n_threads = 2;
int max_iters = 1000;	
// Variables to perform SAXPY operation
double* X;
double a;
double* Y;
double* Y_avgs;
int i, it;

struct rangos{
	int inicio;
	int cant;
};

void* do_saxpy(void* args){
	struct rangos *p = (struct rangos *) args;
	//SAXPY iterative SAXPY mfunction
	double *Y_sum = (double*) malloc(sizeof(double) * max_iters);
	int i;
	int fin = p -> inicio + p -> cant;
	for(int iter = 0; iter < max_iters; iter++){
		for(i = p -> inicio; i < fin; i++){
			Y[i] = Y[i] + a * X[i];
			Y_sum[iter] += Y[i];
		}
	}
	return Y_sum;
}

int main(int argc, char* argv[]){
	// Variables to get execution time
	struct timeval t_start, t_end;
	double exec_time;

	// Getting input values
	int opt;
	while((opt = getopt(argc, argv, ":p:s:n:i:")) != -1){  
		switch(opt){  
			case 'p':  
			printf("vector size: %s\n", optarg);
			p = strtol(optarg, NULL, 10);
			assert(p > 0 && p <= 2147483647);
			break;  
			case 's':  
			printf("seed: %s\n", optarg);
			seed = strtol(optarg, NULL, 10);
			break;
			case 'n':  
			printf("threads number: %s\n", optarg);
			n_threads = strtol(optarg, NULL, 10);
			break;  
			case 'i':  
			printf("max. iterations: %s\n", optarg);
			max_iters = strtol(optarg, NULL, 10);
			break;  
			case ':':  
			printf("option -%c needs a value\n", optopt);  
			break;  
			case '?':  
			fprintf(stderr, "Usage: %s [-p <vector size>] [-s <seed>] [-n <threads number>]\n", argv[0]);
			exit(EXIT_FAILURE);
		}  
	}  
	srand(seed);

	printf("p = %d, seed = %d, n_threads = %d, max_iters = %d\n", \
	 p, seed, n_threads, max_iters);	

	// initializing data
	X = (double*) malloc(sizeof(double) * p);
	Y = (double*) malloc(sizeof(double) * p);
	Y_avgs = (double*) malloc(sizeof(double) * max_iters);

	for(i = 0; i < p; i++){
		X[i] = (double)rand() / RAND_MAX;
		Y[i] = (double)rand() / RAND_MAX;
	}
	for(i = 0; i < max_iters; i++){
		Y_avgs[i] = 0.0;
	}
	a = (double)rand() / RAND_MAX;

#ifdef DEBUG
	printf("vector X= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ",X[i]);
	}
	printf("%f ]\n",X[p-1]);

	printf("vector Y= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ", Y[i]);
	}
	printf("%f ]\n", Y[p-1]);

	printf("a= %f \n", a);	
#endif

	pthread_t hilos[n_threads];
	struct rangos args[n_threads];
	int cant_por_hilo = p / n_threads;
	int overflow = p % n_threads;
	int finAnt = 0;
	double* results[n_threads];

	for(int t = 0 ; t < n_threads ; t++){
		if(t < overflow)
			args[t].cant = cant_por_hilo + 1;
		else
			args[t].cant = cant_por_hilo;
		
		args[t].inicio = finAnt;
		finAnt += args[t].cant;
	}

	gettimeofday(&t_start, NULL);

	for(int t = 0 ; t < n_threads ; t++){
		pthread_create(&hilos[t], NULL, &do_saxpy, &args[t]);
	}

	for(int t = 0 ; t < n_threads ; t++){
		pthread_join(hilos[t], (void*) &results[t]);
	}

	for(i = 0 ; i < max_iters ; i++){
		for(int t = 0 ; t < n_threads ; t++){
			Y_avgs[i] += results[t][i];
		}
		Y_avgs[i] = Y_avgs[i] / p;
	}

	gettimeofday(&t_end, NULL);

#ifdef DEBUG
	printf("RES: final vector Y= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ", Y[i]);
	}
	printf("%f ]\n", Y[p-1]);
#endif
	
	// Computing execution time
	exec_time = (t_end.tv_sec - t_start.tv_sec) * 1000.0;  // sec to ms
	exec_time += (t_end.tv_usec - t_start.tv_usec) / 1000.0; // us to ms
	printf("Execution time: %f ms \n", exec_time);
	printf("Last 3 values of Y: %f, %f, %f \n", Y[p-3], Y[p-2], Y[p-1]);
	printf("Last 3 values of Y_avgs: %f, %f, %f \n", Y_avgs[max_iters-3], Y_avgs[max_iters-2], Y_avgs[max_iters-1]);
	return 0;
}	