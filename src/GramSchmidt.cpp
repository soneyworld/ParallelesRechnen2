/*
 * GramSchmidt.cpp
 *
 *      Author: Till Lorentzen
 */
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include "config.h"
#include "Benchmark.h"
#include "GramSchmidt.h"

#define ORTHOGONAL
#define RANDOMMATRIX

#define MATRIXDIMENSION 1000
#define W(i,j) w[(i)*lda+(j)]
#define V(i,j) v[(i)*lda+(j)]

#ifdef OPENMP
#include <omp.h>
#endif

//int addition = 0;
//int multiplication = 0;
//int division = 0;

NUMBER dot(NUMBER *a, NUMBER *b, unsigned int length) {
	//	multiplication++;
	NUMBER sum = a[0] * b[0];
	for (int i = 1; i < length; i++) {
		//	addition++;
		//	multiplication++;
		sum += a[i] * b[i];
	}
	return sum;
}

void gs(NUMBER *w, NUMBER *scalar_v, unsigned int length, unsigned int lda) {
	NUMBER *v = w;
	int i, j, k;
	NUMBER scalar;
	NUMBER temp;
#ifdef ORTHOGONAL
	NUMBER ortho = 0;
#endif
	for (i = 0; i < length; i++) {
		for (k = i; k > 0; k--) {
			//	division++;
			scalar = dot(&W(i,0), &V(k-1,0), length) / scalar_v[k - 1];
			for (j = 0; j < length; j++) {
				//			addition++;
				//			multiplication++;
				temp = W(i,j) - scalar * V(k,j);
				V(i,j) = temp;
#ifdef ORTHOGONAL
				ortho += temp * temp;
#endif
			}
		}
#ifdef ORTHOGONAL
		ortho = sqrt(ortho);
		for (j = 0; j < length; j++) {
			V(i,j) = V(i,j) / ortho;
		}
#endif
		scalar_v[i] = dot(&V(i,0), &V(i,0), length);
	}
}

void initOrthoMatrix(NUMBER *v, unsigned int length, unsigned int lda) {
#ifdef RANDOMMATRIX
	srand(time(0));
	double random;
#endif
	for (int i = 0; i < length; i++) {
#ifdef RANDOMMATRIX
		random = (rand() + 1.0) / (double) RAND_MAX;
		random *= 2.0;
#endif
#ifdef OPENMP
#pragma omp parallel for
#endif
		for (int j = 0; j < length; j++) {
			if (i == j) {
				V(i,j) = 1;
#ifdef RANDOMMATRIX
			} else if (i < j) {
				V(i,j) = (NUMBER) (random * V(i,j-1));
#endif
			} else {
				V(i,j) = 0;
			}
		}
	}
}

void printMatrix(NUMBER *v, unsigned int length, unsigned int lda) {
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			printf("%f ", V(j,i));
		}
		printf("\n");
	}
}

long calcCPUOperations(int n, int add, int mult, int div) {
	double scalarproduct = n * mult + (n - 1) * add;
	double littlegaus = (n * (n - 1)) / 2;
#ifdef ORTHOGONAL
	double sqrt = div + n * div;
	double result = n * (scalarproduct + sqrt) + (littlegaus * (scalarproduct
			+ div + n * 2 * (add + mult)));
#else
	double result = n * (scalarproduct) + (littlegaus * (scalarproduct + div
					+ n * (add + mult)));
#endif
	return (long) result;
}

long calcAdd(int n) {
	return calcCPUOperations(n, 1, 0, 0);
}

long calcMult(int n) {
	return calcCPUOperations(n, 0, 1, 0);
}

long calcDiv(int n) {
	return calcCPUOperations(n, 0, 0, 1);
}

void printOperations(double time) {
	long calcadd = calcAdd(MATRIXDIMENSION);
	long calcmult = calcMult(MATRIXDIMENSION);
	long calcdiv = calcDiv(MATRIXDIMENSION);
	long allops = calcCPUOperations(MATRIXDIMENSION, 1, 1, 1);
	printf(
			"calculated:\n n: %d\n   Addtionen: %ld\n   Multiplikationen: %ld\n   Divisionen: %ld\n   all Operations: %ld\n\n",
			MATRIXDIMENSION, calcadd, calcmult, calcdiv, allops);
	double flops = allops / time;
	printf("flops: %f\n", flops);
	printf("mflops: %f\n", flops / 1000000);
	printf("gflops: %f\n", flops / 1000000000);
	//	printf("counted:\n n: %d\n   Addtionen: %d\n   Multiplikationen: %d\n   Divisionen: %d\n   all Operations: %ld\n", MATRIXDIMENSION, addition, multiplication, division,addition+multiplication+3*division);
}

int main(int argc, char **argv) {
#ifdef _OPENMP
	printf("OPENMP is activated\n");
#endif
	double time;
	int lda = MATRIXDIMENSION;
	NUMBER *w = (NUMBER*) malloc(sizeof(NUMBER) * MATRIXDIMENSION * lda);
	NUMBER *v = w;
	NUMBER *scalar_v = (NUMBER*) malloc(sizeof(NUMBER) * MATRIXDIMENSION);
	printf("start init:\n");
	printf("sizeof double: %d\n", (int)sizeof(double));
	printf("sizeof float: %d\n", (int)sizeof(float));
	printf("using precision: %d\n", (int)sizeof(NUMBER));
	initOrthoMatrix(v, MATRIXDIMENSION, lda);
	printf("start gs\n");

	/*	V(0,0) = 0;
	 V(0,1) = 0;
	 V(0,2) = 1;
	 V(1,0) = 0;
	 V(1,1) = 0;
	 V(1,2) = 1;*/
	/*printMatrix(v, MATRIXDIMENSION, lda);
	 printf("\n");*/

	//gs(w, scalar_v, MATRIXDIMENSION, lda);
	time = getTimeInSeconds();
	//for(int i=0;i<100;i++)
	gs(w, scalar_v, MATRIXDIMENSION, lda);
	time = getTimeInSeconds() - time;
	//printMatrix(v, MATRIXDIMENSION, lda);
	//printf("%2.0f\n", dot(&V(0,0), &V(1,0), MATRIXDIMENSION));
	printOperations(time);
	free(w);//free(w)==free(v)
	free(scalar_v);
	return 0;
}
