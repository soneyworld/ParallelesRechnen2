#include "config.h"
#ifdef OPENMP
#include "Dot.h"
#include "Benchmark.h"
#include <cstdio>
#include <omp.h>

#define MATRIXDIMENSION 1000000
#define W(i,j) w[(i)*lda+(j)]

NUMBER dot(NUMBER *a, NUMBER *b, unsigned int length) {
	NUMBER sum;
	int i;
	for (i = 0; i < length; i++)
		sum += a[i] * b[i];
	return sum;
}

NUMBER paradot(NUMBER *a, NUMBER *b, unsigned int length) {
	NUMBER sum;
	int i,kerne = omp_get_num_threads();
#pragma omp parallel for private(i) reduction(+:sum) num_threads(kerne)
	for (i = 0; i < length; i++)
		sum += a[i] * b[i];
	return sum;
}
#endif

int main(int argc, char **argv) {
#ifdef OPENMP
	int lda = 2;
	NUMBER *w = (NUMBER*) malloc(sizeof(NUMBER) * MATRIXDIMENSION * lda);
	double time;
	double openmptime;
	time = getTimeInSeconds();
	paradot(w,w+MATRIXDIMENSION,MATRIXDIMENSION);
	openmptime = getTimeInSeconds() - time;
	time = getTimeInSeconds();
	dot(w,w+MATRIXDIMENSION,MATRIXDIMENSION);
	time = getTimeInSeconds()-time;
	printf("openmp: %f\n"
			"iterativ: %f\n",openmptime, time);
	printf("Speedup: %f\n",(time/openmptime));
#endif
	return 0;
}
