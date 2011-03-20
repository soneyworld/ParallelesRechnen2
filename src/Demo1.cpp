#include "config.h"
#ifdef OPENMP
#include <omp.h>
#endif
#include <cstdio>

int main(int argc, char **argv) {
#ifdef OPENMP
	int kerne, nummer;
	kerne = omp_get_num_procs();
	omp_set_num_threads(kerne);
#pragma omp parallel private(nummer) num_threads(kerne)
	{
		nummer = omp_get_thread_num();
		printf("Thread: %d\n", nummer);
	}
	int sum;
#pragma omp parallel private(nummer) num_threads(kerne)
	{
		int i;
		nummer = omp_get_thread_num();
#pragma omp for reduction(+:sum) schedule(static,1)
		for (i = 0; i < 100; i++) {
			sum += 1;
			printf("Thread %d: %d\n", nummer, sum);
		}
	}
	printf("Summe: %d\n", sum);
#endif
	return 0;
}
