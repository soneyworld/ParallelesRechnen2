#include "config.h"
#ifdef OPENMP
#include <omp.h>
#endif
#include <cstdio>

int main(int argc, char **argv) {
#ifdef OPENMP
	int kerne,nummer;
	kerne = omp_get_num_procs();
	omp_set_num_threads(kerne);
#pragma omp parallel private(nummer) num_threads(kerne)
	{
		nummer = omp_get_thread_num();
		printf("Thread: %d\n",nummer);
	}

	int i, sum;
#pragma omp parallel for private(i) reduction(+:sum) num_threads(kerne)
	for(i=0;i<100;i++){
		sum +=1;
	}
	printf("Summe: %d\n",sum);
#endif
	return 0;
}
