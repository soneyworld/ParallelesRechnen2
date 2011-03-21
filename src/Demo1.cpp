#include "config.h"
#ifdef OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <unistd.h>
#include "Benchmark.h"

void test1() {
	int kerne, nummer;
	kerne = omp_get_num_procs();
	omp_set_num_threads(kerne);
	omp_set_dynamic(0);
#pragma omp parallel private(nummer)
	{
		nummer = omp_get_thread_num();
		printf("Thread %d: works\n", nummer);
#pragma omp barrier
#pragma omp single
		{
			printf("Thread %d: sums\n", nummer);
		}
		printf("Thread %d: works on\n", nummer);
	}
}

void test2() {
	int kerne, nummer;
	kerne = omp_get_num_procs();
	omp_set_num_threads(kerne);
	omp_set_dynamic(0);
#pragma omp parallel private(nummer) num_threads(kerne)
	{
		nummer = omp_get_thread_num();
		printf("Thread: %d\n", nummer);
	}
}

void test3() {
	int sum;
	int kerne, nummer;
	kerne = omp_get_num_procs();
	omp_set_num_threads(kerne);
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

}
void test4() {
	int kerne, tid;
	kerne = omp_get_num_procs();
	omp_set_num_threads(kerne);
#pragma omp parallel private(tid) num_threads(kerne)
#pragma omp sections
	{
		tid = omp_get_thread_num();
#pragma omp section
		{
			sleep(1);
			printf("Thread %d\n", tid);
		}
#pragma omp section
		{
			printf("Thread %d\n", tid);
		}
	}
}
int fib(int n) {
	int i, j;
	if (n < 2)
		return (n);
	else {
		i = fib(n - 1);
		j = fib(n - 2);
		return i + j;
	}
}
int parafib(int n) {
	int i, j;
	if (n < 2)
		return (n);
	else {
#pragma omp task untied shared(i)
		i = fib(n - 1);
#pragma omp task untied shared(j)
		j = fib(n - 2);
#pragma omp taskwait
		return i + j;
	}
}

void test5(int x) {
	int kerne, tid;
	kerne = omp_get_num_procs();
	omp_set_num_threads(kerne);
	omp_set_dynamic(1);
	int result;
#pragma omp parallel num_threads(kerne)
	{
#pragma omp single
		result = parafib(x);
	}
	printf("fib(%i)=%i \n", x, result);
}
void test6(int x) {
	int result;
	result = fib(x);
	printf("fib(%i)=%i \n", x, result);
}

void test7() {
	int x = 20;
	double time;
	time = getTimeInSeconds();
	test5(x);
	time = getTimeInSeconds() - time;
	printf("%f seconds parallel\n", time);
	time = getTimeInSeconds();
	test6(x);
	time = getTimeInSeconds() - time;
	printf("%f seconds sequenziell\n", time);
}

void test8() {
	double time;
	time = getTimeInSeconds();
#pragma omp parallel num_threads(2)
	{
#pragma omp single
		{
#pragma omp task
			{
				sleep(1);
			}
#pragma omp task
			{
				sleep(1);
			}
#pragma omp taskwait
			{

			}
		}
	}
	time = getTimeInSeconds() - time;
	printf("TEST 8: %f seconds parallel\n", time);
}

int main(int argc, char **argv) {
#ifdef OPENMP
	test1();
	test2();
	test3();
	test4();
	test7();
	test8();
#endif
	return 0;
}
