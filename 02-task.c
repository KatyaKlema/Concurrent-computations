#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

#define PRECISE (double)1e-8

double f(double arg) {
    return 1 / (1 + arg * arg);
}


double criticalsection(char * number) {
    int i = 0;
    double integral = 0;

    int N = atoi(number);
    omp_set_num_threads(N);
    N = omp_get_num_threads();
    int myid = omp_get_thread_num();
#pragma omp parallel shared(integral), private(i, myid)
    {
#pragma omp for
        for (i = 0; i < (int)(1 / PRECISE); i++) {
            double el = (f((i + 1) * PRECISE) + f(i * PRECISE));
#pragma omp critical
            {
                integral += el;
            }

        }
    }
}

double reduction(char * number) {
    int i = 0;
    double integral = 0;

    int num_threads = atoi(number);
    omp_set_num_threads(num_threads);
    num_threads = omp_get_num_threads();

#pragma omp parallel shared(integral), private(i)
    {
#pragma omp for reduction(+: integral)
        for (i = 0; i < (int)(1 / PRECISE); i++) {
            double el = (f((i + 1) * PRECISE) + f(i * PRECISE));
            integral += el;
        }
    }
    return 2 * integral * PRECISE;
}

int main(int argc, char* argv[]) {
#ifdef _OPENMP
    printf("OpenMP is supported! %d \n", _OPENMP);
#endif
    int i = 0;
    if (argc < 2) {
        printf("To few arguments : %d, expected : 2 \n", argc - 1);
        return -1;
    }
    double integral = criticalsection(argv[1]);
    printf("integral = %f\n", 2 * integral * PRECISE);
}
