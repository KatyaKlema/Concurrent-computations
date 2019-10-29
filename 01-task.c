#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
double temp_prev = 0, temp_next = 0, t0 = 0, k = 1, h = 0.02, dt = 0.0002, pi = 3.14159;

void fillAnswer(double *answer, int N){
    int i;
    for(i = 0; i < N; ++i){
        int p;
        double sum = 0;
        for (p = 0; p < 1000; p++) {
            sum += (exp(-1 * k * pi * pi * (2 * p + 1) * (2 * p + 1) * 0.1)) * (sin(pi * (2 * p + 1) * (i) * h)) / (2 * p + 1);
        }
        answer[i] = 4 * sum / pi;
    }
}


void  updateArray(double *arr, int len, double *temp_prev, double *previous){
    int i = 1;
    for (i = 1; i < len - 1; i++) {
        *temp_prev = *previous;
        *previous = arr[i];
        arr[i] = arr[i] + ((k * dt) / (h * h)) * (arr[i + 1] - 2 * arr[i] + *temp_prev);
    }

}


int main(int argc, char** argv) {
  int N = atoi(argv[1]);
  double* arr = malloc(N * sizeof(double));
  double* solution = malloc(N * sizeof(double));
 
  int i;
  for(i = 0; i < N; ++i){
        arr[i] = 1;
  }
    
  fillAnswer(solution, N);
  
  int rank, size;
  MPI_Status Status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
  int difference = N / size, remainder = N % size;
  int len = difference;
  if(rank < remainder)
      len++;

  float t = 0;
  double previous = 0;
  for (t = 0; t <= 0.1; t += 0.0002) {
    if(rank == size - 1){
        MPI_Send(&arr[0], 1, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD);
        MPI_Recv(&temp_prev, 1, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD, &Status);
        previous = arr[0];
        arr[0] = arr[0] + ((k * dt) / (h * h)) * (arr[1] - 2 * arr[0] + temp_prev);
        updateArray(arr, len, &temp_prev, &previous);
        arr[len - 1] = arr[len - 1] + ((k * dt) / (h * h)) * (0 - 2 * arr[len - 1] + previous);
    }
    else if (rank == 0) {
      MPI_Send(&arr[len - 1], 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
      MPI_Recv(&temp_next, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);
      previous = arr[0];
      arr[0] = arr[0] + ((k * dt) / (h * h)) * (arr[1] - 2 * arr[0] + 0);
      updateArray(arr, len, &temp_prev, &previous);
      arr[len - 1] = arr[len - 1] + ((k * dt) / (h * h)) * (temp_next - 2 * arr[len - 1] + previous);
    }
    else {
      MPI_Send(&arr[len - 1], 1, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD);
      MPI_Send(&arr[0], 1, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD);
      MPI_Recv(&temp_prev, 1, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD, &Status);
      MPI_Recv(&temp_next, 1, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD, &Status);
      previous = arr[0];
      arr[0] = arr[0] + ((k * dt) / (h * h)) * (arr[1] - 2 * arr[0] + temp_prev);
      updateArray(arr, len, &temp_prev, &previous);
      arr[len - 1] = arr[len - 1] + ((k * dt) / (h * h)) * (temp_next - 2 * arr[len - 1] + previous);
    }
  }
  if(rank == 0){
      for (i = 1; i < remainder; i++) {
          MPI_Recv(&arr[(difference + 1) * i], difference + 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
      }
      if (remainder == 0)
          i = 1;
      else
          i = remainder;
      
      for (; i < size; i++) {
          MPI_Recv(&arr[difference * i + remainder], difference, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
      }
      for (i = 0; i < N; i += 5) {
          printf("%.6f  %.6f\n", arr[i], solution[i]);
      }
  }
  else {
    MPI_Send(&arr[0], len, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
    
  free(arr);
  MPI_Finalize();
  return 0;
}

