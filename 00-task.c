#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
  int i, j;
  int N = atoi(argv[1]);
  int* arr = malloc(N * sizeof(int));
  int current_sum = 0, ans_sum = 0, previous = 0;
  int rank, size;
  MPI_Status Status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int distance  = N / (size - 1), remainder  = N % (size - 1);
  
  if(rank != 0){
        int to_recv = distance;
        if (rank <= remainder) {
            to_recv = distance + 1;
        }
        MPI_Recv(&arr[0], to_recv, MPI_INT, 0, rank, MPI_COMM_WORLD, &Status);
        for (j = 0; j < to_recv; j++) {
            current_sum += arr[j];
        }
        printf("For rank=%d and size=%d, answer sum = %d\n", rank, size, current_sum);
        MPI_Send(&current_sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
  }
  else{
    for(i = 0; i < N; i++) {
      arr[i] = i;
      previous += i;
    }
    int to_send = distance + 1;
    int start = 0;
    for (i = 1; i < size; i++) {
      MPI_Send(&arr[start], to_send, MPI_INT, i, i, MPI_COMM_WORLD);
      if (i <= remainder) {
        start += (distance + 1);
      }
      else {
        start += distance;
      }
      if (i >= remainder) {
        to_send = distance;
      }
    }
    for (i = 1; i < size; i++) {
      MPI_Recv(&current_sum, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &Status);
      ans_sum += current_sum;
    }
    printf("For rank = %d and size = %d, , current sum = %d; answer sum = %d\n", rank, size, ans_sum, previous);
  }
  
  free(arr);
  MPI_Finalize();
  return 0;
}
