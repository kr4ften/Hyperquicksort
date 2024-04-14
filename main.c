#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "hqs.h"
#include <math.h>

// Return current time in seconds, for use in benchmarking
static double curr_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds  = tv.tv_sec + (double) tv.tv_usec / 1000000;
    return seconds;
}

// Generate an array of random ints
int* gen_arr(int sze) {
    
    int* arr = (int*) malloc(sze* sizeof(int));

    for (int i = 0; i < sze; i++) 
        arr[i] = rand() % 1000000;
    
    return arr;
}

// Print array
void print_arr(int* arr, int sze) { 
    // For testing on small lists
    printf("[");

    for (int i = 0; i < sze - 1; i++) {
        printf("%d, ", arr[i]);
    }

    printf("%d]\n", arr[sze - 1]);
}

// Test if array is sorted
void test_sorted(int* arr, int sze) {
    int test = 1;
    for (int i = 0; i < sze - 1; i++) {
        if (arr[i] > arr[i + 1]) {
            printf("%d => %d\n", arr[i], arr[i + 1]);
            test = 0;
        }
    }

    if (test == 1) {
        printf("Sorted ok!\n");
    } else {
        printf("Sorted not ok!\n");
    }
}

int main (int argc, char* argv[]) {
    size_t sze;
    size_t thread_count;
    if (argc != 3) {
        printf("Incorrect number of arguments! Usage: ./quicksort N thread_count\n");   
        return -1;
    } else {
        sze = atoi(argv[1]);
        thread_count = atoi(argv[2]);
        printf("Running hyper quicksort with %ld elements, and %ld threads...\n", sze, thread_count);
    }

    printf("thread_count = %ld \n", thread_count);

    // Check input is power of two (and not zero)
    if ((thread_count & (thread_count - 1)) != 0) {
        printf("Incorrect usage, thread_count must be a power of 2.\n");
        return -1;
    }

    // Set random seed
    srand(42);

    int* arr = gen_arr(sze);

    double start_time = curr_time();

    hqs(arr, sze, thread_count);

    printf("time: %f\n", curr_time() - start_time);

    test_sorted(arr, sze);

    //print_arr(arr, sze);

    free(arr);

    return 0;

}

