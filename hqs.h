// hqs.h - A header only library for hyper quicksort (a.k.a. parallel quicksort)
// By Oliver Kraft - oliver@kraften.net
#include <stdlib.h>
#include <string.h>
#include <omp.h>

# define dtype int

// Sequential pivot selection strategy. 
// 1 = Final element - sufficent in most cases 
// 2 = Random element
// 3 = Random 3-point median
#define sequential_pivot 1

// Parallel pivot selection strategy. 
// 1 = Median of thread 0
// 2 = Mean of medians
// 3 = Median of sorted medians
#define global_pivot 1


// Print array
void debug_print(int* arr, size_t sze, size_t id, size_t pool_sze) { 
    // For testing on small lists
    printf("%ld/%ld - [", id, pool_sze);

    for (size_t i = 0; i < sze - 1; i++) {
        printf("%d, ", arr[i]);
    }

    printf("%d]\n", arr[sze - 1]);
}


// Pregenerate random numbers
unsigned int* gen_random(size_t sze) {
    unsigned int* rns = (unsigned int*) malloc(sze * sizeof(unsigned int));

    for (size_t i = 0; i < sze; i++) {
        rns[i] = rand();
    }

    return rns;
}


// Quicksort
size_t partition(dtype* arr, size_t sze, size_t pi) {
    dtype tmp;

    dtype pivot = arr[pi];

    arr[pi] = arr[sze - 1];
    arr[sze - 1] = pivot;

    // Perform partition
    size_t i = 0;
    for (size_t j = 0; j < sze - 1; j++) {
        if (arr[j] < pivot) {
            tmp = arr[i];
            arr[i] = arr[j]; // Move j to unused leftmost position
            arr[j] = tmp; // Put i where j was
            i++;
        }
    }

    arr[sze - 1] = arr[i];
    arr[i] = pivot;

    return i;
}

void qs(dtype* arr, size_t sze, unsigned int* rns, size_t* r, size_t rns_sze) {
    if (sze >= 2) { // Base case: do nothing

        #if sequential_pivot == 1
        size_t pi = sze - 1; // Selects a last element as pivot index
        
        // METHOD 2: Random element
        #elif sequential_pivot == 2
        size_t pi = rns[*r] % sze;
        *r += 1;
        if (*r + 1 > rns_sze) {
            //printf("WARNING OUT OF RANDOM NUMBERS!");
            *r = 0;
        }
        
        // METHOD 3: Three point median
        #elif sequential_pivot == 3
        int A = rns[*r] % sze;
        int B = rns[*r + 1] % sze;
        int C = rns[*r + 2] % sze;
        *r += 3;
        if (*r + 3 > rns_sze) {
            //printf("WARNING OUT OF RANDOM NUMBERS!");
            *r = 0;
        }
        size_t pi;
        if ((arr[A] < arr[B]) ^ (arr[A] < arr[C])) { // XOR rules out A
            pi = A;
        } else if ((arr[B] < arr[A]) ^ (arr[B] < arr[C])) { // XOR rules out B
            pi = B;
        } else {
            pi = C;
        }
        #endif

        pi = partition(arr, sze, pi); // Partition and update pivot index

        qs(arr, pi, rns, r, rns_sze); // Quicksort < pivot
        qs(&arr[pi + 1], sze - pi - 1, rns, r, rns_sze); // Quicksort > pivot
    }
}

// Hyper quicksort
typedef struct thread_data_t {
    size_t id;
    size_t pool_sze; // Number of threads in current pool

    dtype** data; // First element in large array
    size_t* szes;

    size_t* pis; // Pivot indices for threads in current pool
    
    dtype* pivot; // Global pivot value
    dtype* pivots; // Pivots for current pool

    struct thread_data_t** all_threads;

    #if glob_pivot != 1
    size_t* r; // current random number
    unsigned int* rns;
    size_t rns_sze;
    #endif
} thread_data_t;

void _hqs(thread_data_t* t_data) {
    // Unpack struct for readability
    size_t id = t_data->id;
    size_t pool_sze = t_data->pool_sze;

    int** data = t_data->data;
    size_t* szes = t_data->szes;
    
    size_t* pis = t_data->pis;
    
    dtype* pivot = t_data->pivot;
    dtype* pivots = t_data->pivots;
    
    thread_data_t** all_threads = t_data->all_threads;

    #if global_pivot == 3
    size_t* r = t_data->r;
    unsigned int* rns  = t_data->rns;
    size_t rns_sze = t_data->rns_sze;
    #endif

    // Data for current thread
    int* arr = data[id];
    size_t sze = szes[id];
    
    /*
    // FOR DEBUG
    if (id == 0) {
        printf("Debug:\n");
    }
    #pragma omp barrier
    #pragma omp critical
    {
    debug_print(arr, (int) sze, id + 1, pool_sze);
    }*/

    if (pool_sze > 1) { // Base case only one thread in pool => sorted
        
        #if global_pivot != 1
        // Get pivot elements
        if (sze > 1) {
            pivots[id] = arr[sze/2]; // Median 
        } else {
            pivots[id] = -1; // Thread empty exclude
        }

        #pragma omp barrier // Wait for all threads before selecting pivot
        #endif

        if (id == 0) { // Thread 0 selects global pivot

            #if global_pivot == 1 // Median of thread 0
            *pivot = arr[sze/2]; 
            
            #elif global_pivot == 2 // Mean of medians
            dtype sum = 0; // Add overflow check to loop here
            size_t excluded = 0;
            for (size_t i = 0; i < pool_sze; i++) {
                if (pivots[i] != -1) {
                    sum += pivots[i];
                } else {
                    excluded++;
                }
            }

            *pivot = sum/(pool_sze - excluded);

            // Approach 3:
            #elif global_pivot == 3 // Median of medians

            dtype tmp_pivots[pool_sze];
            size_t tmp_sze = 0;
            for (size_t i = 0; i < pool_sze; i++){
                if (pivots[i] != -1) {
                    tmp_pivots[tmp_sze] = pivots[i];
                    tmp_sze++;
                }
            }
            
            qs(tmp_pivots, tmp_sze, rns, r, rns_sze);

            if (tmp_sze % 2 == 0) {
                *pivot = (pivots[tmp_sze/2]/2 + pivots[tmp_sze/2 - 1]/2);
            } else {
                *pivot = pivots[pool_sze/2];
            }
            #endif

        }
        
        #pragma omp barrier // Wait for global pivot selection

        // Get indices of pivot elements
        size_t i = 0;
        for (i = 0; i < sze; i++) { // Replace with binary search
            if (arr[i] > *pivot) {
                pis[id] = i;
                break;
            }
        }
        if (i == sze) pis[id] = i;

        #pragma omp barrier // Wait for all threads to obtain indices

        int half_pool = pool_sze/2;

        dtype* k1; // Elements smaller than pivot
        dtype* k2;

        int s1; // No. of smaller elements
        int s2;

        if (id < half_pool) { 
            int p1 = pis[id];
            int p2 = pis[id + half_pool];

            s1 = p1;
            s2 = p2;

            k1 = &data[id][0];
            k2 = &data[id + half_pool][0]; 

        }

        if (id >= half_pool) {
            int p1 = pis[id];
            int p2 = pis[id - half_pool];

            s1 = szes[id] - p1;
            s2 = szes[id - half_pool] - p2;

            k1 = &data[id][p1];
            k2 = &data[id - half_pool][p2];
        }

        dtype* merged = (int*) malloc((s1 + s2) * sizeof(int)); // TODO: Replace with realloc
        size_t merged_sze;

        {
            size_t i, j1, j2;
            i = 0;
            j1 = 0;
            j2 = 0;
            while (j1 < s1 && j2 < s2) {
                //printf("%d/%d\n", j1, j2);
                if (k1[j1] < k2[j2]) {
                    merged[i] = k1[j1];
                    j1++;
                } else {
                    merged[i] = k2[j2];
                    j2++;
                }
                i++;
            }
            if (j1 == s1 && j2 < s2) {
                memcpy(&merged[i], &k2[j2], (s2 - j2) * sizeof(dtype));
                i += s2 - j2;
            } else if (j1 < s1) {
                memcpy(&merged[i], &k1[j1], (s1 - j1) * sizeof(dtype));
                i += s1 - j1;
            }
            merged_sze = i;
        }

        #pragma omp barrier // Wait for distribution before freeing old data arrays

        free(data[id]);
        data[id] = merged;
        szes[id] = merged_sze;

        // Create new pools!

        if (id==half_pool) { // Update threads that will deal with larger elements
            // Allocate new arrays

            dtype** new_data = &data[half_pool];   //(int**) malloc(half * sizeof(int*));
            size_t* new_szes = &szes[half_pool];    //(int*) malloc(half * sizeof(int*));

            size_t* new_pis = &pis[half_pool];

            dtype* new_pivot = (int*) malloc(sizeof(int));
            dtype* new_pivots = &pivots[half_pool];  //(int*) malloc(half * sizeof(int));

            thread_data_t** new_all_threads = &all_threads[half_pool];  //(thread_data_t**) malloc(half * sizeof(thread_data_t*));

            // Update threads that will deal with larger elements
            int offset = half_pool;
            for (int j = half_pool - 1; j >= 0; j--) {  
                new_data[j] = data[j + offset];
                new_szes[j] = szes[j + offset];
                
                thread_data_t* thread_i = t_data->all_threads[j + offset];

                //printf("#%d : id = %d, new_sze = %d\nNew_szes: ", id, thread_i->id, new_szes[j]);
                
                thread_i->id = j;
                thread_i->pool_sze = half_pool;

                thread_i->data = new_data;
                thread_i->szes = new_szes;

                thread_i->pis = new_pis;

                thread_i->pivot = new_pivot;
                thread_i->pivots = new_pivots;
                
                new_all_threads[j] = thread_i; // Add thread to new threadpool
                thread_i->all_threads = new_all_threads; // Replace pointer to current thread pool with pointer to new threadpool
            }
        }
        if (id == 0) { // Update threads that will deal with smaller elements
            
            for (int j = half_pool - 1; j >= 0; j--) {
                thread_data_t* thread_i = t_data->all_threads[j];

                thread_i->id = j;
                thread_i->pool_sze = half_pool;
            }
        }

        #pragma omp barrier // Wait for new pools to be created before proceeding
        
        _hqs(t_data); // Recursion call!

    }

}

void hqs(dtype* arr, size_t total_sze, size_t nthreads) {

    size_t* szes = (size_t*) malloc(nthreads * sizeof(size_t));
    dtype** data = (dtype**) malloc(nthreads * sizeof(dtype*));

    size_t* pis = (size_t*) malloc(nthreads * sizeof(size_t));

    dtype* pivot = (dtype*) malloc(sizeof(dtype));
    dtype* pivots = (dtype*) malloc(nthreads * sizeof(dtype));

    thread_data_t* all_threads = (thread_data_t*) malloc(nthreads * sizeof(thread_data_t));
    thread_data_t** ptrs = (thread_data_t**) malloc(nthreads * sizeof(thread_data_t*));

    #pragma omp parallel num_threads(nthreads)
    {

        int id = omp_get_thread_num();

        if (id == 0) {
            printf("Number of threads: %d\n", omp_get_num_threads());
        }

        #if glob_pivot != 1 || seq_pivot != 1
        // Get some random numbers for current thread
        size_t rns_sze = 1024; // Could probably be higher
        unsigned int* rns = (unsigned int*) malloc(rns_sze * sizeof(unsigned int));
        size_t r = 0; // Current random number
        #pragma omp critical
        {
            for (int i = 0; i < rns_sze; i++) rns[i] = rand();
        }
        #else
        // No random numbers needed, placeholder vars
        int rng_max = 0; // Could probably be higher
        int* rng = 0;
        int r = 0;
        #endif

        #if glob_pivot != 1 || seq_pivot != 1
        all_threads[id].rns = rns;
        all_threads[id].r = &r;
        all_threads[id].rns_sze = rns_sze;
        #endif

        // Get data for current thread
        dtype* current_data = (dtype*) malloc((total_sze + nthreads)/nthreads * sizeof(dtype));
        size_t sze = 0;
        
        #pragma omp for schedule(static)
        for (int i = 0; i < total_sze; i++) {
            //printf("%d : %d\n", id,i);
            current_data[sze] = arr[i];
            sze++;
        }
        
        data[id] = current_data;
        all_threads[id].id = id;
        all_threads[id].pool_sze = nthreads;
        all_threads[id].data = data;
        all_threads[id].szes = szes;
        szes[id] = sze;

        all_threads[id].pis = pis;

        all_threads[id].pivot = pivot;
        all_threads[id].pivots = pivots;
        
        ptrs[id] = &all_threads[id];
        all_threads[id].all_threads = ptrs;


        qs(current_data, sze, rns, &r, rns_sze); 

        #pragma omp barrier // Wait for all threads to perform local quicksort

        _hqs(&all_threads[id]);

        #pragma omp barrier
        #pragma omp single
        {
            int index = 0;

            for (int i = 0; i < nthreads; i++) {
                int s = all_threads[i].szes[0];
                int* d = all_threads[i].data[0];
                memcpy(&arr[index], d, s * sizeof(int));
                index += all_threads[i].szes[0];

                // Free memory                
                free(d); // Free actual data array
                free(all_threads[i].pivot);
                
            }
            free(data); // Free array of data arrays
            free(szes); // Free size array 
            free(pis); // Free pivot index array
            free(pivots); // Free pivots
            free(ptrs); // Free array of thread pointers
            free(all_threads); // Free thread array
            
        }
        #if glob_pivot != 1 || seq_pivot != 1
        free(rns);
        #endif
    }
}

