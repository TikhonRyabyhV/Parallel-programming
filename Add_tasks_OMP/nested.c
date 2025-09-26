#include <stdio.h>
#include <omp.h>

#define MAX_LEVELS 3
#define THREADS_PER_LEVEL 2

int main() {
    omp_set_nested(1);
    omp_set_max_active_levels(MAX_LEVELS);

    // 1st level
    #pragma omp parallel num_threads(THREADS_PER_LEVEL)
    {
        int level1_thread = omp_get_thread_num();
        int level1_num_threads = omp_get_num_threads();
        int level = omp_get_level();

        printf("Level %d, thread %d: summary %d threads, active level: %d\n",
               level, level1_thread, level1_num_threads, omp_get_active_level());

	// 2nd level
        #pragma omp parallel num_threads(THREADS_PER_LEVEL)
        {
            int level2_thread = omp_get_thread_num();
            int level2_num_threads = omp_get_num_threads();
            int level = omp_get_level();
            int parent_level1 = omp_get_ancestor_thread_num(1);

            printf("Level %d, thread %d (parent: level 1, thread %d): "
                   "summary %d threads, level 1: %d, active level: %d\n",
                   level, level2_thread, parent_level1, level2_num_threads,
                   omp_get_team_size(1), omp_get_active_level());

	    // 3d level
            #pragma omp parallel num_threads(THREADS_PER_LEVEL)
            {
                int level3_thread = omp_get_thread_num();
                int level3_num_threads = omp_get_num_threads();
                int level = omp_get_level();
                int parent_level1 = omp_get_ancestor_thread_num(1);
                int parent_level2 = omp_get_ancestor_thread_num(2);

                printf("Level %d, thread %d (parents: level 2, thread %d; level 1, thread %d): "
                       "summary %d threads, level 2: %d, level 1: %d, active level: %d\n",
                       level, level3_thread, parent_level2, parent_level1,
                       level3_num_threads, omp_get_team_size(2), omp_get_team_size(1),
                       omp_get_active_level());
            }
        }
    }

    return 0;
}
