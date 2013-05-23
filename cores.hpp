#ifndef volmodels_cores_hpp
#define volmodels_cores_hpp

#ifdef _OPENMP
    #include <omp.h>
#endif

int numberOfCoresUsed()
{
    int cores = 1;
    
    #ifdef _OPENMP

    #pragma omp parallel 
    {
       #pragma omp master
       cores = omp_get_num_threads();
	}
    
    #endif
    
    return cores;
}

#endif
