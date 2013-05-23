// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

void checkOpenMP()
{
#ifdef _OPENMP
    
#pragma omp parallel 
    {
#pragma omp master
        std::cout << "OpenMP enabled. Threads: " << omp_get_num_threads() << std::endl;
	}
    
#else
    
    std::cout << "OpenMP disabled." << std::endl;
    
#endif            
}

int main(int argc, char *argv[]) {
    checkOpenMP();
}
