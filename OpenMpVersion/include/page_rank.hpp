#include <iostream>
#include <omp.h>
#include "matrix.hpp"
#
#include <sstream>
#include <fstream>
#include <unistd.h>


#define ALPHA 1
#define TOL 1e-6
#define MAX_ITER 2000