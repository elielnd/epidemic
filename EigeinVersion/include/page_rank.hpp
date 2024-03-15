#include <fstream>
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <random>
#include "../eigen/Eigen/Sparse"



#define ALPHA 0.85
#define NU 0.2
#define DELTA 0.24
#define INITIAL_INFECTED_PERCENTAGE 0.2
#define VACCINATION_PERCENTAGE 0.35
#define TOL 1e-6
#define RUN 5
#define OUTPUT_DIR "output/"
#define MAX_ITER 1000
#define WO_VACCINE_FILE "wo_vaccine.csv"
#define W_VACCINE_FILE "w_vaccine.csv"
#define USING_INFECTED_VECTOR_FILE "infected.csv"
#define OUTPUT_FILE "eigenvector.txt"
#define INF_RANK "infected_rank.txt"
#define TIME 100
#define BELLKIS 0.01

Eigen::SparseMatrix<double> readAdjacencyMatrix(const std::string &filename);
void writeMatrixToFile(const std::string &filename, const Eigen::VectorXd &vector);
std::pair<Eigen::VectorXd, double> powerIteration(const Eigen::SparseMatrix<double> &transitionMatrix, double alpha, double tol, int max_iter, int &iterations);
Eigen::SparseMatrix<double> buildTransitionMatrix(const Eigen::SparseMatrix<double> &inputMatrix);