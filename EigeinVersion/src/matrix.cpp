#include "../include/page_rank.hpp"

using namespace Eigen;

// Fonction pour lire une matrice d'adjacence à partir d'un fichier
Eigen::SparseMatrix<double> readAdjacencyMatrix(const std::string &filename)
{
    std::ifstream file(filename);
    int numNodes = 0;
    std::vector<Triplet<double>> triplets;

    if (file.is_open())
    {
        std::string line;
        while (std::getline(file, line))
        {
            int fromNode = 0, toNode = 0;
            std::istringstream iss(line);
            try
            {
                if (iss >> fromNode >> toNode)
                {
                    // Code to execute if the line inputFile >> fromNode >> toNode) is successfulll
                    triplets.push_back(Triplet<double>(fromNode, toNode, 1));
                    numNodes = std::max({numNodes, fromNode, toNode});
                }
            }
            catch (const std::exception &e)
            {
                continue;
            }
        }
        file.close();

        if (numNodes == 0)
        {
            std::cerr << "No nodes found in the input file" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // Création d'une matrice creuse pour stocker la matrice d'adjacence
    Eigen::SparseMatrix<double> adjacencyMatrix(numNodes + 1, numNodes + 1);
    adjacencyMatrix.setFromTriplets(triplets.begin(), triplets.end());

    return adjacencyMatrix;
}

void writeMatrixToFile(const std::string &filename, const VectorXd &vector)
{
    std::ofstream file(filename);
    if (file.is_open())
    {
        for (int i = 0; i < vector.size(); i++)
        {
            file << vector(i) << std::endl;
        }
        file.close();
    }
    else
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Fonction pour construire une matrice de transition à partir d'une matrice creuse
Eigen::SparseMatrix<double> buildTransitionMatrix(const Eigen::SparseMatrix<double> &inputMatrix)
{
    Eigen::SparseMatrix<double> inputMatrixTranpose = inputMatrix.transpose();
    int numRows = inputMatrixTranpose.rows();
    int numCols = inputMatrixTranpose.cols();

    // Création d'une matrice creuse pour stocker la matrice de transition
    Eigen::SparseMatrix<double> transitionMatrix(numRows, numCols);

    // Parcours de chaque colonne de la matrice creuse d'entrée
    for (int col = 0; col < numCols; ++col)
    {
        // Extraction de la colonne courante
        Eigen::SparseVector<double> currentColumn = inputMatrixTranpose.col(col);

        // Calcul de la somme des éléments de la colonne courante
        double sum = currentColumn.sum();

        // Si la somme est non nulle, division de la colonne par la somme
        if (sum != 0.0)
        {
            currentColumn /= sum;
        }
        // Assignation de la colonne modifiée à la matrice de transition
        transitionMatrix.col(col) = currentColumn;
    }

    return transitionMatrix;
}

Eigen::SparseMatrix<double> buildMatrixG(const Eigen::SparseMatrix<double> &A)
{
    // Création d'une matrice creuse B avec des éléments non nuls égaux à 1/n
    Eigen::SparseMatrix<double> B(A.rows(), A.cols());
    B.reserve(VectorXi::Constant(A.cols(), A.rows()));

    int n = A.rows();
    double value = 1.0 / n;
    for (int k = 0; k < A.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
        {
            B.coeffRef(it.row(), it.col()) = value;
        }
    }

    B.makeCompressed(); // Finalisation de la construction de la matrice creuse B

    return B;
}

std::pair<VectorXd, double> powerIterationbis(const SparseMatrix<double> &A, double tol, int max_iter, int &iterations)
{
    int n = A.rows();
    VectorXd x = VectorXd::Random(n);
    x = x.cwiseAbs(); // Prendre la valeur absolue de chaque élément du vecteur
    x /= x.norm();

    VectorXd y = A * x;
    double lambda = x.dot(y);

    VectorXd r = y - lambda * x;
    int iter = 0;

    while (r.norm() > tol && iter < max_iter)
    {
        x = y / y.norm();
        y = A * x;
        lambda = x.dot(y);
        r = y - lambda * x;
        iter++;
    }

    iterations = iter;
    return std::make_pair(x, lambda);
}

std::pair<VectorXd, double> powerIteration(const SparseMatrix<double> &transitionMatrix, double alpha, double tol, int max_iter, int &iterations)
{
    int n = transitionMatrix.rows();
    double beta = (1 - alpha) / n;
    double lambda = 0;
    iterations = 0;

    // Vecteur initial x
    VectorXd x = VectorXd::Random(n);
    x = x.cwiseAbs(); // Prendre la valeur absolue de chaque élément du vecteur
    x /= x.norm();

    // Creez un vecteur v rempli de 1
    VectorXd v = VectorXd::Ones(n);

    VectorXd y = VectorXd::Zero(n);

    VectorXd r = x - y;

    while (r.norm() > tol && iterations < max_iter)
    {
        double sum = x.sum();
        y = alpha * transitionMatrix * x + beta * sum * v;
        lambda = x.dot(y);
        r = y - x;
        x = y / y.norm();
        iterations++;
    }

    return std::make_pair(x, lambda);
}