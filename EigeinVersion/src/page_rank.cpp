#include "../include/page_rank.hpp"

using namespace Eigen;


int main(int argc, char *argv[])
{
    double alpha = ALPHA;
    double tol = TOL;
    std::string eigenvector_filename = OUTPUT_FILE;

    // Vérifier le nombre d'arguments
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <adjacency_matrix_file>" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (argc >= 3)
    {
        alpha = std::atof(argv[2]);
    }

    if (argc >= 4)
    {
        tol = std::atof(argv[3]);
    }

    if (argc >= 5)
    {
        eigenvector_filename = std::atoi(argv[4]);
    }

    //***********************************************
    // Lire la matrice d'adjacence à partir du fichier
    //***********************************************
    std::string filename = argv[1];
    Eigen::SparseMatrix<double> adjacencyMatrix = readAdjacencyMatrix(filename);
    int nonZeroCount = adjacencyMatrix.nonZeros();
    
    //***********************************************
    // Construire la matrice de transition P
    //***********************************************
    SparseMatrix<double> transitionMatrixP = buildTransitionMatrix(adjacencyMatrix);

  

    //***********************************************
    // Appliquer la méthode de la puissance pour obtenir le vecteur de PageRank
    //***********************************************
    int iter = 0;
    auto start = std::chrono::high_resolution_clock::now();
    auto result = powerIteration(transitionMatrixP, alpha, tol, MAX_ITER, iter);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // Récupérer les résultats
    VectorXd x = result.first;
    double lamda = result.second;

    // Recuperer les valeurs minimales et maximales du vecteur de PageRank ainsi que leurs indices
    double max_element = x.maxCoeff();
    double max_index = -1;
    double min_element = x.minCoeff();
    double min_index = -1;

    for (int i = 0; i < x.size(); i++)
    {
        if (x(i) == max_element)
        {
            max_index = i;
        }
        if (x(i) == min_element)
        {
            min_index = i;
        }
    }

    // Affichage des résultats
    std::cout << "--- Résultats de l'Algorithme de PageRank --- \n"
              << "\nConfiguration : \n"
              << "\t-Taille de la matrice d'adjacence : " << adjacencyMatrix.rows() << "*" << adjacencyMatrix.cols() << "\n"
              << "\t-Nombre d'éléments non nuls : " << nonZeroCount << "\n"
              << "\nParamètres de l'algorithme : \n"
              << "\t-Tolerance : " << tol << "\n"
              << "\t-Damping factor(alpha) : " << alpha << "\n"
              << "\nResultats : \n"
              << "Nombre d'itérations effectuées : " << iter << "\n"
              << "Valeur propre dominante : " << lamda << "\n"
              << "Vecteur propre correspondant: \n"
              << "\t-Maximum value: " << max_element << " (Index : " << max_index << ")" << std::endl;

        //***********************************************
    // Sauvegarder les résultats dans un fichier
    //***********************************************

    x.normalize();
    eigenvector_filename = OUTPUT_DIR + eigenvector_filename;
    writeMatrixToFile(eigenvector_filename, x);
    std::cout << "Ecriture du vecteur propre dans le fichier : " << eigenvector_filename << std::endl;

    std::cout << "\nTemps d'exécution de power_iteration : " << duration << " ms" << std::endl;
    return 0;
}

