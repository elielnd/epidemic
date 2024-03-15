#include "../include/page_rank.hpp"

// Fonction pour lire le fichier et remplir la matrice
Matrix readMatrixFromFile(const std::string &filename)
{
    std::ifstream inputFile(filename);
    if (!inputFile.is_open())
    {
        std::cerr << "Impossible d'ouvrir le fichier. \n";
        exit(EXIT_FAILURE);
    }

    int numNodes = 0, fromNode = 0, toNode = 0;
    std::string line;
    while (getline(inputFile, line))
    {
        if (line.empty() || line[0] == '#')
        {
            // Ignorer les lignes vides ou celles commençant par #
            continue;
        }

        // Extraire les nœuds des arêtes pour déterminer le nombre total de nœuds
        std::istringstream iss(line);
        if (iss >> fromNode >> toNode)
        {
            numNodes = std::max({numNodes, fromNode, toNode});
        }
    }

    // Réinitialiser la position du curseur dans le fichier
    inputFile.clear();
    inputFile.seekg(0, std::ios::beg);

    numNodes += 1;

    // Initialiser la matrice d'adjacence avec des zéros
    double *A_data = new double[numNodes * numNodes];
    std::fill(A_data, A_data + numNodes * numNodes, 0.0);

    while (std::getline(inputFile, line))
    {
        int fromNode = 0, toNode = 0;
        std::istringstream iss(line);
        try
        {
            iss >> fromNode >> toNode;
            A_data[(fromNode * numNodes) + toNode] = 1;
        }
        catch (const std::exception &e)
        {
            continue;
        }
    }

    // Fermer le fichier
    inputFile.close();

    Matrix A(A_data, numNodes, numNodes);
    return A;
}

std::pair<double, Vector> power_iteration(Matrix &transitionMatrix, double alpha, double tol, int it_max, int &iterations)
{
    // Matrice
    int n = transitionMatrix.getRow();
    int m = transitionMatrix.getColumn();
    double *transitionMatrixData = transitionMatrix.getData();

    double beta = (1 - alpha) / n;
    double lambda = 0;
    iterations = 0;

    // Initialisation du Vecteur propre initial
    double *x_data = new double[n];
    Vector x(x_data, n);
    x.fillRandomValue();
    x.normalize(); // Normalisation du vecteur propre initial

    // Cree un vecteur v rempli de 1
    double *v_data = new double[n];
    std::fill(v_data, v_data + n, 1.0);
    Vector v(v_data, n);

    // Creer un vecteur y
    double *y_data = new double[n];
    std::fill(y_data, y_data + n, 0.0);
    Vector y(y_data, n);

    // Creer un vecteur r
    Vector r = x.minus(y);

    while (r.norm2() > tol && iterations < it_max)
    {
        double sum = x.sum();

        Vector y_buff = transitionMatrix.multiplyByVector(x, alpha);
        y = y_buff.plus(v.multiplyByScalar(beta * sum));

        lambda = x.dotProduct(y);

        r = y.minus(x);

        x = y.multiplyByScalar(1.0 / y.norm2());

        iterations++;
    }

    return std::make_pair(lambda, x);
}

int main(int argc, char *argv[])
{
    // double data[] = {0, 1, 1, 1, 0,
    //                  1, 0, 1, 0, 1,
    //                  0, 1, 0, 0, 1,
    //                  0, 0, 0, 0, 1,
    //                  0, 0, 1, 0, 0};

    double alpha = ALPHA;
    double tol = TOL;

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

    //***********************************************
    // Lire la matrice d'adjacence à partir du fichier
    //***********************************************
    std::string filename = argv[1];
    Matrix adjacencyMatrix = readMatrixFromFile(filename);

    //***********************************************
    // Construire la matrice de transition P
    //***********************************************
    Matrix transitionMatrixP = adjacencyMatrix.buildTransitionMatrix();

    //***********************************************
    // Appliquer la méthode de la puissance pour obtenir le vecteur de PageRank
    //***********************************************

    int n = adjacencyMatrix.getRow();
    int m = adjacencyMatrix.getColumn();
    double *adjacencyMatrixData = adjacencyMatrix.getData();
    int max_index = -1;
    double max_element = -1;
    int min_index = -1;
    double min_element = -1;
    int iterations = -1;

    // Mesurer le temps d'exécution de power_iteration
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::pair<double, Vector> page_rank = power_iteration(transitionMatrixP, alpha, tol, MAX_ITER, iterations);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // Calculer la durée d'exécution en millisecondes
    double execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    // Afficher le temps d'exécution
    double lambda = page_rank.first;
    Vector x = page_rank.second;
    x.findMinMaxIndex(max_index, max_element, min_index, min_element);

    // Affichage des résultats
    std::cout << "--- Résultats de l'Algorithme de PageRank --- \n"
              << "\nConfiguration : \n"
              << "\t-Taille de la matrice d'adjacence : " << n << "*" << m << "\n"
              << "\nParamètres de l'algorithme : \n"
              << "\t-Tolerance : " << tol << "\n"
              << "\t-Damping factor(alpha) : " << alpha << "\n"
              << "\nResultats : \n"
              << "Nombre d'itérations effectuées : " << iterations << "\n"
              << "Valeur propre dominante : " << lambda << "\n"
              << "Vecteur propre correspondant: \n"
              << "\t-Maximum value: " << max_element << " (Index : " << max_index << ")" << std::endl;
    //   << "\t-Minimum value: " << min_element << " (Index : " << min_index << ")" << std::endl;

    std::cout << "\nTemps d'exécution de power_iteration : " << execution_time << " ms" << std::endl;

    // Écriture du vecteur propre dans un fichier
    std::string eigenvector_filename = "eigenvector.txt";
    std::cout << "Ecriture du vecteur propre dans le fichier : " << eigenvector_filename << std::endl;
    x.writeDataToFile(eigenvector_filename);

    return 0;
}
