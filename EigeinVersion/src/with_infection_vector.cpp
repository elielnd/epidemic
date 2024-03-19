#include "../include/page_rank.hpp"

using namespace Eigen;

std::vector<int> spreadEpidemicInfectionVector(Eigen::SparseMatrix<double> &socialGraph, Eigen::VectorXd infection_vector, double initialInfectedPercentage, double vaccinationPercentage, double nu, double alpha, double delta, double bellkis, int time)
{ 
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    int numIndividuals = socialGraph.rows();

    // Creer un vecteur pour enregistrer les individus guéris
    VectorXd curedIndividuals = VectorXd::Zero(numIndividuals);

    // Initialize infected individuals
    int numInitialInfected = std::ceil(numIndividuals * initialInfectedPercentage);

    for (int i = 0; i < numInitialInfected; i++)
    {
        int randomIndividual = std::uniform_int_distribution<>(0, numIndividuals - 1)(gen);
        socialGraph.coeffRef(randomIndividual, randomIndividual) = 1;
    }

    // Vaccinate individuals
    int numVaccinated = std::ceil(numIndividuals * vaccinationPercentage);

    // Trié le vecteur d'infection par ordre décroissant
    std::vector<std::pair<int, double>> infection_vector_sorted;
    for (int i = 0; i < numIndividuals; i++)
    {
        infection_vector_sorted.push_back(std::make_pair(i, infection_vector(i)));
    }
    std::sort(infection_vector_sorted.begin(), infection_vector_sorted.end(), [](const std::pair<int, double> &a, const std::pair<int, double> &b)
              { return a.second > b.second; });

    for (int i = 0; i < numVaccinated; i++)
    {

        socialGraph.coeffRef(infection_vector_sorted[i].first, infection_vector_sorted[i].first) = 0;
        curedIndividuals(infection_vector_sorted[i].first) = 1;
    }

    // Iterate
    int iterations = 0;

    // Creer un vecteur pour enregistrer le nombre d'individus infectés à chaque étape
    std::vector<int> infectedIndividuals;

    while (iterations < time)
    {
        std::cout << "Iteration: " << iterations << std::endl;
        Eigen::SparseMatrix<double> newSocialGraph = socialGraph;

        for (int i = 0; i < numIndividuals; i++)
        {
            for (int j = 0; j < numIndividuals; j++)
            {
                if (socialGraph.coeff(i, j) == 1 && i == j)
                {
                    // Individual infects each of his/her neighbors
                    for (int k = 0; k < numIndividuals; k++)
                    {
                        double random = dis(gen);
                        if ((socialGraph.coeff(j, k) == 1) && curedIndividuals[k] == 1 && random < bellkis && k != j)
                        {
                            newSocialGraph.coeffRef(k, k) = 1;
                        }

                        if (socialGraph.coeff(j, k) == 1 && curedIndividuals[k] == 0 && dis(gen) < nu && k != j)
                        {
                            newSocialGraph.coeffRef(k, k) = 1;
                        }
                    }

                    double random = dis(gen);
                    // Individual tries to infect a non-neighbor individual
                    if (random < (1 - alpha))
                    {
                        int randomIndividual = std::uniform_int_distribution<>(0, numIndividuals - 1)(gen);

                        if (curedIndividuals[randomIndividual] == 1 && dis(gen) < bellkis)
                            newSocialGraph.coeffRef(randomIndividual, randomIndividual) = 1;

                        if (curedIndividuals[randomIndividual] == 0)
                        {
                            newSocialGraph.coeffRef(randomIndividual, randomIndividual) = 1;
                        }
                    }

                    random = dis(gen);
                    // Individual is cured
                    if (random < delta)
                    {
                        curedIndividuals[i] = 1;
                        newSocialGraph.coeffRef(i, i) = 0;
                    }
                }
            }
        }


        // Enregistrer le nombre d'individus infectés
        int numFinalInfected = newSocialGraph.diagonal().sum();
        infectedIndividuals.push_back(numFinalInfected);

        socialGraph = newSocialGraph;
        iterations++;
    }

    return infectedIndividuals;
}

int main(int argc, char *argv[])
{
    
    double alpha = ALPHA;
    double tol = TOL;
     double nu = NU;
    double delta = DELTA;
    double bellkis = BELLKIS;
    double initialInfectedPercentage = INITIAL_INFECTED_PERCENTAGE;
    double vaccinationPercentage = VACCINATION_PERCENTAGE;
    int time = TIME;
    std::string outputFile = USING_INFECTED_VECTOR_FILE;
    std::string eigenvector_filename = INF_RANK;

    // Vérifier le nombre d'arguments
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <adjacency_matrix_file>" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (argc >= 3)
    {
        outputFile = argv[2];
    }

    outputFile = OUTPUT_DIR + outputFile;

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
    auto result = powerIteration(transitionMatrixP, alpha, tol, MAX_ITER, iter);

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

    //***********************************************
    // Spread epidemic without vaccination
    //***********************************************

    int run = RUN;

    std::vector<std::vector<int>> infectedIndividualsPerTimeStep(run);

    for (int r = 0; r < run; r++)
    {
        std::vector<int> infectedIndividualsWithoutVaccination = spreadEpidemicInfectionVector(adjacencyMatrix, x, initialInfectedPercentage, vaccinationPercentage, nu, alpha, delta, bellkis, time);
        infectedIndividualsPerTimeStep[r] = infectedIndividualsWithoutVaccination;
    }

    std::vector<int> averageInfectedIndividuals(time, 0);

    for (int t = 0; t < time; t++)
    {
        for (int r = 0; r < run; r++)
        {
            averageInfectedIndividuals[t] += infectedIndividualsPerTimeStep[r][t];
        }
        averageInfectedIndividuals[t] /= run;
    }

    std::cout << "\n**********************************" << std::endl;
    std::cout << "*** Spread epidemic using the infection vector ***" << std::endl;

    std::cout << "\n----------------------------------" << std::endl;
    std::cout << "--- Configuration: ---\n"
              << "\t-Time : " << time << "\n"
              << "\t-Alpha : " << alpha << "\n"
              << "\t-Nu : " << nu << "\n"
              << "\t-Delta : " << delta << "\n"
              << "\t-RUN : " << run << "\n"
              << "\t-bellkis : " << bellkis << "\n"
              << "\t-Initial Infected Percentage : " << initialInfectedPercentage << "\n"
              << "\t-Vaccination Percentage : " << vaccinationPercentage << "\n"
              << std::endl;

    std::cout << "----------------------------------\n"
              << std::endl;

    // Enregistrer le nombre d'individus infectés à chaque étape dans un fichier csv
    std::ofstream file(outputFile);
    file << "Time,Infected Individuals" << std::endl;
    for (int i = 0; i < averageInfectedIndividuals.size(); i++)
    {
        file << i << "," << averageInfectedIndividuals[i] << std::endl;
    }


    std::cout << "... Writing to file " << outputFile << " ... \n";
    std::cout << "Results have been saved in: " << outputFile << std::endl;

    return 0;
}