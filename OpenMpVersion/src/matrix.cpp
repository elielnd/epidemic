#include "../include/matrix.hpp"

//***************
// Matrix
//***************
Matrix::Matrix(double *data, int n, int m)
{
    this->data = data;
    this->n = n;
    this->m = m;
}

Vector Matrix::multiplyByVector(Vector vector, double alpha)
{
    assert(m == vector.getRow());
    double *result = new double[n];

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        result[i] = 0;
#pragma omp parallel for reduction(+ : result[i])
        for (int j = 0; j < m; ++j)
        {
            result[i] += alpha * data[i * m + j] * vector.getData()[j];
        }
    }

    Vector vectorRes(result, n);
    return vectorRes;
}

double Matrix::sum()
{
    double sum = 0;
#pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < n * m; ++i)
    {
        sum += data[i];
    }
    return sum;
}

Matrix::~Matrix()
{
}

Matrix Matrix::transpose()
{
    double *transposedMatrix = new double[n * m];

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; j++)
        {
            transposedMatrix[j * n + i] = data[i * m + j];
        }
    }
    Matrix matrixRes(transposedMatrix, m, n);
    return matrixRes;
}

double Matrix::columnSum(int j)
{
    double sum = 0;
#pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < n; ++i)
    {
        sum += data[i * m + j];
    }
    return sum;
}

void Matrix::multiplyColumnByScalar(int j, double scalar)
{
// #pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        data[i * m + j] *= scalar;
    }
}

Matrix Matrix::buildTransitionMatrix()
{
    // Transposer la matrice
    Matrix transposed = this->transpose();
    double *transposedData = transposed.getData();

    for (int i = 0; i < m; ++i)
    {
        // Calculer la somme de la ligne avec CBLAS
        double colSum = transposed.columnSum(i);

        if (colSum != 0.0)
        {
            // cblas_dscal(n, 1/colSum, transposedData + i, m);
            // Multiplier la colonne par 1/colSum
            transposed.multiplyColumnByScalar(i, 1 / colSum);
        }
    }

    return transposed;
}

int Matrix::getRow()
{
    return n;
}

int Matrix::getColumn()
{
    return m;
}

double *Matrix::getData()
{
    return data;
}

void Matrix::setData(double *data)
{
    this->data = data;
}

// Méthode pour afficher le vecteur
void Matrix::print() const
{
    std::cout << "Martix: \n";
    for (int i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; j++)
        {
            std::cout << data[i * n + j] << " ";
        }

        std::cout << "\n";
    }
    std::cout << "\n";
}

//***************
// Vector
//***************
Vector::Vector(double *data, int n)
{
    this->data = data;
    this->n = n;
}

Vector::~Vector()
{
}

void Vector::fillRandomValue()
{
    for (int i = 0; i < n; ++i)
    {
        data[i] = static_cast<double>(rand()) / RAND_MAX;
    }
}

double Vector::norm2() const
{
    double result = 0.0;

// Utilisation de la directive parallel for reduction(+:result) pour paralléliser la sommation
#pragma omp parallel for reduction(+ : result)
    for (int i = 0; i < n; ++i)
    {
        result += data[i] * data[i];
    }

    return std::sqrt(result);
}

void Vector::normalize()
{
    double norm = this->norm2();

    if (norm != 0.0)
    {
        multiplyByScalar(1.0 / norm);
    }
}



Vector Vector::multiplyByScalar(double alpha)
{
    double *result = new double[n];

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        result[i] = alpha * data[i];
    }

    Vector vectorRes(result, n);
    return vectorRes;
}

Vector Vector::minus(Vector vector)
{
    assert(n == vector.getRow());
    double *result = new double[n];

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        result[i] = data[i] - vector.getData()[i];
    }

    Vector vectorRes(result, n);
    return vectorRes;
}

Vector Vector::plus(Vector vector)
{
    assert(n == vector.getRow());
    double *result = new double[n];

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        result[i] = data[i] + vector.getData()[i];
    }

    Vector vectorRes(result, n);
    return vectorRes;
}

double Vector::dotProduct(Vector vector)
{
    assert(n == vector.getRow());
    double dotProduct = 0;

#pragma omp parallel for reduction(+:dotProduct)
    for (int i = 0; i < n; ++i)
    {
        dotProduct += data[i] * vector.getData()[i];
    }

    return dotProduct;
}

double Vector::sum()
{
    double sum = 0;
#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; ++i)
    {
        sum += data[i];
    }
    return sum;
}

void Vector::findMinMaxIndex(int &max_index, double &max_element, int &min_index, double &min_element)
{
    max_index = -1;
    min_index = -1;
    max_element = -1;
    min_element = -1;

    // std::max_element pour obtenir un itérateur pointant vers l'élément maximum
    auto min_max_iterator = std::minmax_element(data, data + n);

    min_element = *min_max_iterator.first;
    max_element = *min_max_iterator.second;

    // Obtenir l'indice du maximum en soustrayant l'itérateur du début du tableau.

    min_index = std::distance(data, min_max_iterator.first);
    max_index = std::distance(data, min_max_iterator.second);
}

void Vector::writeDataToFile(const std::string &filename)
{
    // Ouvrir un fichier en écriture
    std::ofstream outputFile(filename);

    // Vérifier si l'ouverture du fichier a réussi
    if (outputFile.is_open())
    {
        // Écrire chaque élément du tableau dans le fichier, un élément par ligne
        for (int i = 0; i < n; ++i)
        {
            outputFile << data[i] << std::endl;
        }

        // Fermer le fichier
        outputFile.close();
        // std::cout << "Données du vecteur propre écrites avec succès dans le fichier: eigenvector.txt" << std::endl;
    }
    else
    {
        // Si l'ouverture du fichier échoue, afficher un message d'erreur
        std::cerr << "Erreur : Impossible d'ouvrir le fichier." << std::endl;
    }
}

// Méthode pour afficher le vecteur
void Vector::print() const
{
    std::cout << "Vector: ";
    for (int i = 0; i < n; ++i)
    {
        std::cout << data[i] << " ";
    }
    std::cout << "\n";
}
int Vector::getRow()
{
    return n;
}

double *Vector::getData()
{
    return data;
}

void Vector::setData(double *data)
{
    this->data = data;
}