#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cassert>

class Vector;

class Matrix
{
private:
    double *data; // Tableau de données 1D
    int n;        // Nombre de lignes
    int m;        // Nombre de colonnes

public:
    /**
     * @brief Constructeur de la classe Matrix
     *
     * @param data Tableau de données 1D
     * @param n Nombre de lignes
     * @param m Nombre de colonnes
     */
    Matrix(double *data, int n, int m);

    /**
     * @brief Get the Data object
     *
     * @return double*
     */
    double *getData();

    /**
     * @brief Set the Data object
     *
     * @param data
     */
    void setData(double *data);

    /**
     * @brief Get the Row object
     *
     * @return int
     */
    int getRow();

    /**
     * @brief Get the Column object
     *
     * @return int
     */
    int getColumn();

    /**
     * @brief Multiplication de la matrice par un vecteur et un scalaire : alpha * A * vector
     *
     * @param vector Vecteur
     * @param alpha scalaire
     * @return Vector
     */
    Vector multiplyByVector(Vector vector, double alpha);

    /**
     * @brief Somme de tous les éléments de la matrice
     *
     * @return double
     */
    double sum();

    /**
     * @brief Calcul la Transposé de la matrice
     *
     * @return Matrix
     */
    Matrix transpose();

    /**
     * @brief Affiche la matrice
     *
     */
    void print() const;

    /**
     * @brief Construit la matrice de transition
     *
     * @return Matrix
     */
    Matrix buildTransitionMatrix();

    /**
     * @brief Calcul la somme des éléments de la colonne j
     *
     * @param j Numéro de la colonne
     * @return double
     */
    double columnSum(int j);

    /**
     * @brief Multiplie la colonne j par un scalaire
     *
     * @param j  Numéro de la colonne
     * @param scalar
     */
    void multiplyColumnByScalar(int j, double scalar);

    /**
     * @brief Destructeur de la classe Matrix
     *
     */
    ~Matrix();
};

class Vector
{
private:
    double *data; // Tableau de données 1D
    int n;        // Nombre de lignes

public:
    /**
     * @brief Constructeur de la classe Vector
     *
     * @param data Tableau de données 1D
     * @param n Nombre de lignes
     */
    Vector(double *data, int n);

    /**
     * @brief Get the Row object
     *
     * @return int
     */
    int getRow();

    /**
     * @brief Set the Row object
     *
     * @param n
     */
    void setRow(int n);

    /**
     * @brief Get the Data object
     *
     * @return double*
     */
    double *getData();

    /**
     * @brief Set the Data object
     *
     * @param data
     */
    void setData(double *data);

    /**
     * @brief Normalise le vecteur
     *
     */
    void normalize();

    /**
     * @brief Calcul la norme euclidienne du vecteur
     *
     * @return double
     */
    double norm2() const;
    // double Vector::norm2() const

    /**
     * @brief Calcul la somme des éléments du vecteur
     *
     * @return double
     */
    double sum();

    /**
     * @brief Soustraction de deux vecteurs
     *
     * @param vector
     * @return Vector
     */
    Vector minus(Vector vector);

    /**
     * @brief Addition de deux vecteurs
     *
     * @param vector
     * @return Vector
     */
    Vector plus(Vector vector);

    /**
     * @brief Multiplication de vecteur par un scalaire
     *
     * @param alpha
     * @return Vector
     */
    Vector multiplyByScalar(double alpha);

    /**
     * @brief Produit scalaire de deux vecteurs
     *
     * @param vector
     * @return double
     */
    double dotProduct(Vector vector);

    /**
     * @brief Trouve l'index et la valeur maximale et minimale du vecteur
     *
     * @param max_index index de la valeur maximale
     * @param max_element valeur maximale
     * @param min_index index de la valeur minimale
     * @param min_element valeur minimale
     */
    void findMinMaxIndex(int &max_index, double &max_element, int &min_index, double &min_element);

    /**
     * @brief Rempli le vecteur avec des valeurs aléatoires
     *
     */
    void fillRandomValue();

    /**
     * @brief Ecrit les données du vecteur dans un fichier
     *
     * @param filename Chemin du fichier
     */
    void writeDataToFile(const std::string &filename);

    /**
     * @brief Affiche le vecteur
     *
     */
    void print() const;

    /**
     * @brief Destructeur de la classe Vector
     *
     */
    ~Vector();
};
