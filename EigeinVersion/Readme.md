# PageRank
Implémentation de l’algorithme de PageRank avec la méthode des puissance en C++ en 


## Table des matières

1. [Introduction](#introduction)
2. [Dépendances](#dépendances)
3. [Installation](#installation)
4. [Utilisation](#utilisation)
5. [Exemples](#exemples)


## Introduction
L'algorithme de PageRank est une technique d'analyse de liens qui assigne à chaque page web un score de pertinence numérique. L'idée fondamentale derrière PageRank est que les pages importantes sont celles qui sont liées par d'autres pages importantes. En d'autres termes, une page est considérée comme importante si elle est pointée par d'autres pages importantes.

L'implémentation de cet algorithme utilisera la méthode des puissances pour calculer les scores de PageRank. La méthode des puissances est une technique itérative qui converge vers le vecteur propre principal de la matrice d'adjacence du graphe web.


## Installation

make

## Utilisation 
./page_rank DATA ALPHA TOL

DATA: Chemin vers les données à utiliser, un exemple de format de données est contenu dans le dossier data/

TOL: Tolérance 

ALPHA: Damping factor

exemple : 
    ./page_rank "data/email-Eu-core.txt" 0.85 1e-08


## Exemples
un exemple de sortie est le suivant :

--- Résultats de l'Algorithme de PageRank --- 

Configuration : 
        -Taille de la matrice d'adjacence : 1005*1005

Paramètres de l'algorithme : 
        -Tolerance : 1e-06
        -Damping factor(alpha) : 0.85

Resultats : 
Nombre d'itérations effectuées : 2000
Valeur propre dominante : 0.969566
Vecteur propre correspondant: 
        -Maximum value: 0.276819 (Index : 1)
Ecriture du vecteur propre dans le fichier : eigenvector.txt


### Vecteur propre 
Un fichier contenant le vecteur propre est crée 


