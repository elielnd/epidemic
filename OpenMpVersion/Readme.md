# PageRank
Implémentation de l’algorithme de PageRank avec la méthode des puissance en C++


## Table des matières

1. [Introduction](#introduction)
2. [Dépendances](#dépendances)
3. [Installation](#installation)
4. [Utilisation](#utilisation)
5. [Exemples](#exemples)


## Introduction
L'algorithme de PageRank est une technique d'analyse de liens qui assigne à chaque page web un score de pertinence numérique. L'idée fondamentale derrière PageRank est que les pages importantes sont celles qui sont liées par d'autres pages importantes. En d'autres termes, une page est considérée comme importante si elle est pointée par d'autres pages importantes.

L'implémentation de cet algorithme utilisera la méthode des puissances pour calculer les scores de PageRank. La méthode des puissances est une technique itérative qui converge vers le vecteur propre principal de la matrice d'adjacence du graphe web.

## Dépendances

Ce projet dépend de OpenMP, il faut vous rassurer que la librairies soit présent avant de compiler le 
projet


## Installation

make

## Utilisation 
./page_rank DATA ALPHA TOL

DATA: Chemin vers les données à utiliser, un exemple de format de données est contenu dans le dossier data/

ALPHA: Damping factor

TOL: Tolérance 

exemple : 
    ./page_rank "data/email-Eu-core.txt" 0.85 1e-6


## Exemples
un exemple de sortie est le suivant :

--- Résultats de l'Algorithme de PageRank --- 

Configuration : 
- Taille de la matrice d'adjacence : 1005*1005
- Damping factor(alpha) : 0.85

Nombre d'itérations effectuées : 56
Tolerance : 1e-06

Resultats : 
Valeur propre dominante : 0.997827
Vecteur propre correspondant: 
        -Maximum value: 0.244645 (Index : 160)
        -Minimum value: 0.00324071 (Index : 78)
Données du vecteur propre écrites avec succès dans le fichier: eigenvector.txt



### Vecteur propre 
Un fichier contenant le vecteur propre est crée 


