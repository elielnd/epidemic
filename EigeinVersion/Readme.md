# PageRank
Implémentation de l’algorithme de PageRank avec la méthode des puissance en C++ en 


## Table des matières

1. [Introduction](#introduction)
2. [Dépendances](#dépendances)
3. [Installation](#installation)
4. [Utilisation](#utilisation)
5. [Exemples](#exemples)


## Introduction
L'algorithme de PageRank est une méthode d'analyse de liens qui attribue à chaque page web un score de pertinence numérique. Son principe est que les pages importantes sont celles qui sont liées par d'autres pages importantes. Cette implémentation utilise la méthode des puissances pour calculer les scores de PageRank, convergent vers le vecteur propre principal de la matrice d'adjacence du graphe web.

Dans ce projet, l'algorithme de PageRank est adapté pour simuler la propagation d'épidémies. Trois versions sont fournies : une sans vaccination, une avec vaccination, et une autre utilisant un vecteur d'infection.

## Dépendances
Ce projet dépend de la bibliothèque Eigen, fournie dans le package et ne nécessitant normalement aucune autre installation.

## Installation
Exécutez la commande `make`.

## Utilisation 
Utilisez la commande suivante pour exécuter le programme :
./page_rank DATA FILE

- DATA : Chemin vers les données à utiliser. Un exemple de format de données est contenu dans le dossier `data/`.
- FILE : (optionnel) Fichier de sortie.

Exemple : 

./page_rank "data/email-Eu-core.txt" 

./page_rank DATA ALPHA TOL

DATA: Chemin vers les données à utiliser, un exemple de format de données est contenu dans le dossier data/

TOL: Tolérance 

ALPHA: Damping factor

exemple : 
    ./page_rank "data/email-Eu-core.txt" 0.85 1e-08



## Exemples
Voici un exemple de sortie en exécutant la commande :

./without_vaccine ../data/Wiki-Vote.txt wiki_without_vaccine.csv

**********************************
*** Spread epidemic without vaccination ***

----------------------------------
--- Configuration: ---
        -Time : 100
        -Alpha : 0.85
        -Nu : 0.2
        -Delta : 0.24
        -RUN : 1
        -prob d'intection après guérison : 0.01
        -Initial Infected Percentage : 0.2
        -Vaccination Percentage : 0

----------------------------------

... Writing to file output/wiki_without_vaccine.csv ... 
Results have been saved in: output/wiki_without_vaccine.csv