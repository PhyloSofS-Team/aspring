#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 12:58:02 2023

@author: diego
"""

#%%

#Exercice 1 :
#Écrire un programme qui affiche le triangle suivant en utilisant une boucle for :
#*
#**
#***
#****
#*****
#******


for i in range(1, 7):
    print('*' * i)

#%%

#Exercice 2 : Poly-A
#Créez une chaîne de caractères représentant un brin d'ADN poly-A (c'est-à-dire composé
#uniquement de bases A) de longueur 20 en utilisant une boucle for, sans avoir à entrer
#manuellement chaque base. Évitez d'utiliser la multiplication de caractères pour répéter la
#base A.
    
poly_a = ''

for i in range(20):
    poly_a += 'A'   

#%%

#Exercice 3 : Codons
#Donnée une chaîne de caractères représentant une séquence codante, par exemple
#'ATGCAGTAC', écrire un script qui affiche les codons contenus dans la séquence en
#utilisant une boucle for et la fonction range().
    
seq = 'ATGCAGTAC'
for i in range(0, len(seq), 3):
    print(seq[i:i+3])

#%%

#Exercice 4 : Noms
#Écrire un script pour convertir les noms dans une liste en majuscules.
#Entrée : Une liste de noms avec la première lettre en majuscules et le reste en minuscules
#(chaînes de caractères)
#Exemple d’entrée: ["Jacques", "Marie", "Dupont"]
#Sortie : Afficher les noms dans la liste en majuscules (chaînes de caractères)
#Exemple de sortie : ["JACQUES", "MARIE", "DUPONT"]
#Procédure :
#1. Pour chaque nom dans la liste de noms,
#a. Convertir chaque caractère du nom en majuscules en utilisant la méthode de
#conversion de chaîne de caractères appropriée
#b. Remplacer le nom dans la même position dans la liste avec le nom converti en
#majuscules
#2. Afficher chaque nom dans la liste en majuscules.

noms = ["Jacques", "Marie", "Dupont"]
for i in range(len(noms)):
    noms[i] = noms[i].upper()
print(noms)

#%%

#Exercice 5 : Afficher les nombres pairs
#Créez un script qui parcourt une liste de nombres entiers et affiche uniquement les nombres
#pairs. Exemple d’entrée : [1, 4, 6, 8, 9, 11]

#%%

#Exercice 6 : Séquence complémentaire d’un brin d’ADN
#La liste ci-dessous représente la séquence d’un brin d’ADN :
#["A", "C", "G", "T", "T", "A", "G", "C", "T", "A", "A", "C", "G"]
#Créez un script qui transforme cette séquence en sa séquence complémentaire. Rappel : la
#séquence complémentaire s’obtient en remplaçant A par T, T par A, C par G et G par C.

dna = ["A", "C", "G", "T", "T", "A", "G", "C", "T", "A", "A", "C", "G"]
for i in range(len(dna)):
    current_base = dna[i]
    if current_base == 'A':
        dna[i] = 'T'
    elif current_base == 'T':
        dna[i] = 'A'
    elif  current_base == 'C':
        dna[i] = 'G'
    else:
        dna[i] = 'C'

print(dna)

#%%

#Exercice 7 : Minimum d’une liste
#La fonction min() de Python renvoie l’élément le plus petit d’une liste constituée de valeurs
#numériques ou de chaînes de caractères. Sans utiliser cette fonction, créez un script qui
#détermine le plus petit élément de la liste [8, 4, 6, 1, 5].


liste = [8, 4, 6, 1, 5]

min_val = liste[0]

for val in liste[1:]:
    if val < min_val:
        min_val = val

print(min_val)

#%%

#Exercise 8 : Trouver le premier codon "ATG" dans une liste de codons
#La liste ci-dessous représente une séquence de codons d'ADN :
#["TAC", "CAT", "TGC", "ATG", "TTA", "ATA"]
#Créez un script qui trouve et affiche l'index du premier codon "ATG" dans cette liste en
#utilisant une boucle while. Quel est l'avantage d'utiliser while plutôt que for ?

codons = ["TAC", "CAT", "TGC", "ATG", "TTA", "ATA"]

index = 0
while index < len(codons) and codons[index] != "ATG":
    index += 1

print(index)

#%%

#Exercice 9 :
#Écrire une fonction qui permet de déterminer la température en celsius à partir de la
#température en fahrenheit selon la formule :
#temp_celsius = (temp_fahrenheit − 32) ×5/9

def fahrenheit_to_celsius(temp_fahrenheit):
    return (temp_fahrenheit - 32) * 5/9

fahrenheit_to_celsius(451)

#%%

#Exercice 10 :
#Créez une fonction calc_distance() qui calcule la distance euclidienne en trois dimensions
#entre deux atomes. Testez votre fonction sur les 2 points A [0, 0, 0] et B [1, 1, 1].
#Trouvez-vous bien √3 ?
#On rappelle que la distance euclidienne d entre deux points A et B de coordonnées
#cartésiennes respectives (xA, yA, zA) et (xB, yB, zB) se calcule comme suit :
#d = √((xB - xA) 2 + (yB - yA) 2 + (zB - zA) 2 )

import math

A = [0, 0, 0]
B = [1, 1, 1]

def calc_distance(A, B):
    return math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2 + (B[2] - A[2])**2)

print(calc_distance(A, B))


#%%

#Exercice 11 : Calculer le poids moléculaire d'une molécule
#Écrire une fonction qui permet de calculer le poids moléculaire d’une molécule à partir de
#sa formule brute.
#Entrée : Une formule brute d'une molécule sous forme de chaîne de caractères, par exemple
#"C2H5OH"
#Sortie : Le poids moléculaire de la molécule, par exemple 46 g/mol
#Description :
#Le poids moléculaire d'une molécule peut être calculé en utilisant le nombre d'atomes de
#chaque élément présent dans la molécule et leur poids atomique respectif. Pour cet exercice,
#nous nous limiterons aux éléments suivants et leurs poids atomiques respectifs :
#Hydrogène (H) : 1 g/mol
#Carbone (C) : 12 g/mol
#Azote (N) : 14 g/mol
#Oxygène (O) : 16 g/mol


def poids_element(elem):
    poids = 0
    if elem == 'C':
        poids = 12
    elif elem == 'H':
        poids = 1
    elif elem == 'N':
        poids = 14
    elif elem == 'O':
        poids = 16
    return poids

molecule = "C2H5OH"

def poids_molecule(molecule):
    molecule = list(molecule)
    n = len(molecule)
    poids = 0
    i = 0
    while i < n:
        poids_elem = poids_element(molecule[i])
        if poids_elem == 0:
            nombre = int(molecule[i])
            poids += poids_element(molecule[i - 1]) * (nombre - 1)
        else:
            poids += poids_elem
        i += 1
    return poids

print(poids_molecule(molecule))
