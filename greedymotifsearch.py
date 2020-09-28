#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The basic idea of the greedy motif search algorithm is to find the set of motifs across a number of DNA sequences that match each other  most closely.
"""
# GreedyMotifSearch(Dna, k, t)
#        BestMotifs ← motif matrix formed by first k-mers in each string from Dna
#        for each k-mer Motif in the first string from Dna
#            Motif1 ← Motif
#            for i = 2 to t
#                form Profile from motifs Motif1, …, Motifi - 1
#                Motifi ← Profile-most probable k-mer in the i-th string in Dna
#            Motifs ← (Motif1, …, Motift)
#            if Score(Motifs) < Score(BestMotifs)
#                BestMotifs ← Motifs
#        return BestMotifs

import profile_matrix as pm
import hammingdist as hd
import math
import os

def generate_profile(motifs):
    k = len(motifs[0])
    profile = {'A': [0] * k, 'C': [0] * k, 'G': [0] * k, 'T': [0] * k}
    div = float(len(motifs))
    for i in range(k):
        for motif in motifs:
            profile[motif[i]][i] += 1
        for key in profile:
            profile[key][i] /= div
    return profile

def window(s, k):
    for i in range(1 + len(s) - k):
        yield s[i:i+k]
        
def Score(motifs):
    """
    Finds score of motifs relative to the consensus sequence

    motifs: a list of given motifs (list)

    Returns: score of given motifs (int)
    """
    consensus = FindConsensus(motifs)
    score = 0.0000
    for motif in motifs:
        score += hd.hammingDistance(consensus, motif)
    #print(score)
    return round(score, 4)

def FindConsensus(motifs):
    """
    Finds a consensus sequence for given list of motifs

    motifs: a list of motif sequences (list)

    Returns: consensus sequence of motifs (str)
    """
    consensus = ""
    for i in range(len(motifs[0])):
        countA, countC, countG, countT = 0, 0, 0, 0
        for motif in motifs:
            if motif[i] == "A":
                countA += 1
            elif motif[i] == "C":
                countC += 1
            elif motif[i] == "G":
                countG += 1
            elif motif[i] == "T":
                countT += 1
        if countA >= max(countC, countG, countT):
            consensus += "A"
        elif countC >= max(countA, countG, countT):
            consensus += "C"
        elif countG >= max(countC, countA, countT):
            consensus += "G"
        elif countT >= max(countC, countG, countA):
            consensus += "T"
    return consensus

print(FindConsensus())

def GreedyMotifSearch(Dna, k, t):
    bestMotifs = [stringDna[:k] for stringDna in Dna]
    bestScore = math.inf
    base = Dna[0]
    
    for i in window(base, k):
        newMotifs = [i]
        for j in range(1, len(Dna)):
            profile = generate_profile(newMotifs)
            probable = pm.MostProbableProfile(Dna[j], k, profile)
            newMotifs.append(probable)

        if Score(newMotifs) < bestScore:
            bestScore = Score(newMotifs)
            bestMotifs = newMotifs
    return (' ').join(bestMotifs)

file = os.path.expanduser("~/Documents/Bella/Bioinformatics/dataset_159_5.txt")
data = open(file, "r")
reader = data.read().split()

sequence = reader[2:]
k = int(reader[0])
t = int(reader[1])
#indexes = []
#i = 0
#while i <= k*4-k:
#    indexes.append(i)
#    i += k
#    
#all_profiles_values = reader[2:]
#
#profileMatrixNew = {}
#
#nucleotides = ['A', 'C', 'G', 'T']
#
#j = 0
#for index_val in indexes:
#    profileMatrixNew[nucleotides[j]] = [float(x) for x in all_profiles_values[index_val:index_val+k]]
#    j += 1
    
#print(reader)
print(GreedyMotifSearch(sequence, k, t))