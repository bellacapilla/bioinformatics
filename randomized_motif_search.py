#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 09:02:15 2020

@author: isabellac

RandomizedMotifSearch(Dna, k, t)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
        BestMotifs ← Motifs
        while forever
            Profile ← Profile(Motifs)
            Motifs ← Motifs(Profile, Dna)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
            else
                return BestMotifs
                
"""
import os
import random

N = 1000

def RandomizedMotifSearch(Dna, k, t):
    M = random_motifs(Dna, k, t)
    bestMotifs = M
    while True:
        profile = profile_with_pseudocounts(M)
        M = _motifs(profile, Dna)
        if _score(M) < _score(bestMotifs):
            bestMotifs = M
        else:
            return bestMotifs
        
def random_motifs(Dna, k, t):
    randMotifs = []

    for i in range(t):
        x = random.randint(0, t)
        randMotifs.append(Dna[i][x:x+k])

    return randMotifs
        
def profile_with_pseudocounts(motifs):
    profile = {}
    t = len(motifs)
#    k = len(motifs[0])
    countMotifs = count_with_pseudocounts(motifs)

    for symbol in "ACGT":
        profile[symbol] = []

    for x in countMotifs:
        for y in countMotifs[x]:
            z = y/float(t+4)
            profile[x].append(z)

    return profile

def count_with_pseudocounts(motifs):
    count = {}
    pseudocounts = {}
    t = len(motifs)
    k = len(motifs[0])

    for symbol in "GACT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1

    for symbol in "GACT":
        pseudocounts[symbol] = []

    for x in count:
        for y in count[x]:
            z = y + 1
            pseudocounts[x].append(z)

    return pseudocounts

def _score(motifs):
    count = 0
    k = len(motifs[0])
    t = len(motifs)
    consensusMotif = _consensus(motifs)

    for i in range(t):
        for j in range(k):
            if motifs[i][j] != consensusMotif[j]:
                count += 1

    return count

def _motifs(profile, Dna):
    motifs = []
    t = len(Dna)
    k = len(profile['A'])

    for i in range(t):
        motifs.append(profile_most_probable_kmer(Dna[i], k, profile))

    return motifs

def _consensus(motifs):
    k = len(motifs[0])
    count = count_with_pseudocounts(motifs)
    consensus = ""

    for j in range(k):
        M = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > M:
                M = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus

def profile_most_probable_kmer(text, k, profile):
    mostProbVal = -1
    mostProbKmer = ''

    for i in range(0, 1 + len(text) - k):
        kmer = text[i:i+k]
        probKmerVal = _pr(kmer, profile)
        if probKmerVal > mostProbVal:
            mostProbVal = probKmerVal
            mostProbKmer = kmer

    return mostProbKmer

def _pr(text, profile):
    P = 1

    for i in range(len(text)):
        P = P * profile[text[i]][i]

    return P




file = os.path.expanduser("~/Documents/Bella/Bioinformatics/dataset_161_5.txt")
data = open(file, "r")
reader = data.read().split()

Dna = reader[2:]
k = int(reader[0])
t = int(reader[1])

#Dna=['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
#k=8
#t=5
M = RandomizedMotifSearch(Dna, k, t)
bMotifs = M

for i in range(N+1):
    M = RandomizedMotifSearch(Dna, k, t)
    if _score(M) < _score(bMotifs):
         bMotifs = M
    else:
        bestMotifs = bMotifs

print ('\n'.join(bestMotifs))

#print(RandomizedMotifSearch(Dna, k, t))