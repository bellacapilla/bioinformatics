#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:41:12 2020

@author: isabellac
"""
import os
import operator
import kmers as km

profileMatrix = {
        'A': [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
        'C': [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
        'G': [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
        'T': [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]
        
        }

def ProfileMatrixCalc(pattern, matrixProfile):
    sum = 1
    for i in range(len(pattern)):
        sum *= matrixProfile[pattern[i]][i]
        
    return sum

print(ProfileMatrixCalc('TCGGTA', profileMatrix))


def MostProbableProfile(sequence, size, matrix):
    kmers = km.getKmers(sequence, size)
    results = {}
    
    for kmer in kmers:
        
        results[kmer] = ProfileMatrixCalc(kmer, matrix)
    
    return max(results.items(), key=operator.itemgetter(1))[0]
    

file = os.path.expanduser("~/Documents/Bella/Bioinformatics/dataset_159_3.txt")
data = open(file, "r")
reader = data.read().split()

#print(MostProbableProfile('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, profileMatrix))
sequence = reader[0]
k = int(reader[1])
indexes = []
i = 0
while i <= k*4-k:
    indexes.append(i)
    i += k
    
all_profiles_values = reader[2:]

profileMatrixNew = {}

nucleotides = ['A', 'C', 'G', 'T']

j = 0
for index_val in indexes:
    profileMatrixNew[nucleotides[j]] = [float(x) for x in all_profiles_values[index_val:index_val+k]]
    j += 1

test = MostProbableProfile(sequence, k, profileMatrixNew)
#print(test)