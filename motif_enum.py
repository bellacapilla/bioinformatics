#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:57:50 2020
@author: isabellac

In a collection of strings Dna and a int d, a k-mer is a (k,d)-motif if it appears 
in every string from Dna with at most d mismatches.

The function MotifEnumeration is a brute force (exhaustive search) problem-solving technique
that finds all (k, d)-motifs in a collection of strings by exploring all possible soltion candidates and
checking whether each candidate solves the problem. 

Input: collection of strings Dna, and integers k and d.
Output: All (k,d)-motifs in Dna

"""
import os
import kmers as km
import hammingdist as hd

def MotifEnumeration(Dna, k, d):
    Patterns = set()
    candidates = {}
    kmers = list(km.getKmersAndNeighborhood(Dna[0], k,d))

    for j in range(len(kmers)):
        pattern_one = kmers[j]
        
        for z in range(1, len(kmers)):
            pattern_two = kmers[z]
            hamming_value = hd.hammingDistance(pattern_one, pattern_two)
            
            if hamming_value <= d and  pattern_one not in candidates:
                candidates[pattern_one] = set({0})
                
    for candidate in candidates: 
        
        for i in range(1, len(Dna)):
            kmers_to_compare = km.getKmers(Dna[i], k)
            
            for kmer_two in kmers_to_compare:
                difference = hd.hammingDistance(candidate, kmer_two)
                
                if difference <= d:
                    candidates[candidate].add(i)
                    
    #print(candidates)            
    for candidate in candidates:
        if len(candidates[candidate]) == len(Dna):
            Patterns.add(candidate)
        
    return (' ').join(sorted(Patterns))

file = os.path.expanduser("~/Documents/Bella/Bioinformatics/dataset_156_8.txt")
data = open(file, "r")
reader = data.read().split()

Dna = []
for u in range(2, len(reader)):
    Dna.append(reader[u])
test = MotifEnumeration(Dna, int(reader[0]), int(reader[1]))
#    


#Dna = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
#k = 3
#d = 1
#test = MotifEnumeration(Dna, k, d)

#print(test)
