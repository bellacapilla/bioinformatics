#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:45:43 2020

@author: isabellac
"""

from hammingdist import hammingDistance
from kmers import getKmersAndComplement, getKmers

## Approximate Pattern Matching
# Returns the indexes of kmers matching with a maximum k of hamming distance within a sequence
def approximatePatternMatch(pattern, sequence, k):
    kmer_array = getKmers(sequence, len(pattern))
    indexes = []
        
    for kmer_index in range(len(kmer_array)):
        hamming_diff = hammingDistance(kmer_array[kmer_index], pattern)
        if hamming_diff <= int(k):
            indexes.append(kmer_index)
            
    return len(indexes)

## Get Frequent Matches 
def kmersFrequencyMatches(sequence, kmer_size, limit_hamming):
    kmers = getKmersAndComplement(sequence, kmer_size)
    matches = []
    
    
    for i in range(len(kmers)):
        kmer = kmers[i:i+kmer_size]
        
        for j in range(i, len(kmers)):
            match = kmers[j:j+kmer_size]
            mismatch = hammingDistance(kmer, match)
            
            if mismatch <= int(limit_hamming):
                matches.append(kmer)
                
#        for reversed_kmer in reversed_kmers:
#            mismatch = hammingDistance(reversed_kmer, pattern)
#            if mismatch <= int(limit_hamming):
#                matches.append(reversed_kmer)
    
    return matches


def mostFrequentKmerLimitHamming(text, kmer_size, mismatch):
    kmers_array = (getKmers(text, kmer_size))
    dic = {}
    
    for kmer in kmers_array: 
        print(kmer)
        freq = kmersFrequencyMatches(text, kmer, mismatch)
        for match in freq:
            if match not in dic:
                dic[match] = 0
            else:
                
                dic[match] = dic[match] + 1
        
    
    return dic