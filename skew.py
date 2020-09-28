#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:34:08 2020

@author: isabellac
"""
sequence = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
 
## Get Skew Genome Points
def skewCalc(sequence):
    count = 0
    result = [0]
    for nucleotide in sequence:
        if nucleotide == 'C':
            count -= 1
        elif nucleotide == 'G':
            count += 1
        
        result.append(count)
    return result


## Get min points of Skew Genome
def getMinimumPoint(pattern):
    skew_genome = skewCalc(pattern)
    min_value = min(skew_genome)
    indexes = []
    
    for j in range(len(skew_genome)-1):
        if skew_genome[j] == min_value:
            indexes.append(j)
    return indexes

print(getMinimumPoint(sequence))