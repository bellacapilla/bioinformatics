#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 13:29:48 2020

@author: isabellac
"""

import os
import kmers as km
import hammingdist as hd

def MedianString(Dna, k):
    distance = float('inf')
    kmers = km.arrangements(k)
    Median = ''
    
    for kmer in kmers:
        
        dist_compare = sum([min([hd.hammingDistance(Dna[i:i+k], kmer) for i in range(len(Dna) - k + 1)]) for dna in Dna])
        

        if distance > dist_compare:
            distance = dist_compare
            Median = kmer
                
    return Median
            
            
#file = os.path.expanduser("~/Documents/Bella/Bioinformatics/dataset_158_9.txt")
#data = open(file, "r")
#reader = data.read().split()
#Dna = []
#for u in range(1, len(reader)):
#    Dna.append(reader[u])
#test = MedianString(Dna, int(reader[0]))
    

#Dna = ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTTCGGGACAG']
#k = 3
#
#test = MedianString(Dna, k)
print(MedianString('CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCCGCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTCGGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG', 7))