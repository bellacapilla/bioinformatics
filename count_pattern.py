#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:35:20 2020

@author: isabellac
"""
import os
from hammingdist import hammingDistance
import kmers as km

sequence = "atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"
  
# Count how many times each possible kmer and its reverse complement appear within Text
def Count_pattern(text, k, d):
    # Generate the "possible kmers"
    kmers = km.arrangements(k)
    
    countDict = {}  # initiate dictionary to count kmers and reverse complements (rc) found, as well as their mismatches
    RCDict = km.RevComplement(kmers) # dictionary with kmer and rc pairs
    kmerList = list(RCDict.keys())  # list of kmers in order
    rcList = list(RCDict.values())  # list of rc's in corresponding order

    for kmer in kmerList:
        countDict[kmer] = 0
    for i in range(len(kmerList)):
        for j in range(len(text)-k+1):
            word = text[j:j+k]
            
            # kmer and kmer mismatches
            hamming1 = hammingDistance(word, kmerList[i])
            if hamming1 <= d:
                forCount1 = (kmerList[i])
                countDict[forCount1] += 1

            # rc and rc mismatches
            hamming2 = hammingDistance(word, rcList[i])
            if hamming2 <= d:
                forCount2 = (kmerList[i])
                countDict[forCount2] += 1
                
    # most frequent kmer and rc
    maxCount = max(countDict.values())
    returnPair = []
    for key in countDict.keys():
        if maxCount == 0:
            return('no frequent kmer or reverse complement!')
        else:
            if countDict[key] == maxCount:
                returnKmer = key
                returnRC = RCDict[key]
                if returnKmer not in returnPair and returnRC not in returnPair:
                    returnPair.append(returnKmer)
                    returnPair.append(returnRC)
                    
    return returnPair
    
print(Count_pattern(sequence,3, 1))