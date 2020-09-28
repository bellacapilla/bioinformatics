#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:25:33 2020

@author: isabellac
"""
from itertools import product
#from neighbor import Neighbors

## Get kmers
def getKmers(sequence, kmer_size):
    kmer_array = []
    for j in range(len(sequence)-kmer_size+1):
        kmer_array.append(sequence[j:j+kmer_size])
    return kmer_array

## Get possible kmers and its complement
def getKmersAndComplement(sequence, kmer_size):
    kmer_array = []
    for j in range(len(sequence)-kmer_size+1):
        kmer = sequence[j:j+kmer_size]
        kmer_array.append(kmer)
        kmer_array.append(kmer[::-1])
    return kmer_array

def getKmersAndNeighborhood(sequence, kmer_size, mismatch):
    kmer_array = []
    for j in range(len(sequence)-kmer_size+1):
        kmer_array.append(Neighbors(sequence[j:j+kmer_size], mismatch))
        kmer_array.append([sequence[j:j+kmer_size]])
    
    # Unifies list of lists into single list
    kmer_array = [item for sublist in kmer_array for item in sublist]
    return set(kmer_array)

def arrangements(k):
    result = product('ACGT', repeat = k)
    result = map(list, result)
    words = []
    for item in result:
        word = ''
        for letter in item:
            word += letter
        words.append(word)
        
    return words

print(arrangements(8))
    
# Generate the reverse complements of "possible kmers"
nucs = { 'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}

RCDict = {}
RCList = []

def RevComplement(kmerList):
    for kmer in kmerList:
        c = ''
        for letter in kmer:
            c += nucs[letter]
        rc = c[::-1]
        RCList.append(rc)
        RCDict[kmer] = rc
    return(RCDict)
    