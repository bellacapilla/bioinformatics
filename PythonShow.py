# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 12:49:19 2020

@author: Isabella Capilla

The Python Show
Exploring the DNA with Python
"""
import os


## Basic Sequence Slide 6
sequence = ("atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc").upper() 
kmer_example_one = "atgatcaag"





## Get possible kmers
def getKmers(sequence, kmer_size):
    kmer_array = []
    for j in range(len(sequence)-kmer_size+1):
        kmer_array.append(sequence[j:j+kmer_size])
    return kmer_array

print(getKmers(sequence, 3))
    








################################    
# Generate the reverse complements of "possible kmers"
nucs = { 'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}

RCDict = {}
RCList = []

##########################
## Get possible kmers and its complement

def RevComplement(kmerList):
    for kmer in kmerList:
        c = ''
        for letter in kmer:
            c += nucs[letter]
        rc = c[::-1]
        RCList.append(rc)
        RCDict[kmer] = rc
    return(RCDict)

#print(RevComplement(getKmers(sequence, 9)))
    







##############################
### Hamming Distance Calculation
def hammingDistance(pattern_one, pattern_two):
    count = 0
    for i in range(len(pattern_one)):
        if pattern_one[i] != pattern_two[i]:
            count += 1
    return count








## Count reappearances of kmer 
def Count_pattern(text, k, d):
    
    # Generate the "possible kmers"
    kmers = getKmers(sequence, k)
    #print(kmers)

    
    countDict = {}  # initiate dictionary 
    
    RCDict = RevComplement(kmers) # dictionary with kmer and rc pairs
    kmerList = list(RCDict.keys())  # list of kmers in order
    rcList = list(RCDict.values())  # list of rc's in corresponding order

    # Initializate counts as zero
    for kmer in kmerList:
        countDict[kmer.upper()] = 0
        
    
    for i in range(len(kmerList)):
        for j in range(len(text)-k+1):
            word = text[j:j+k]
            
            
            # kmer matching and sum count
            if (word == kmers[i]):
                countDict[kmers[i]] += 1
                
            #kmer and kmer mismatches
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
    returnList = []
    

    for key in countDict.keys():
        if maxCount == 0:
            return('no frequent kmer')
        else:
            if countDict[key] == maxCount:
                returnKmer = key
                returnRC = RCDict[key]
                #if returnKmer not in returnList:
                if returnKmer not in returnList and returnRC not in returnList:
                    returnList.append(returnKmer)
                    returnList.append(returnRC)
                    
    
    return(returnList)
    
print(Count_pattern(sequence.upper(), 5, 0))





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

print((skewCalc(sequence)))
print(getMinimumPoint(sequence))
print(len(sequence))



