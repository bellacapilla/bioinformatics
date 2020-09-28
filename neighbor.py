#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:27:33 2020

@author: isabellac
"""

from hammingdist import hammingDistance

def Neighbors(Pattern, d):
    if d == 0:
        return Pattern
    elif len(Pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    Neighborhood = []
    SuffixNeighbors = Neighbors(Pattern[1:], d)
    for Text in SuffixNeighbors:
        nucleotides = {'A', 'C', 'G', 'T'}
        if hammingDistance(Pattern[1:], Text) < d:
            for x in nucleotides:
                Neighborhood.append(x + Text)
        else:
            Neighborhood.append(Pattern[0] + Text)
    return Neighborhood