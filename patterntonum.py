#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:51:59 2020

@author: isabellac
"""

map_dic = {
        'A':0,
        'C':1,
        'G':2,
        'T':3}
def patternToNumber(pattern_dna):
    length_pattern = 1
    calculation = 0
    for elem in pattern_dna:
        if elem in map_dic:
            calculation += (map_dic[elem]) * (4**(len(pattern_dna)-length_pattern))
            length_pattern += 1
    return calculation