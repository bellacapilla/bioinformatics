#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:57:02 2020

@author: isabellac
"""
import os
import motif_enum as mn

## Open file for exercise
file = os.path.expanduser("~/Documents/Bella/Bioinformatics/dataset_9_8.txt")
data = open(file, "r")
reader = data.read().split()

Dna = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
k = 3
d = 1

test = mn.MotifEnumeration(Dna, k, d)
print(test)




