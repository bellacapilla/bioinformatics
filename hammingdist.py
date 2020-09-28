#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:30:08 2020

@author: isabellac
"""

### Hamming Distance Calculation
def hammingDistance(pattern_one, pattern_two):
    count = 0
    for i in range(len(pattern_one)):
        if pattern_one[i] != pattern_two[i]:
            count += 1
    return count