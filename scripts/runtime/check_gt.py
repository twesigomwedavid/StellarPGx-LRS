#!/usr/bin/env python3


import os
import sys

infile = sys.argv[1]

# outfile = sys.argv[2]

list1 = []

for line in open(infile, "r"):

    list2 = []
    
    if line.startswith("#"):
        print(line.strip())
        
    else:
        line = line.strip().split()
        list1.append(line)

        
for i in list1:
    
    if i[6] == "RefCall": # and float(i[5]) > 0.1:
        i[6] = "PASS"
        gt_field = i[9].split(":")
        ad = gt_field[3].split(",")
        dp = int(ad[0]) + int(ad[1])
        
        if 0.25 < int(ad[1])/dp < 0.75:
            gt_field[0] = "0/1"
            
        elif int(ad[1])/int(dp) > 0.75:
            gt_field[0] = "1/1"

        else:
            gt_field[0] = "0/0"
            
        gt_field = ":".join(gt_field)

        i[9] = gt_field
        
        print("\t".join(i))

        
    elif i[6] == "PASS":

        gt_field = i[9].split(":")
        ad = gt_field[3].split(",")
        dp = int(ad[0]) + int(ad[1])

        if 0.25 < int(ad[1])/dp < 0.75:
            gt_field[0] = "0/1"

        elif int(ad[1])/int(dp) > 0.75:
            gt_field[0] = "1/1"

        else:
            gt_field[0] = "0/0"

        gt_field = ":".join(gt_field)

        i[9] = gt_field

        print("\t".join(i))

    else:
        pass
