#!/bin/bash

for line in open("alleles.dbs", "r"):
    print(line.strip().split())
