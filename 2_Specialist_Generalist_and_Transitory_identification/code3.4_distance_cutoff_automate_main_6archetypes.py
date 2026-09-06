#!/usr/bin/env python
# coding: utf-8
import re
import glob
import os
from pathlib import Path

def find_files(root_folder):
    base_path = Path(root_folder)
    matched_files = []
    for file_path in base_path.rglob("*.txt"):
        if "distribution" in file_path.name.lower():         
            matched_files.append(str(file_path))
    return matched_files

## Please mention your sample folder that extracted from distance distribution
folder_to_search = "../../pipeline_testing/distance_distribution/SMC08-N/"   
input_file = find_files(folder_to_search)
input_file.sort()

cut_off=[]

for file in input_file:
    with open(file) as f:
        data_mod=[]
        for line in f:
            for word in ['mu1', 'mu2', 'mu3', 'mu4', 'mu5', 'mu6', 'sigma1', 'sigma2', 'sigma3', 'sigma4', 'sigma5', 'sigma6']:
                  if word in line:
                        data_mod.append(re.findall(r'\d+\.\d+', line)[0])
        f.close()

    data=[]
    for i in data_mod:
        data.append(float(i))

    x = open(file)
    for line in x:
        if line.startswith("    Model(tri"):
            model = line[10:18]
        elif line.startswith("    Model(bi"):
            model= line[10:17]
        elif line.startswith("    Model(tetra"):
            model= line[10:20]
        elif line.startswith("    Model(penta"):
            model= line[10:20]
        elif line.startswith("    Model(hexa"):
            model= line[10:19]    
        elif line.startswith("    Model(gauss"):
            model= line[10:15]
    x.close()
    
      #for unimodal model
    if (model=="gauss"):
        mean=data[0]
        sigma= -(data[1])
    
    #for bimodal model
    if (model=="bimodal"):
        mean=min(data[0],data[2])
        if(mean==data[0]):
            sigma=data[1]
        else:
            sigma=data[3]

    #for trimodal model
    if(model=="trimodal"):
        mean=min(data[0],data[2],data[4])
        if(mean==data[0]):
            sigma=data[1]
        elif(mean==data[2]):
            sigma=data[3]
        else:
            sigma=data[5]

    #for tetramodal model
    if(model=="tetramodal"):
        mean=min(data[0],data[2],data[4],data[6])
        if(mean==data[0]):
            sigma=data[1]
        elif(mean==data[2]):
            sigma=data[3]
        elif(mean==data[4]):
            sigma=data[5]
        else:
            sigma=data[7]

    if(model=="pentamodal"):
        mean=min(data[0],data[2],data[4],data[6], data[8])
        if(mean==data[0]):
            sigma=data[1]
        elif(mean==data[2]):
            sigma=data[3]
        elif(mean==data[4]):
            sigma=data[5]
        elif(mean==data[6]):
            sigma=data[7]
        else:
            sigma=data[9]

    if(model=="hexamodal"):
        mean=min(data[0],data[2],data[4],data[6], data[8], data[10])
        if(mean==data[0]):
            sigma=data[1]
        elif(mean==data[2]):
            sigma=data[3]
        elif(mean==data[4]):
            sigma=data[5]
        elif(mean==data[6]):
            sigma=data[7]
        elif(mean==data[8]):
            sigma=data[9]
        else:
            sigma=data[10]

    cut_off.append(float(mean) + float(sigma))

print(cut_off[0])
print(cut_off[1])
print(cut_off[2])
print(cut_off[3])
print(cut_off[4])
print(cut_off[5])