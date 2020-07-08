import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt


from myclass import *


def log_lmp(fins):
    all_data = []
    for fin in fins:
        with open(fin,'r') as ff1:
            n = 0
            lines = ff1.readlines()
            while "Step Elaps" not in lines[n]:
                n+=1
            head = lines[n].split()
            data = [ [] for h in head ]
            n+=1
            while "Loop time" not in lines[n]:
                line = lines[n].split()
                for i in range(len(head)):
                    data[i].append( float(line[i]) )
                n+=1
                if n==len(lines):
                    break
        all_data.append(data)
    if len(all_data)==1:
        all_data = all_data[0]
    return head,all_data

