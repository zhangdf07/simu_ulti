import mdtraj as md
import numpy as np

from myclass import *


def read_lmp(fin):
    with open(fin, 'r') as ff1:
        lines = ff1.readlines()
        bbox = [ [0,0,0],[0,0,0] ]
        n = 0
        while "Atoms" not in lines[n]:
            if "atoms" in lines[n]:
                natom = int(lines[n].split()[0])
            elif "atom types" in lines[n]:
                ntyp = int(lines[n].split()[0])
            elif "xlo" in  lines[n]:
                bbox[0][0]=float(lines[n].split()[0])
                bbox[1][0]=float(lines[n].split()[1])
            elif "ylo" in  lines[n]:
                bbox[0][1]=float(lines[n].split()[0])
                bbox[1][1]=float(lines[n].split()[1])
            elif "zlo" in  lines[n]:
                bbox[0][2]=float(lines[n].split()[0])
                bbox[1][2]=float(lines[n].split()[1])
            elif "yz" in  lines[n]:
                tilt=[ float(t) for t in lines[n].split()[0:3] ]
            n+=1
        mybox = cell(natom,ntyp,bbox,tilt)
        
        n+=2
        data = []
        for i in range(natom):
            line = lines[n+i].split()
            data.append(atom(int(line[0]),int(line[1]),float(line[2]),[float(line[j]) for j in range(3,6)]))
        return mybox,data


def write_lmps(fname, mybox, data):
    with open(fname, 'w') as ff2:
        ff2.write("DF-converted \n")
        ff2.write(' \n')
        ff2.write(f" {mybox.nat} atoms " +' \n')
        ff2.write(f"{mybox.nty} atom types " + ' \n')
        ff2.write(' \n')
        ff2.write(f" {mybox.box[0][0]}  {mybox.box[1][0]}  xlo xhi " + ' \n')
        ff2.write(f" {mybox.box[0][1]}  {mybox.box[1][1]}  ylo yhi " + ' \n')
        ff2.write(f" {mybox.box[0][2]}  {mybox.box[1][2]}  zlo zhi " + ' \n')
        ff2.write(f" {mybox.tlt[0]}  {mybox.tlt[1]}  {mybox.tlt[2]}  xy xz yz " + ' \n')
        ff2.write(' \n')
        ff2.write('Atoms \n')
        ff2.write(' \n')
        for i in range(mybox.nat):
            stmp = f"  {data[i].idx}  {data[i].typ}  {data[i].chg}  {data[i].xyz[0]}  {data[i].xyz[1]}  {data[i].xyz[2]}"
            ff2.write(stmp + '\n')
    print('Finished writing to: ' + fname)



def gro2lmp(fname):
    traj = md.load(fname)
    topo = traj.topology
    bbox = [ [0,0,0], [b*10 for b in traj.unitcell_lengths[0]] ]
    natom = topo.n_atoms
    
    ntyp = 4  
    mybox = cell(natom,ntyp,bbox,[0.0,0.0,0.0])
    
    data = []
    for at in topo.atoms:
        aid = at.index + 1 
        #print(at.name)
        if 'SP' in at.name:
            aty = at.name[-1]
        elif "O" in at.name or "H" in at.name:
            if "O" in at.name:
                aty = 4
            elif "H" in at.name:
                aty = 3
            else:
                exit("Wrong water")
        else:
            
            exit("Wrong")
        chg = 0.0
        data.append(atom(aid,aty,chg,[round(j*10,3) for j in traj.xyz[0,at.index,:]]  ))
        if aid%1000==0: print(aid)
        
    return mybox,data
    
    
def write_gro(fname, elelist, cell,data):
    with open(fname, 'w') as ff2:
        ff2.write("DF-converted \n")
        ff2.write(str(cell.nat) + ' \n')
        for n in range(cell.nat):
            fmt = "{:>5}{:<5}{:>5}{:>5}{:>8}{:>8}{:>8} \n"
            stmp = fmt.format(data[n].idx,
                                'MATER',
                                elelist[data[n].typ-1], #+ str(data[n].typ),# 'SP'+str(data[n].typ), 
                                data[n].idx, 
                                round(data[n].xyz[0]*0.1,3),
                                round(data[n].xyz[1]*0.1,3),
                                round(data[n].xyz[2]*0.1,3) )
            ff2.write(stmp)
        ff2.write(f"   {round(cell.box[1][0]*0.1,4)}  {round(cell.box[1][1]*0.1,4)}  {round(cell.box[1][2]*0.1,4)}  ")
    print('Finished writing to: ' + fname)
    

def get_vdist(a,b):
    vdist = np.array(a) - np.array(b)
    return vdist
    

def shf_box(im, mybox,data):
    if im==0:
        shf = [ b for b in mybox.box[0] ]
        mybox.box[0] = [mybox.box[0][n] - shf[n] for n in range(3) ]
        mybox.box[1] = [mybox.box[1][n] - shf[n] for n in range(3) ]
        for n in range(mybox.nat):
            data[n].xyz = [ data[n].xyz[i] - shf[i] for i in range(3) ]
    elif im==1:
        print("add")
    else:
        exit("Wrong method")
    return mybox,data
    

def lmp_sort(data):
    new_data = [ 0 for n in range(len(data))]
    for d in data:
        i = d.idx-1
        new_data[i] = d
    return new_data   


def combine_lmps(b1,d1,b2,d2):
    for n in range(b2.nat):
        d2[n].idx += b1.nat
        d1.append( d2[n] )
    b1.nat += b2.nat
    return b1,d1
    

def read_xyz(fname,eledict):
    with open(fname, 'r') as ff1:
        data = []
        lines = ff1.readlines()
        for n in range(int(lines[0])):
            line = lines[n+2].split()
            data.append(atom(n+1,eledict[line[0]],0.0,[float(line[j]) for j in range(1,4)]))
    return data
    

def a2k_gro(fname):
    traj = md.load(fname)
    topo = traj.topology
    natom = topo.n_atoms
    nresd = topo.n_residues
    bbox = traj.unitcell_lengths[0]
    
    #print(set([ at.residue.name for at in topo.atoms ]))
    #exit()
    
    #a2k = [ n for n in range(natom) ]
    #a2k = [ n for n in range(natom) if topo.atom(n).residue.name not in ['POL']]
    a2k = [ n for n in range(natom) if topo.atom(n).residue.name in ['POL'] and traj.xyz[0,n,2]<4]
    #a2k = a2k[:10]
    traj.restrict_atoms(a2k)
    ##fram.save(fout2)
    return traj

