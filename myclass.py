import numpy as np


class atom:
    def __init__(self,idx,typ,chg,xyz):
        self.idx = idx
        self.typ = typ
        self.chg = chg
        self.xyz = xyz

class cell:
    def __init__(self,nat,nty,box,tlt):
        self.nat = nat
        self.nty = nty
        self.box = box
        self.tlt = tlt


class pdb_atom:
    def __init__(self,aid,anm,loc,rnm,cid,rid,xa,ya,za,ele,cod,occ,fac,sid):
        self.aid =aid 
        self.anm =anm 
        self.loc =loc 
        self.rnm =rnm 
        self.cid =cid 
        self.rid =rid 
        self.xa = xa 
        self.ya = ya 
        self.za = za 
        self.ele =ele 
        self.cod =cod 
        self.occ =occ 
        self.fac =fac 
        self.sid =sid 

    