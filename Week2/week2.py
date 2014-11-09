# -*- coding: utf-8 -*-
from Bio.Seq import Seq
from Bio.Alphabet import *

mass = {pro:int(mass) for pro,mass in (line.split() for line in open("integer_mass_table.txt"))}

def split_by_frame(seq,frame_len):
    seq=seq.rstrip('\n')
    ret = [ seq[i:i+frame_len] for i in range(len(seq)-frame_len+1)]
    return ret

def print_by_line(lst):
    for item in lst:
        print(item)

def print_peptides(lstlst):
    for lst in lstlst:
        print('-'.join(map(str,lst)))

def rc(seq):
    return str(Seq(seq).reverse_complement())

def protein_translation(rna):
    return str(Seq(rna,generic_rna).translate())

def peptide_encoding(rna,pro):
    pro = pro.rstrip('\n')
    frame_len = len(pro)*3
    rna_div = split_by_frame(rna,frame_len)
    return [rna_seg for rna_seg in rna_div if protein_translation(rna_seg)==pro 
    or protein_translation(rc(rna_seg))==pro ]

def cyclospectrum(pro,dict_=lambda x:mass[x],linear=False):
    ret =[0]    
    for i in range(len(pro)):
        t = 0
        for j in range(i+1):
            t+=dict_(pro[j])
        ret.append(t)
        tmp = pro+pro[:i]
        if linear or i == len(pro)-1:
            tmp=pro
        for j in range(i+1,len(tmp)):
            t+=dict_(tmp[j])
            t-=dict_(tmp[j-(i+1)])
            ret.append(t)
    return sorted(ret)


def cyclopeptide_sequencing(spec):
    total = max(spec)
    specset = set(spec)
    candidate = {v for k,v in mass.items() if v in specset}
    ret = []
    def add_to_set(s,x):
        r=s.copy()
        for item in s:
            r.add(x+item)
        r.add(x)
        return r
    
    def dfs(remain=total,peptides=[]):
        if remain==0:
            if cyclospectrum(peptides[:],lambda x:x)==spec:
                ret.append(peptides[:])
            return
        if remain<0:
            return
        tspec= cyclospectrum(peptides,lambda x:x,linear=True)
        if all(item in specset for item in tspec):        
            for item in candidate:
                peptides.append(item)
                dfs(remain-item,peptides)
                peptides.pop()
    dfs(total)
    return ret