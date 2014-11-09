# -*- coding: utf-8 -*-
from Bio.Seq import Seq
from collections import Counter,deque
import itertools


mass = {pro:int(mass) for pro,mass in (line.split() for line in open("integer_mass_table.txt"))}
rmass = {y:x for x,y in mass.items()}

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
        
def str_to_array(strs):
    return [int(x) for x in strs.split()]

def rc(seq):
    return str(Seq(seq).reverse_complement())

def protein_translation(rna):
    return str(Seq(rna).translate())

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
            if cyclospectrum(peptides[:],lambda x:x,linear=False)==spec:
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
    
    
def cycloscoring(seq,spec,dict_=lambda x:mass[x],linear=False):
    spec_real = cyclospectrum(seq,dict_,linear=linear)
    c1 = Counter(spec)
    c2 = Counter(spec_real)
    ret = 0
    for key in c1:
        if key in c2:
            ret+= min(c1[key],c2[key])
    return ret
    

def spectral_convlution(seq):
    return sorted(abs(a-b) for a,b in itertools.combinations(seq,2) if a!=b)


def gen_candidate(conv,m):
    ret = []
    for x,occu in Counter(conv).most_common():
        if 57<=x<=200:
            ret.append(x)
    return ret[:m]

    
def leaderboard_cyclo_sequencing(spec,n,m=20,conv=False):
    total = max(spec)
    candidate = set(mass.values())    
    if conv:    
        candidate =gen_candidate(spectral_convlution(spec),m)
    q = deque()
    toadd = []
    q.append(((),0))
    
    leader, leader_score = (),0
    iter_no=0
    while q:
        iter_no+=1
        print('Searching... %d '%(iter_no))
        while q:
            now,now_score = q.popleft()
            for item in candidate:
                tmp = now+(item,)
                tmp_score = 0
                stmp = sum(tmp)
                
                if stmp>total:
                    continue
                elif stmp==total:
                    tmp_score = cycloscoring(tmp,spec,lambda x:x,False)
                else:
                    tmp_score = cycloscoring(tmp,spec,lambda x:x,True)
                if tmp_score>leader_score:
                    leader_score = tmp_score
                    leader= tmp
                toadd.append((tmp,tmp_score))
            
        addsorted = sorted(toadd,key=lambda x:-x[1])
        if len(addsorted)>n:
            q.extend( (x,score) for x,score in addsorted if score>=addsorted[n-1][1])
        else:
            q.extend(addsorted)
            
        toadd=[]
    return leader,leader_score
        
    
    
def quiz():
    rnas = ['CCAAGUACAGAGAUUAAC','CCCAGUACCGAGAUGAAU','CCGAGGACCGAAAUCAAC','CCCAGGACUGAGAUCAAU']
    print(list(map(protein_translation,rnas)))
    peps = 'TMLA MLAT TLAM MIAT IAMT MAIT'.split()
    spec1 = str_to_array('0 71 101 113 131 184 202 214 232 285 303 315 345 416')
    print(list(filter(lambda x: cyclospectrum(x)==spec1,peps)))
    spec2 = str_to_array('0 57 71 71 71 104 131 202 202 202 256 333 333 403 404')
    print(cycloscoring('MAMA',spec2))
    spec3 = str_to_array('0 97 97 129 129 194 203 226 226 258 323 323 323 355 403 452')
    print(cycloscoring('PEEP',spec3,linear=True))
    spec4 = str_to_array('0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333')
    peps2 = 'AVQ CTV TCE TCQ VAQ ETC'.split()
    print(list(map(lambda x:cycloscoring(x,spec4,linear=True),peps2)))
    print(Counter(spectral_convlution(str_to_array('0 57 118 179 236 240 301'))).most_common(10))
    
    