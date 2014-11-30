# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp

POS_DICT = {'A':0,'C':1,'G':2,'T':3}

def dist(strlong,strshort):
    arr1 = np.char.array(strlong).view('S1')
    arr2 = np.char.array(strshort).view('S1')
    n = len(arr1)
    m = len(arr2)
    return min(sum(arr1[i:i+m] != arr2) for i in range(n-m+1))

def is_d_motif(dnas,d,motif):
    return all(dist(dna,motif)<=d for dna in dnas)
    
def dist_dna_to_motif(dnas,motif):
    return sum(dist(dna,motif) for dna in dnas)
    
def all_dna_of_length(k,tmp=[]):
    if k==0:
        yield ''.join(tmp)
        return
    for c in 'ACGT':
        tmp.append(c)
        for dna in all_dna_of_length(k-1,tmp):
            yield dna
        tmp.pop()

def motif_enumeration(dnas,k,d):
    for tmp_motif in all_dna_of_length(k):
        if is_d_motif(dnas,d,tmp_motif):
            yield tmp_motif


def median_string(dnas,d):
    return min((tmp_motif for tmp_motif in all_dna_of_length(d)),\
    key = lambda motif:dist_dna_to_motif(dnas,motif))


def profile_most_probable(dna,mat):
    k = len(mat[0])
    logprob = -10000;
    ans = dna[:k]    
    
    for i in range(len(dna)-k+1):
        tmpprob=0        
        for p,c in enumerate(dna[i:i+k]):
            tmpprob+=np.log(mat[POS_DICT[c]][p]+1e-10)
            if mat[POS_DICT[c]][p]==0:
                tmpprob = -10000
                break
        if tmpprob>logprob:
            logprob=tmpprob
            ans = dna[i:i+k]
    return ans

def profile_rand_gen(dna,mat):
    k = len(mat[0])
    C = len(dna)-k+1
    logprobs = np.zeros(C)
    for i in range(C):
        tmpprob=0        
        for p,c in enumerate(dna[i:i+k]):
            tmpprob+=np.log(mat[POS_DICT[c]][p]+1e-20)
        logprobs[i]=tmpprob
    logprobs -= np.max(logprobs)
    probs = np.exp(logprobs)
    probs /= np.sum(probs)
    idx = np.random.choice(np.arange(C),p=probs)
    return dna[idx:idx+k]        
        
def motifs_score(motifs):
    motifs_arr = np.array(motifs).view('S1').reshape(len(motifs),len(motifs[0]))
    motif,counts = sp.stats.mode(motifs_arr,axis=0)
    return len(motifs)*len(motifs[0])-np.sum(counts,dtype=int)


def add_to_count(mat,motif,weight=1):
    if mat is None:
        mat = np.zeros([4,len(motif)])
    for i,c in enumerate(motif):
        mat[POS_DICT[c]][i]+= weight
    return mat

def greedy_motif_search(dnas,k):
    n = len(dnas)
    m = len(dnas[0])
    best_motif = [dna[:k] for dna in dnas]
        
    for p1 in range(m-k+1):

        motifs = [dnas[0][p1:p1+k]]
        print p1,motifs
        profile = add_to_count(np.ones([4,k]),motifs[0])
        for p2 in range(1,n):
            motifs.append(profile_most_probable(dnas[p2],profile))
            add_to_count(profile,motifs[-1])
        
        #if tmp_score==best_score and motifs<best_motif:
        #    best_motif = motifs
        if motifs_score(motifs)<motifs_score(best_motif):
            print motifs_score(motifs)
            best_motif = motifs
            
    return best_motif
        
def randomized_motif_search(dnas,k,initial=None):
    n = len(dnas)
    m = len(dnas[0])
    if initial is None:
        bestmotifs = [dna[i:i+k] for i,dna in zip(np.random.randint(low=0,high=m-k,size=n),dnas)]
    else:
        bestmotifs = initial
        
    def Profile(motifs):
        profile = np.ones([4,k])
        for motif in motifs:
            add_to_count(profile,motif)
        return profile
    
    def Motifs(profiles):
        return [profile_most_probable(dna,profiles) for dna in dnas]
    
    while True:
        profiles = Profile(bestmotifs)
        motifs = Motifs(profiles)
        if motifs_score(motifs)<motifs_score(bestmotifs):
            bestmotifs = motifs
        else:
            return motifs_score(bestmotifs),bestmotifs
            
def gibbs_sampler(dnas,k,initial=None,maxiters=1000):
    def Profile(motifs):
        profile = np.ones([4,k])
        for motif in motifs:
            add_to_count(profile,motif)
        return profile
    n = len(dnas)
    m = len(dnas[0])
    if initial is None:
        bestmotifs = [dna[i:i+k] for i,dna in zip(np.random.randint(low=0,high=m-k,size=n),dnas)]
    else:
        bestmotifs = initial
        
    nowiters =0
    nowmotifs = bestmotifs[:]
    nowprofile = Profile(nowmotifs)
    while nowiters < maxiters:
        if nowiters%10==0:        
            print nowiters
        nowiters+=1
        p = np.random.randint(0,n)
        
        add_to_count(nowprofile, nowmotifs[p], -1)
        nowmotifs[p] = profile_rand_gen(dnas[p],nowprofile)
        add_to_count(nowprofile, nowmotifs[p], 1)
        
        best_score = motifs_score(bestmotifs)
        if motifs_score(nowmotifs)<best_score:
            print 'ans=',motifs_score(nowmotifs)
            bestmotifs = nowmotifs[:]        
    return motifs_score(bestmotifs),bestmotifs

