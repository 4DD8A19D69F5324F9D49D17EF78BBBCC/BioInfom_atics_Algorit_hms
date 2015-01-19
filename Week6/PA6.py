import sys
import re
from collections import defaultdict

def parse_array_str(array_str):
    spt = array_str[1:-1].split()
    return [int(x) for x in spt]


def parse_cycle_str(cycle_str):
    return [[int(x.lstrip('(').rstrip(')')) for x in s.split()] for s in cycle_str.split(')(')]
    
def cycle_arr_to_str(cycle_arr):
    ret = []
    for lst in cycle_arr:
        def f(x):
            if x>0:
                return '+'+str(x)
            else:
                return str(x)
        tmp = [f(x) for x in lst]
        ret.append('('+' '.join(tmp)+')')
    return ''.join(ret)

def greedy_sorting(seq,out=sys.stdout):
    n = len(seq)
    steps = 0
    def seq_to_string(seq):
        tmp = ['+'+str(x) if x>0 else str(x) for x in seq]
        return '('+' '.join(tmp)+')'
    for i in range(n):
        fixed = False
        if abs(seq[i]) != i+1:
            p = 0
            for j in range(i+1,n):
                if abs(seq[j]) == i+1:
                    p = j
            for j in range(i,p+1):
                seq[j]*=-1
            seq[i:p+1] = seq[i:p+1][::-1]
            fixed = True
        if fixed:
            print >>out,seq_to_string(seq)
            steps += 1
        if seq[i]<0:
            seq[i]*=-1
            print >>out,seq_to_string(seq)
            steps += 1
    return steps
def counting_breaks(seq):
    seqp = [0] + seq + [len(seq)+1]
    ret = 0
    for x,y in zip(seqp,seqp[1:]):
        ret += ((y-x) !=1)
    return ret


def buildG(lstlst,color):
    ret = {}
    for lst in lstlst:
        lp = lst+[lst[0]]
        for x,y in zip(lp,lp[1:]):
            ret[(-x,color)]=y
            ret[(y,color)]=-x
    return ret

def mergeG(g1,g2):
    r = defaultdict(set)
    for k in g1:
        r[k] = g1[k]
    for k in g2:
        r[k] = g2[k]
    return r

def ccs(g):
    uf = {k:k for k,label in g}
    cc = len(uf)
    def find(x):
        if x==uf[x]:
            return x
        else:
            uf[x] = find(uf[x])
            return uf[x]
    for s,label in g:
        t = g[(s,label)]
        if find(s) != find(t):
            uf[find(s)] = uf[find(t)]
            cc-=1
    return cc



def perform_two_break(g,a,b,c,d):
    g[(a,0)]=b
    g[(b,0)]=a
    g[(c,0)]=d
    g[(d,0)]=c

def two_break_step(g):
    for x,label in g:
        if label ==1:
            y = g[(x,1)]
            u = g[(x,0)]
            v = g[(y,0)]
            if len(set([x,y,u,v]))==4:
                perform_two_break(g, x, y, u, v)
                return


def extract_seq(g,label):
    edge = {}
    vis = set()

    for x,_label in g:
        if _label == label:
            y = g[(x,label)]
            edge[x]= y
      
    def dfs(x,route):
        vis.add(x)
        y = edge[x]
        route.append(-y)
        if -y in vis:
            return
        else:
            vis.add(y)
            dfs(-y,route)
        
    ret = []
    for x in range(1,len(edge)/2+1):
        tmp = []
        if x not in vis:
            dfs(x,tmp)
        if tmp:
            ret.append(tmp[::-1])
    return ret
        

def shared_k_mers(k,s1,s2):
    
    
    def rc(s):
        m = {'A':'T','T':'A','G':'C','C':'G'}
        return ''.join([m[c] for c in s.strip()])[::-1]
    
    def splitby(k,s):
        ret = []
        for i in range(len(s)-k+1):
            ret.append(s[i:i+k])
        return ret
    
    def group_by_index(lst):
        ret = defaultdict(list)
        for i,x in enumerate(lst):
            ret[x].append(i)
            ret[rc(x)].append(i)
        return ret
    op = lambda x: group_by_index(splitby(k, x))
    d1 = op(s1)
    d2 = op(s2)
    ret =[]
    for k in d1:
        if k in d2:
            for x in d1[k]:
                for y in d2[k]:
                    ret.append((x,y))
    return set(ret)

print greedy_sorting(parse_array_str('(+20 +7 +10 +9 +11 +13 +18 -8 -6 -14 +2 -4 -16 +15 +1 +17 +12 -5 +3 -19)'))
print counting_breaks(parse_array_str('(+20 +8 +9 +10 +11 +12 +18 -7 -6 -14 +2 -17 -16 -15 +1 +4 +13 -5 +3 -19)'))
print len(shared_k_mers(3, 'TGGCCTGCACGGTAG', 'GGACCTACAAATGGC'))