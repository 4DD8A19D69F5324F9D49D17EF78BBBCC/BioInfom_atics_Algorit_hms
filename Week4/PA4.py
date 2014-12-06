# -*- coding: utf-8 -*-
from collections import defaultdict
import numpy as np
import sys

sys.setrecursionlimit(10000000)

def load_strings(fname):
    return [s.rstrip('\n') for s in open(fname).readlines()]

def load_graph(fname):
    strs = np.loadtxt(fname,dtype=str)
    graph = defaultdict(list)
    for line in strs:
        for t in line[2].split(','):
            graph[line[0]].append(t)
    return graph

def join_strings(strs):
    ret = [strs[0]]
    for c in strs[1:]:
        ret.append(c[-1])
    return ''.join(ret)


def compositions(k,text):
    return sorted(text[i:i+k] for i in range(len(text)-k+1))


def overlap_graph(strs):
    predict = defaultdict(list)
    for s in strs:
        predict[s[:-1]].append(s)
    edges = defaultdict(list)
    for s in strs:
        if s[1:] in predict:
            for s2 in predict[s[1:]]:
                edges[s].append(s2)
    return edges

def debruijn(strs):
    edges = defaultdict(list)
    for s in strs:
        edges[s[:-1]].append(s[1:])
    return edges

def euler(graph):
    d = defaultdict(int)
    ec = 0
    for x in graph:
        for y in graph[x]:
            d[y]+=1
        d[x]+=len(graph[x])
        ec+=len(graph[x])
    odds = [x for x in d if d[x]%2==1]
    def dfs(x,result,vis=None):
        if vis is None:
            vis = {}
        for i,v in enumerate(graph[x]):
            if (x,i) not in vis:
                vis[(x,i)]=1
                dfs(v,result,vis)
        result.append(x)

    if len(odds)==0:
        odds.append(list(d.keys())[0])

    if len(odds)>2:
        return '2 More Odds.'
    else:
        for s in odds:
            result = []
            dfs(s,result)
            if len(result)>=ec:
                return result[::-1]
    return 'NotFound'


def binarystrings(k):
    for i in range(1<<k):
        yield bin(i)[2:].zfill(k)

class pairedstr(object):
    def __init__(self,s1,s2):
        self.s = (s1,s2)
    def __str__(self):
        return '|'.join(self.s)
    def __repr__(self):
        return self.__str__()
    def __eq__(self,b):
        return self.s == b.s
    def __hash__(self):
        return self.s.__hash__()

    def __getitem__(self,val):
        return pairedstr(self.s[0][val],self.s[1][val])
    def to_tuple(self):
        return self.s


def decode_pairs(pairs,d):
    tups = [pair.to_tuple() for pair in pairs]
    fst,snd = list(zip(*tups))
    k = len(fst[0])
    return join_strings(fst+snd[-d-k-1:])


def maximal_non_branching_path(graph):
    paths = []
    indegrees = defaultdict(int)
    def is_one(node_no):
        return indegrees[node_no]==1 and node_no in graph and len(graph[node_no])==1
    for v in graph:
        for v2 in graph[v]:
            indegrees[v2]+=1
    
    vis = {}
    for v in graph:
        if not is_one(v):
            vis[v]=1
            for w in graph[v]:
                vis[w]=1                
                cur_path = [v,w]
                now = w
                while is_one(now):
                    vis[now]=1
                    cur_path.append(graph[now][0])
                    now = graph[now][0]
                paths.append(cur_path)
    #cycle
    for v in graph:
        if is_one(v) and not vis[v]:
            cur_path = [v]
            vis[v] = 1
            now = graph[v][0]
            while now !=v and is_one(now):
                cur_path.append(now)
                now = graph[now][0]
            paths.append(cur_path)
    
    #self cycle
    for v in graph:
        if v in graph[v]:
            cur_path.append([v,v])

    return sorted(join_strings(path) for path in paths)