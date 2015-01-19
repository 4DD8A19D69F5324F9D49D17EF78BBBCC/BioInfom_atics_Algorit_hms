import re
from collections import defaultdict,deque
import sys
sys.setrecursionlimit(10**8)

def parse_array(arraystr):
    sp = list(filter(len,arraystr.split('\n')))
    return [[int(x) for x in sub.split()] for sub in sp]
def parse_dag(dagstr):
    sp = dagstr.strip().split('\n')
    pat = re.compile('->|:')
    graph_array = [tuple(pat.split(sub.strip())) for sub in sp]
    graph_dict = defaultdict(list)
    for s,t,_w in graph_array:
        w = int(_w)
        graph_dict[s.strip()].append((t.strip(),w))
    return graph_dict 

def rev_graph(graph_dict):
    ret = defaultdict(list)
    for s in graph_dict:
        for t,w in graph_dict[s]:
            ret[t].append((s,w))
    return ret

def dp_change(n,a):
    dp = [10000]*(n+1)
    dp[0]=0
    for item in a:
        value = item
        sz = 1
        while value<=n:
            ndp = dp[:]
            for i in range(item,n+1):
                ndp[i] = min(ndp[i],dp[i-value]+sz)
            dp = ndp
            value*=2
            sz*=2
    print dp[n]

def dp_manhattan(down,right):
    n = len(down)
    m = len(down[0])-1
    dp = [[0]*(m+1) for i in range(n+1)]
    for i in range(0,n+1):
        for j in range(0,m+1):
            if i:
                dp[i][j] = max(dp[i-1][j]+down[i-1][j],dp[i][j])
            if j:
                dp[i][j] = max(dp[i][j-1]+right[i][j-1],dp[i][j])
    return dp[n][m]

def longest_path(start,end,graph):
    d = defaultdict(lambda : -1000000)
    p = {}
    topo = []
    que = deque([start])
    d[start]=0
    d_in = defaultdict(int)
    for s in graph:
        for t,w in graph[s]:
            d_in[t]+=1
    while que:
        s = que.popleft()
        topo.append(s)
        for t,w in graph[s]:
            d_in[t]-=1
            if d_in[t]==0:
                que.append(t)
    print topo
    for s in topo:
        for t,w in graph[s]:
            if d[s]+w>d[t]:
                d[t]=d[s]+w
                p[t]=s
    p[start]='end','end'
    now = end
    ret = []
    while now in p:
        ret.append(now)
        now = p[now]
    
    #ret.pop()
    return ret[::-1]
        


def longpath_topo(graph,start,end):
    d = defaultdict(lambda : -1000000)
    p = {}
    topo = []
    que = deque([start])
    d[start]=0
    d_in = defaultdict(int)
    for s in graph:
        for t,w,label in graph[s]:
            d_in[t]+=1
    while que:
        s = que.popleft()
        topo.append(s)
        for t,w,label in graph[s]:
            d_in[t]-=1
            if d_in[t]==0:
                que.append(t)
    for s in topo:
        for t,w,label in graph[s]:
            if d[s]+w>d[t]:
                d[t]=d[s]+w
                p[t]=s,label
    p[start]='end','end'
    now = t
    ret = []
    while now in p:
        now,label = p[now]
        ret.append(label)
    ret.pop()
    return d[end],ret





def lcs_table(s1,s2,table,localedge=False):
    graph = defaultdict(list)
    graph[(0,0)]=[]
    n = len(s1)
    m = len(s2)
    indel = -5
    for i in range(n+1):
        for j in range(m+1):
            if i-1>=0 :
                graph[(i,j)].append(((i-1,j),indel, (s1[i-1],'-') ))
            if j-1>=0 :
                graph[(i,j)].append(((i,j-1),indel, ('-',s2[j-1]) ))
            if i-1>=0 and j-1>=0:
                graph[(i,j)].append(((i-1,j-1),table[(s1[i-1],s2[j-1])] , (s1[i-1],s2[j-1])))
    
    if localedge:
        for i in range(n+1):
            for j in range(m+1):
                if i!=n or j!=m:
                    graph[(n,m)].append(((i,j),0,('','')))
        for i in range(n+1):
            for j in range(m+1):
                if i!=0 or j!=0:
                    graph[(i,j)].append(((0,0),0,('','')))
                    
    
    print 'Building Graph Complete!'
    val,path = longpath_topo(graph,(n,m),(0,0))
    zp = zip(*path)
    r1 = ''.join(zp[0])
    r2 = ''.join(zp[1])
    
    return val,r1,r2

def lcs(s1,s2):
    graph = defaultdict(list)
    graph[(0,0)]=[]
    n = len(s1)
    m = len(s2)
    indel = -0
    for i in range(n+1):
        for j in range(m+1):
            if i-1>=0:
                graph[(i,j)].append(((i-1,j),indel, (s1[i-1],'-') ))
            if j-1>=0:
                graph[(i,j)].append(((i,j-1),indel, ('-',s2[j-1]) ))
            if i-1>=0 and j-1>=0  :
                graph[(i,j)].append(((i-1,j-1),1 if s1[i-1]==s2[j-1] else indel , (s1[i-1],s2[j-1])))                 
    val,path = longpath_topo(graph,(n,m),(0,0))
    zp = zip(*path)
    r1 = ''.join(zp[0])
    r2 = ''.join(zp[1])
    
    return val,r1,r2

def lcs_fitting(s1,s2):
    graph = defaultdict(list)
    graph[(0,0)]=[]
    n = len(s1)
    m = len(s2)
    indel = -1
    for i in range(n+1):
        for j in range(m+1):
            if i-1>=0:
                graph[(i,j)].append(((i-1,j),indel, (s1[i-1],'-') ))
            if j-1>=0:
                graph[(i,j)].append(((i,j-1),indel, ('-',s2[j-1]) ))
            if i-1>=0 and j-1>=0  :
                graph[(i,j)].append(((i-1,j-1),1 if s1[i-1]==s2[j-1] else indel , (s1[i-1],s2[j-1])))                 
    
    
    for i in range(n+1):
        graph[(i,0)].append(((0,0),0,('','')))
        graph[(n,m)].append(((i,m),0,('','')))
    
    val,path = longpath_topo(graph,(n,m),(0,0))
    zp = zip(*path)
    r1 = ''.join(zp[0])
    r2 = ''.join(zp[1])
    
    return val,r1,r2

def lcs_overlap(s1,s2):
    graph = defaultdict(list)
    graph[(0,0)]=[]
    n = len(s1)
    m = len(s2)
    indel = -2
    for i in range(n+1):
        for j in range(m+1):
            if i-1>=0:
                graph[(i,j)].append(((i-1,j),indel, (s1[i-1],'-') ))
            if j-1>=0:
                graph[(i,j)].append(((i,j-1),indel, ('-',s2[j-1]) ))
            if i-1>=0 and j-1>=0  :
                graph[(i,j)].append(((i-1,j-1),1 if s1[i-1]==s2[j-1] else indel , (s1[i-1],s2[j-1])))                 
    
    
    for i in range(n+1):
        graph[(i,0)].append(((0,0),0,('','')))
    for j in range(m+1):    
        graph[(n,m)].append(((n,j),0,('','')))
    
    val,path = longpath_topo(graph,(n,m),(0,0))
    zp = zip(*path)
    r1 = ''.join(zp[0])
    r2 = ''.join(zp[1])
    
    return val,r1,r2



def lcs_affine_gap(s1,s2,table):
    n = len(s1)
    m = len(s2)
    eps = 1
    delta = 11
    class Graph(object):
        def __iter__(self):
            for i in range(n+1):
                for j in range(m+1):
                    for k in range(3):
                        yield (i,j,k)
        def __getitem__(self,key):
            i,j,k = key
            if k==0 and i!=0:
                yield (i-1,j,0),-eps,(s1[i-1],'-')
                yield (i-1,j,1),-delta,(s1[i-1],'-')
            elif k==2 and j!=0:
                yield (i,j-1,2),-eps,('-',s2[j-1])
                yield (i,j-1,1),-delta,('-',s2[j-1])
            elif k==1:
                yield (i,j,0),0,('','')
                yield (i,j,2),0,('','')
                if i!=0 and j!=0:
                    yield (i-1,j-1,1),table[(s1[i-1],s2[j-1])],(s1[i-1],s2[j-1])
    val,path = longpath_topo(Graph(), (n,m,1), (0,0,0))
    zp = zip(*path)
    r1 = ''.join(zp[0])
    r2 = ''.join(zp[1])
    return val,r1,r2

def lcs_3d(s1,s2,s3):
    n = len(s1)
    m = len(s2)
    r = len(s3)
    ss = [s1,s2,s3]
    class Graph(object):
        def __iter__(self):
            for i in range(n+1):
                for j in range(m+1):
                    for k in range(r+1):
                        yield (i,j,k)
        def __getitem__(self,key):
            i,j,k = key
            for state in range(1,1<<3):
                pre = [i,j,k]
                gain = 0
                label = ['-']*3
                illegal = False
                for mask in range(3):
                    if state & (1<<mask):
                        pre[mask]-=1
                        if pre[mask]>=0:
                            label[mask] = ss[mask][pre[mask]] 
                        else:
                            illegal = True
                            break
                if state==7:
                    gain =1 if s1[i-1] == s2[j-1] == s3[k-1] else 0
                if illegal:
                    continue
                yield tuple(pre),gain,tuple(label)
        
    val,path = longpath_topo(Graph(), (n,m,r), (0,0,0))
    zp = zip(*path)
    return val,''.join(zp[0]),''.join(zp[1]),''.join(zp[2])

s1 = 'DACAGFEWRWVGLTDEWIMCEGQPIGSGWHHNFHYHQRPKRVFQKVRRKYYIFWRKTEKIMLSQGKSRLGLMFCMYYHPVMWLVKQFHWS'
s2 = 'DACAGFEMRMPDDWIPCEHHFTHAETAQSQYRVFRRKFWRKTNEESCFIKKIMLSQGHSRLGLMFDMYYHPVNTLVKQFHTS'
tablestr = '''
X  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7
'''.strip().split('\n')



dagstr = '''
a -> b: 5
a -> c: 6
a -> d: 5
b -> c: 2
b -> f: 4
c -> e: 4
c -> f: 3
c -> g: 5
d -> e: 6
d -> f: 8
e -> g: 2
f -> g: 1
'''


table_arr =  [line.split() for line in tablestr]

table_n = len(table_arr)
table_dict = {(table_arr[i][0],table_arr[0][j]):int(table_arr[i][j]) for j in range(1,table_n) for i in range(1,table_n)}


print parse_dag(dagstr)

print longest_path('a', 'g', parse_dag(dagstr))


