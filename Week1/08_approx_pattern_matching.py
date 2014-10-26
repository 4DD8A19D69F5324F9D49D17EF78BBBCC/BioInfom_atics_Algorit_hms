def dist(s1,s2):
	return sum(c1!=c2 for c1,c2 in zip(s1,s2))

def group_by(text,k):
	return [text[i:i+k] for i in range(len(text)-k+1)]

def approx_matching(pat,text,d):
	grouped = group_by(text,len(pat))
	return [i for i,s in enumerate(grouped) if dist(s,pat)<=d]

if __name__ == '__main__':
	pat = input()
	text = input()
	d = int(input())
	print(' '.join(map(str,approx_matching(pat,text,d))))