from collections import defaultdict

def clump_finding(text,k,L,t):
	grouped = [text[i:i+k] for i in range(len(text)-k+1)]
	ct = defaultdict(int)
	result = set()
	for word in grouped[:L-k+1]:
		ct[word]+=1
		if ct[word]>=t:
			result.add(word)

	for i,word in enumerate(grouped[L-k+1:]):
		preword = grouped[i]
		ct[preword]-=1
		ct[word]+=1
		if ct[word]>=t:
			result.add(word)
	return result



if __name__ == '__main__':
	text = input()
	k,L,t = [ int(x) for x in input().split()]
	print(' '.join(clump_finding(text,k,L,t)))