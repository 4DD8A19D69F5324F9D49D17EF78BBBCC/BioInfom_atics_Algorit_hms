def hamming_dist(s1,s2):
	return sum(c1!=c2 for c1,c2 in zip(s1,s2))
if __name__ == '__main__':
	s1 = input()
	s2 = input()
	print(hamming_dist(s1,s2))