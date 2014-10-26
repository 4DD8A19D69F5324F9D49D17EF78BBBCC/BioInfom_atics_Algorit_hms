
def reverse_complement(text):
	transform_map = {'A':'T','T':'A','C':'G','G':'C'}
	return ''.join([transform_map[c] for c in text[::-1]])

if __name__ == '__main__':
	print(reverse_complement(input()))