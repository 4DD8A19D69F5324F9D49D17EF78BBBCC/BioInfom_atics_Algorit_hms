charset = set(['A','T','G','C'])

def gen_neighbor(text,d):
	if d==0:
		return set([text])
	else:
		result = set()
		result.add(text)
		for i,c in enumerate(text):
			for elem in charset - set([c]):
				result|=gen_neighbor(text[:i]+elem+text[i+1:],d-1)
		return result

def group_by(text,k):
	return [text[i:i+k] for i in range(len(text)-k+1)]


def freq_word_with_mismatch(text,k,d):
	grouped = group_by(text,k)
	result_dict = {}
	neighbors = [gen_neighbor(key,d) for key in grouped]

	for oneset in neighbors:
		for elem in oneset:
			if elem not in result_dict:
				result_dict[elem]=sum([elem in oneset2 for oneset2 in neighbors])

	maxvalue = max(result_dict.values())
	return [key for key in result_dict if result_dict[key]==maxvalue]


if __name__ == '__main__':
	text = input()
	k,d = [int(x) for x in input().split()]
	print(' '.join(freq_word_with_mismatch(text,k,d)))