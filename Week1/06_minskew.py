import numpy as np
def getskew(text):
	val = {'G':1,'C':-1,'A':0,'T':0}
	return np.cumsum([0]+[val[c] for c in text])

def minskew_pos(text):
	skew = getskew(text)
	minskew = min(skew)
	return [i for i,item in enumerate(skew) if item==minskew]

if __name__ == '__main__':
	text = input()
	print(' '.join(map(str,minskew_pos(text))))