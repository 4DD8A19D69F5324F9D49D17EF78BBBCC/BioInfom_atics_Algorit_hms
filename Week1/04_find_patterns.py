import re

def find_patterns(text,pattern):
	return [pat.start() for pat in re.finditer('(?='+pattern+')',text)]

if __name__ == '__main__':
	pattern = input()
	text = input()
	print(' '.join(map(str,find_patterns(text,pattern))))
