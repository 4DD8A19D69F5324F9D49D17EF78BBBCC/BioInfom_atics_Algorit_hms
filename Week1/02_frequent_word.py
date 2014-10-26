# -*- coding: utf-8 -*-

from collections import Counter

def frequent_word(text,k):
	grouped = [text[i:i+k] for i in range(len(text)-k+1)]
	ct_grouped = Counter(grouped)
	ct_max = max(ct_grouped.values())
	return [key for key in ct_grouped if ct_grouped[key]==ct_max]

result = frequent_word(input(),int(input()))
print(' '.join(result))