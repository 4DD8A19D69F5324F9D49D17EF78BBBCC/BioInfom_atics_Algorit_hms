# -*- coding: utf-8 -*-

import re

def pattern_count(text,pattern):
    return len(re.findall('(?='+pattern+')',text))

print(pattern_count(input(),input()))