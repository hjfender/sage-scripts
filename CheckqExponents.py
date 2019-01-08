#!/usr/bin/env sage -python

import re

#############################################################################
cartan = 'G'
n = 2
#############################################################################

filename = "functions/" + cartan + str(n) + ".py"
func_file = open(filename,"r")

line = func_file.readline()
while line[0:3] != "f =":
    line = func_file.readline()

arr = line.replace('(',' ').replace(')',' ').replace('+',' ').replace('-',' ').split()

odd = 0
even = 0
sumsmore = 0
sumsless = 0

for s in arr:
    if "sqrtq" in s:
        if s[0]!= 's':
            pos = re.search("s",s).start()
            s = s[pos:] 
        t = re.findall("\w+",s)
        qpower = int(t[1])
        if (qpower%2 == 0):
            even += 1
        else:
            odd += 1
        sumpowers = 0
        for i in range(2,len(t)):
            if "x" in t[i]:
                try:
                    sumpowers += int(t[i+1])
                except:
                    sumpowers += 1
        if (qpower <= sumpowers):
            sumsless += 1
        else:
            sumsmore += 1

print "Number of non-integer q powers: " + str(odd)
print "Number of integer q powers: " + str(even)
print "Number of qs whose power is more than half the sum of the powers of its following terms: " + str(sumsmore)
print "Number of qs whose power is less than half the sum of the powers of its following terms: " + str(sumsless)

func_file.close()