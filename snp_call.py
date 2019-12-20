#!/usr/bin/python
import os, sys, math, string, operator

def factor(a,b):
	result = 1.0
	for i in range(a-b+1,a+1):
		result = result * float(i)
	for j in range(1,b+1):
		result = result/float(j)
	return result

def snp_call(nt1,nt2,count1,count2,error_rate):
	prior1 = 0.49
	prior2 = 0.02
	prior3 = 0.49
	error1 = error_rate[nt1][nt2]
	error2 = error_rate[nt2][nt1]
	if count1 + count2 < 2:
		return "-,-"
	else:
		prob1 = prior1 * factor(count1+count2,count1) * math.exp(count1*math.log(1-error1) + count2*math.log(error1))
		condition2 = 0
		for i in range(count1 + count2 + 1):
			for j in range(i+1):
				if j <= count2 and j >= i - count1:
					factor1 = factor(count1+count2,i)
					factor2 = factor(i,j)
					factor3 = factor(count1+count2-i,count2-j)
					condition2 += math.exp((count1+count2)*math.log(0.5) + j*math.log(error1) + (i-j)*math.log(1-error1) + (count1-i+j)*math.log(error2) + (count2-j)*math.log(1-error2)) * factor1 * factor2 *factor3
		prob2 = prior2 * condition2
		prob3 = prior3 * factor(count1+count2,count2) * math.exp(count2*math.log(1-error2) + count1*math.log(error2))
		newprob1 = prob1/(prob1 + prob2 + prob3)
		newprob2 = prob2/(prob1 + prob2 + prob3)
		newprob3 = prob3/(prob1 + prob2 + prob3)
		if newprob1 >= 0.8:
			return nt1 + "," + nt1
		elif newprob2 >= 0.8:
			return nt1 + "," + nt2
		elif newprob3 >= 0.8:
			return nt2 + "," + nt2
		else:
			return "-,-"

error_rate = {}
error_rate['A'] = {}
error_rate['A']['A'] = 0.99836
error_rate['A']['C'] = 0.00090
error_rate['A']['T'] = 0.00033
error_rate['A']['G'] = 0.00041
error_rate['C'] = {}
error_rate['C']['A'] = 0.00050
error_rate['C']['C'] = 0.99895
error_rate['C']['T'] = 0.00028
error_rate['C']['G'] = 0.00027
error_rate['T'] = {}
error_rate['T']['A'] = 0.00030
error_rate['T']['C'] = 0.00043
error_rate['T']['T'] = 0.99837
error_rate['T']['G'] = 0.00090
error_rate['G'] = {}
error_rate['G']['A'] = 0.00028
error_rate['G']['C'] = 0.00025
error_rate['G']['T'] = 0.00052
error_rate['G']['G'] = 0.99895

fp = open(sys.argv[1],"r")
rp = open(sys.argv[2],"w")
rp.write("scaffold\tposition\tgenotype\tcoverage\n")
for countl, line in enumerate(fp):
	if countl:
		words = line.split()
		allchars = words[2].split(",")
		if len(allchars) == 1:
			rp.write(words[0] + "\t" + words[1] + "\t-,-\t1\n")
		else:
			nrchars = []
			char_counts = {}
			for char in allchars:
				try:
					char_counts[char] += 1
				except KeyError:
					char_counts[char] = 1
					nrchars.append(char)
			if len(nrchars) == 1:
				rp.write(words[0] + "\t" + words[1] + "\t" + nrchars[0] + "," + nrchars[0] + "\t" + str(len(allchars)) + "\n")
			elif len(nrchars) == 2:
				snp = snp_call(nrchars[0], nrchars[1], char_counts[nrchars[0]], char_counts[nrchars[1]], error_rate)
				rp.write(words[0] + "\t" + words[1] + "\t" + snp + "\t" + str(len(allchars)) + "\n")
			else:
				rp.write(words[0] + "\t" + words[1] + "\t-,-\t" + str(len(allchars)) + "\n")
fp.close()
rp.close()

