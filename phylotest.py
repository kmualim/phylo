#!/usr/bin/env python3

import pandas as pd 
import numpy as np
import UPGMAalgo as up
from UPGMAalgo import UPGMA

def zeros(shape):
	D = []
	for x in range(shape[0]):
	    D.append([])
	    for y in range(shape[1]):
	        D[-1].append(0)
	return D

def affineGap(s1, s2):  # Implementing needleman algorithm with affine gap penalty
	print(s1, s2)
	match = 1
	mismatch =1
	gap_penalty = 4
	gap_extend = 1
	m = len(s1)
	n = len(s2)  # length of two sequences
	# Generate DP table and traceback path pointer matrix
	midscore = zeros((m+1, n+1))      # the DP table
	upperscore = zeros((m+1, n+1))    # initializing upperscore
	lowerscore = zeros((m+1, n+1))    # initializing lowerscore
	midscore[0][0] = 0  # initializing midscore 
	midscore[0][1] = -gap_penalty
	midscore[1][0] = -gap_penalty
	for i in range(2, len(midscore)):
	    midscore[i][0] = midscore[i-1][0] - gap_extend
	for j in range(2, len(midscore[0])):
	    midscore[0][j] = midscore[0][j-1] - gap_extend

	for i in range(m + 1):   # Calculate DP table
	    upperscore[i][0] = -999999
	    upperscore[0][1] = -gap_penalty
	for i in range(2,len(upperscore[0])): 
	    upperscore[0][i] = upperscore[0][i-1] - gap_extend

	for i in range(0, n + 1):
	    lowerscore[0][i] = -999999
	    lowerscore[1][0] = -gap_penalty
	for i in range(2, len(lowerscore)):
	    lowerscore[i][0] = lowerscore[i-1][0] - gap_extend

	backtracktable = zeros((m+1, n+1))
	for i in range(n+1):
	    backtracktable[0][i] = 'e'
	for i in range(m+1):
	    backtracktable[i][0] = 's'
	backtracktable[0][0] = '*'
	purines='A,G'
	pyrimidines='C,T'
	w=0
	for i in range(1, m + 1):   # filling affine gap alignment tables
	    for j in range(1, n + 1):
	        lowerscore[i][j] = max(lowerscore[i-1][j] - gap_extend, midscore[i-1][j] - gap_penalty)
	        upperscore[i][j] = max(upperscore[i][j-1] - gap_extend, midscore[i][j-1] - gap_penalty)
	        score1 = midscore[i-1][j-1]
	        if s1[i-1] == s2[j-1]:
	            w = match
	        if s1[i-1] != s2[j-1]:
	            w -= mismatch
	        midscore[i][j] = max(score1+w, lowerscore[i][j], upperscore[i][j])

	#filling backtrack table 
	        decision = max(score1+w, lowerscore[i][j], upperscore[i][j])
	        if decision == score1+w:
	            backtracktable[i][j] = 'd'
	        elif decision == lowerscore[i][j]:
	            backtracktable[i][j] = 's'
	        elif decision == upperscore[i][j]:
	            backtracktable[i][j] = 'e'

	# Traceback and compute the alignment 
	align1, align2 = '', ''
	i,j = m,n # start from the bottom right cell
	while i > 0 and j > 0: # end touching the top or the left edge
	    if backtracktable[i][j] == 'd':
	        align1+=s1[i-1]
	        align2+=s2[j-1]
	        i-=1
	        j-=1
	    elif backtracktable[i][j] == 's':
	        align1+=s1[i-1]
	        align2+='-'
	        i-=1
	    elif backtracktable[i][j] == 'e':
	        align1+='-'
	        align2+=s2[j-1]
	        j-=1

	# Finish tracing up to the top left cell
	while i > 0:
	    align1 += s1[i-1]
	    align2 += '-'
	    i -= 1
	while j > 0:
	    align1 += '-'
	    align2 += s2[j-1]
	    j -= 1

	#alignment.finalize(align1,align2)
	finalArray = [align1[::-1], align2[::-1]]
	return (finalArray)

#def sumofPairscoring(s1,s2): #sum of pairs scoring
	

def findScore(s1, s2):  # scoring matrix for two sequences
	match,mismatch, gap,extend = 0, 0, 0, 0
	gapMemory, subScore, counts1, counts2 = 0, 0, 0, 0
	#len1 = len(s1) - alignment.newLength(s1) #removes 'X'
	#len2 = len(s2) - alignment.newLength(s2) #removes 'X'
	#print(type(s2))
	#print(s2)
	len1 = len(s1)
	len2 = len(s2)
	for i in range(len1):
	    if (s1[i] or s2[i] == '-'):
	        if(s1[i] != '-' and s2[i] != '-'):
	            if (s1[i] == s2[i]): 
	                match += 1
	            else: 
	                mismatch += 1
	            counts1 = counts1 +1
	            counts2 = counts2 +1 
	            gapMemory=0      
	        elif (s1[i] != '-' and s2[i] == '-'):
	            if (counts2 >= 0 and counts2 < len2):
	                if(gapMemory == 1):
	                    extend = extend+1
	                else: 
	                    gap=gap+1
	                    gapMemory=1
	                counts1= counts1+1
	        elif (s1[i] == '-' and s2[i] != '-'):
	            if(counts1 >= 0 and counts1 < len1):
	                if(gapMemory == 2):
	                    extend = extend+1
	                else: 
	                    gap = gap +1
	                    gapMemory =2 
	                counts2 = counts2+1                    
	    subScore = (match*1) + (mismatch*(-1)) + (gap*(-4)) + (extend*(-1))
	    return subScore


def findTotalScore(newAligned): 
	score = 0
	for i in range(len(newAligned)): 
		for j in range(1,len(newAligned)): 
			score += findScore(newAligned[i],newAligned[j])
	return score
# align progressively 
def compareSeq(clusters):
	clusters = str(clusters)
	print(clusters)
	clusters=clusters.split(',')
	newAligned=[]
	for i in range(0,len(clusters)-1):
		s1 = clusters[i]
		print(s1)
		print(clusters[i+1])
		seq1,seq2=affineGap(s1,clusters[i+1])
		newAligned.append(seq1)
		if ((i+1)==(len(clusters))-1): 
			newAligned.append(clusters[i+1])
	return newAligned

#@staticmethod  
def similarity(arr): # used to find most similar sequence in a list 
	n = len(arr)
	finalScores=[0 for x in range(n)]
	scores = [[0 for x in range(n)] for y in range(n)]
	for i in range(n):
	    sum=0
	    for j in range(n):
	        if i == j: 
	            scores[i][j] = '-'
	        else:
	            scores[i][j] = findScore(arr[i], arr[j])
	            sum += scores[i][j]
	    finalScores[i]=sum
	return finalScores

	
 
def alignSeq(s):
	#print("alignSeq called!")
	s= s.split(',')
	newAligned = [0 for i in range(len(s))]
	tempList = ' '
	Scores =similarity(s)
	#print(Scores)
	sortedScores = sorted(Scores, reverse=True) #similarity is kept in sortedScores via index 
	#print(sortedScores)
	index = Scores.index(sortedScores[0])
	Scores[index] = -99999
	# the problem is that the sequences may repeat and hence only one sequence is called 
	# update the score and make sure that the sequences get placed inside 
	index2 = Scores.index(sortedScores[1])
	tempList = affineGap(s[index], s[index2])
	newAligned[index] = tempList[0]  
	for i in range(1,len(sortedScores)): 
	    index2 = Scores.index(sortedScores[i])
	    Scores[index2] = -9999
	    tempList = affineGap(s[index], s[index2])
	    newAligned[index2] = tempList[1]
	return newAligned        

if __name__=="__main__": 

	#df = pd.read_csv('finaltestbackup.txt', sep=' ', header=None)
	#scores=[]
	#score=0
	#for row in df[0]: 
#		row=str(row)
#		row=row.split(',')
#		score += findTotalScore(row)
#		scores.append(score)
#		score=0
#	print(scores)
	def makesmall(row):
		lists=[]
		for k in range(len(row[0])):
			#print("k",k)
			#print(len(row))
			smalllist=[]
			for j in range(len(row)):
				#print()
				smalllist.append(row[j][k])
				#print(smalllist)
			lists.append(smalllist)

			#print(lists)
		return lists
	#feature 1: dimensions of puzzle
	def countdimensions(row):   
		#row=row.split(',')
		rows=len(row)
		cols=len(row[0])
		dimensions=rows*cols
		return dimensions
	
	def checkEqual(row): 
		return row[1:]==row[:-1]

	#feature 2: gaps in puzzles - might not be a significant feature
	# gaps only registered if all sequences have gaps 
	def gaps(row): 
		count=0
		row=makesmall(row)
		#print(row)
		for i in row:  
			if('-' in i): 
				if(checkEqual(i)):
					count+=1

		return count 

	#feature 3: number of GC content
	def GCcontent(row):
		count=0 
		row=makesmall(row)
		for i in row: 
			for j in range(len(i)): 
				if (i[j] == G or i[j] == C): 
					count+=1
		return count


	csv_path='inputAndOut_uniqueid_onesolution.csv'
	df = pd.read_csv(csv_path, header=None, skiprows=1) 
	#df1= pd.read_excel(csv_path)
	# size of df is 1902 
	#data = df.iloc[0:1902, 4]
	filtered_data=df.iloc[:,1]
	# this filters data by row 
	#print(len(filtered_data))
	#s=['----TA', '---GGG', '------']
	scores1=[]
	count = 0
	score=0
	puzzlesize=[]
	#count=gaps(s)
	#print(count)
	counting=[]
	dimensions=0
	for row in filtered_data: 
		row=row.split(',')
		#print(row)
		count=gaps(row)
		counting.append(count)
		#print("Gaps", count)
		#break
		#or i in row: 
			
			#break
		dimensions=countdimensions(row)
		puzzlesize.append(dimensions)
			#for j in i: 
	print("Dimensions ",puzzlesize)
	print("Counting ",counting)			
		#puzzlesize.append(dimensions)
		
			# 152447 is the length of the file
			#while(row==i): 
				#df2=df.iloc[row,4]
				#df2 = df2.split(',')

		#determine features 



#		row=str(row)
#		row=row.split(',')
#		#print(row)
#		# print(len(row))
#		x= (UPGMA(row, findScore).tree)
#		# print(x)
#		seqs=compareSeq(x)
#			#print(seqs)
#		score+=findTotalScore(seqs)
#		count+=1
#		print("This is running:", count)
#		scores1.append(score)
#		score=0
#	print(scores1)
			# totaldifference=0
			# for j in range(len(scores1)): 
				# totaldifference+=(Outputscores[i-1]-scores1[j])
			# difference.append(totaldifference)
		# print("Done with i:", i)
	# print(difference)
				 

		#finds the score of each row 
		#row = row.split(',')
		#print(row)
		#print(type(row))
		#break

	#for row in data: 
	#	row = row.split(',')
	#	score = findTotalScore(row)
	#	scores1.append(score)
	#print(scores1)
	#df1['aligned sequence']=listseq
	#copy=df1.copy()
	#del copy['original_score']
	#del copy['best_score']
	#del copy['alignment']
	#copy.to_csv('output.csv', index=False, sep=',') 

	#aligned_data=df.iloc[0:1615,5]
	#scores=[]
	#s=0
	#for row in aligned_data: 
	#	s=findTotalScore(row)
	#	scores.append(s)
	#df['output aligned scores']=scores
	

