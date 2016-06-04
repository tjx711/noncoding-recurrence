#!python

########################################################################
# noncoding_mut_stats.py -                                             #
#  a program to generate statistical info about the noncoding that     #
#  includes:                                                           #
#   - mutations recurrency over samples                                #
#   - mutations recurrency within a fixed slide window                 #
#   - mutations recurrency within specific regions                     #
#   - mutations recurrency by chance                                   #
#   - etc..                                                            #
#                                                                      # 
# Author: Tanjin Xu, xut@oregonstate.edu                               #
# Org: Ramsey Lab, http://lab.saramsey.org                             #
#      EECS, Oregon State University                                   #
#                                                                      #
# Usage: noncoding_mut_stats.py                                        #
#            <noncoding_mutation_file>                                 #
#            <window_size>                                             #
#            <out_file>                                                #
#            [with_details]                                            # 
#                                                                      #
# Change History:                                                      #
# 02/24/2016  Initial version                      tanjin              #
########################################################################

import sys
import numpy as np

def printUsage():
  """
  Print usage of the script
  """
  print " Usage: noncoding_mut_stats.py <noncoding_mutation_file> <window_size> <out_file>"
  print "        <noncoding_mutation_file>: name of the feature to be used"
  print "        <window_size>: size of sliding window"
  print "        <out_file>: file for output"
  print "        [with_details]: optional, output matching details"
#//end printUsage

def readFile(fname,mode):
  try:
    infile = open(fname, mode)
  except IndexError:
    print "file read error..."
    sys.exit()
  return infile

#----------------------------------------------------------------#
# Build the mutations map                                        #
#                                                                #
# IN: mutation file name                                         #
# OUT: map of mutations - chr->[mutation,]                       #
#----------------------------------------------------------------#
def bldMutMap(fname):
  muthdl = readFile(fname,'r')
  mutMap = {}

  for line in muthdl:
    mutation = line.strip().split('\t')
    chrom = mutation[0]

    if not mutMap.has_key(chrom):
      mutMap[chrom] = [mutation[1:]]
    else:
      mutMap[chrom].append(mutation[1:])

  muthdl.close()
  return mutMap
#//end bldMutMap

#----------------------------------------------------------------------#
# Summarize mutation recurrencies within a window of specified size    #
#                                                                      #
#  for example, given a mutation point p, and the window size 100, we  #
#  will get the recurrency by summing up the sample recurrencies of all#
#  mutations within the window [p-50,p+50]                             #
#                                                                      #
# Input:                                                               #
#   mutmap - <chrom> -> {[position,(<sampl_recur>,<mut_label>)]}       #
#   wsize - specific window size                                       #
#   outfl - file handle for output                                     #
# Output:                                                              #
#   null                                                               #
#----------------------------------------------------------------------#
def calRecurrencyWithinWindow(mutMap,wsize,outfl):
	chrlst = [str(i) for i in range(1,23,1)] + ["X","Y"]
	for ch in chrlst:
		mutlst = mutMap[ch]
		result = calRecurrencyByChrom(mutlst,wsize)

		#write result to output file
		for res in result:
			restr = ch + "\t" + str(res[0])+ "\t" + str(res[1]) + "\t" + ','.join(str(j) for j in res[2])
			#print restr
			outfl.write(restr+"\n")
		#//end for
	#//end for
#//end calRecurrencyWithinWindow

def calRecurrencyByChrom(mutlst,wsize):
	i = 0
	j = 0
	result = []

	for mut in mutlst:
		#<position>, <sampl_recur>, <mut_label>
		mcord = int(mut[0])
		w_start = mcord - (wsize / 2)
		w_end = mcord + (wsize / 2)

		first = 0
		sum_recur = 0
		matchings = []

		while i < len(mutlst):
			pcord = int(mutlst[i][0])
			if pcord >= w_start and pcord <= w_end:
				# find a mutation localized within the window
				if first == 0:
					j = i
					first = 1
				sum_recur = sum_recur + int(mutlst[i][1])
				matchings.append(pcord)
				i = i + 1
			elif pcord < w_start:
				i = i + 1
				continue
			else:
				break
		#//end while

		i = j
		result.append((mcord,sum_recur,matchings))
		if sum_recur == 0:
			print "no matching found for:",mcord
	#//end for
	return result
#//end function calRecurrencyByChrom
			

#----------------------------------------------------------
# main function                                                 
#----------------------------------------------------------
if len(sys.argv) < 4:
  printUsage()
  sys.exit()

mut_dat = sys.argv[1]
out_dat = sys.argv[2]
win_size = int(sys.argv[3])

with_details = 0
if len(sys.argv) > 4 and sys.argv[4] == "-with_details":
	with_details = 1

mutMap = bldMutMap(mut_dat)
outfl = readFile(out_dat,"w+")

calRecurrencyWithinWindow(mutMap,win_size,outfl)
