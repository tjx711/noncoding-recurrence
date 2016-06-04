#!python

########################################################################
# noncoding_extract_features.py -                                      #
#  a program to extract features for non-coding mutations              #
#  Features include:                                                   #
#   - TFBS (transcription factor binding sites)                        #
#   - ETS family TFBS                                                  #
#   - GC-percentage                                                    #
#   - DHS (DNase I Hypersensitive Sites                                #
#   - TSS (Trancription Starting Sites)                                #
#   - ...                                                              #
#                                                                      # 
# Author: Tanjin Xu, xut@oregonstate.edu                               #
# Org: Ramsey Lab, http://lab.saramsey.org                             #
#      EECS, Oregon State University                                   #
#                                                                      #
# Usage: noncoding_extract_features.py                                 #
#            <feature_name>                                            #
#            <noncoding_mutation_file>                                 #
#            <feature_bed_file>                                        #
#            <out_file>                                                #
#            [with_details]                                            # 
#                                                                      #
# Change History:                                                      #
# 02/01/2016  Initial version                      tanjin              #
########################################################################

import sys
import numpy as np

def printUsage():
  """
  Print usage of the script
  """
  print " Usage: noncoding_extract_features.py <feature_name> <mutation_file> <feature_bed_file> <out_file> [with_details]"
  print "        <feature_name>: name of the feature to be used"
  print "        <mutation_file>: cols(#CHROM, POS, ID, REF, ALT)"
  print "        <feature_bed_file>: bed format of the feature data"
  print "        <out_file>: file for output"
  print "        [with_details]: optional, output matching details"

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

#----------------------------------------------------------------#
# Build the feature map                                          #
#                                                                #
# In:                                                            #
#   fname - the feature file of bed format                       #
#   collst - list of columns to be read from the bed file except #
#            the first 3 columns which is required at any time   #
#                                                                #
# Return:                                                        #
#   *some TFs may have the same sites                            #
#   <chrom_name> -> {<start_pos,end_pos> -> {(value list)}       #
#----------------------------------------------------------------#
def bldBedFeatureMap(fname, collst):
  fbed = readFile(fname,"r")
  fmap = {}

  if len(collst) < 4:
    print "error: the columns <chr,start_pos,end_pos> are mandatory"
    sys.exit()

  for line in fbed:
    f_record = line.strip().split('\t')
    chrom      = f_record[0].replace("chr","")
    pos_start  = int(f_record[1]) + 1 #in bed format, the coordinates start with 0
    pos_end    = int(f_record[2]) + 1

    f_value = []
    for col in collst[3:]:
      f_value.append(f_record[int(col)])

    f_key = (pos_start,pos_end)

    if fmap.has_key(chrom):
      if (fmap[chrom]).has_key(f_key):
        (fmap[chrom])[f_key].append(f_value)
      else:
        (fmap[chrom])[f_key] = [f_value]
    else:
      fmap[chrom] = {f_key: [f_value]}

  fbed.close()
  return fmap
#//end bldBedFeatureMap

def bldTssFeatureMap(fname,collst):
  fens = readFile(fname,"r")
  fmap = {}

  #<chrom,ens_ts_name,tss_pos,strand>
  if len(collst) < 4:
    print "error: the columns <chr,tss_pos,ts_name,strand> must be specified"
    sys.exit()

  for line in fens:
    f_record = line.strip().split('\t')
    chrom = f_record[0].replace("chr","")
    tssite = int(f_record[2])
    tsname = f_record[1]
    strand = int(f_record[3])

    fkey = tssite
    fval = (tsname,strand)

    if fmap.has_key(chrom):
      if (fmap[chrom].has_key(fkey)):
        (fmap[chrom])[fkey].append(fval)
      else:
        (fmap[chrom])[fkey] = [fval]
    else:
      fmap[chrom] = {fkey: [fval]}

  fens.close()
  return fmap
#//end bldTssFeatureMap

'''
def bldGccMap(fname,collst):
  fgcc = readFile(fname,"r")
  fmap = {}

  chrlst = [str(i) for i in range(1,23,1)] + ["X","Y"]

  for line in fgcc:
    fstr = line.strip();
    if "chrom=" in fstr:
      # chromosome starting line
      # i.e. (variableStep chrom=chr1 span=5)
      chrom = (fstr.split(" ")[1])[9:]

      if chrom not in chrlst:
        continue
    else:
      frecord = fstr.split("\t")
      pos_start = int(frecord[0])
      gc_percent = float(frecord[1])

      if fmap.has_key(chrom):
        (fmap[chrom])[pos_start] = gc_percent
      else:
        fmap[chrom] = {pos_start:gc_percent}
  #//end for
  fgcc.close()
  return fmap
#//end bldGccMap
'''

#---------------------------------------------------------------#
# Search mutations within the feature ranges                    #
#                                                               #
# In:                                                           #
#   chrom - chromosome name                                     #
#   mutlist - list of mutations in that chromosome              #
#   fmap - feature map:(start,end)->{<meta_values>}             #
#   fops - function handling for the feature                    #
#                                                               #
# Return:                                                       #
#   features value appended to each mutation                    #
#---------------------------------------------------------------#
def mutSearch(chrom, mutlist, fmap, fops):
  f_ranges = fmap.keys()
  f_ranges.sort(key=lambda tup:tup[0]) #sort by the start coordinate
  
  pi = 0 #the starging search index for each mutation
  pj = 0 #the ending search index for the mutation

  result = []

  for mut in mutlist:
    mcord = int(mut[0])

    max_score = 0
    avg_score = 0
    total_score = 0
    total_count = 0
    total_srccount = 0
    matchings = []
    
    while pi < len(f_ranges) and mcord > f_ranges[pi][1]:
      pi = pi + 1

    if pi == len(f_ranges) or (pi < len(f_ranges) and mcord < f_ranges[pi][0]):
      # the mutation position out of the sites max range or in-between
      # two adjacent ranges
      result.append([0,0,0,0,""]) #(total_count, total_srcCount, max_score, avg_score, matchings)
      continue
    #else:
      # Found a possible match

    # memorize the current first matching position
    pj = pi

    while pj < len(f_ranges) and mcord >= f_ranges[pj][0]:
      if mcord <= f_ranges[pj][1]:
        #found a site overlapping with the mutation point
        [total_score,max_score,total_count,total_srccount] = \
           fops(f_ranges[pj], fmap, total_score,max_score,total_count,total_srccount,matchings)
      #else:
        #not a match and continue
      pj = pj + 1
    #//end while

    if total_count > 0:
      avg_score = total_score / total_count

    result.append([total_count, total_srccount, max_score, avg_score, matchings])
  #//end for
  return result
#//end mutSearch


#------------------------------------------------------------------#
# Feature function handling for each mutation                      #
#------------------------------------------------------------------#
def dhs_ops(frange, fmap, total_score, max_score, total_count, total_srccount, matchings):
  dhsvallst = fmap[frange]

  for dhsval in dhsvallst:
    total_count = total_count + 1

    if dhsval[2] != ".":
      total_srccount = total_srccount + int(dhsval[2]) #non-uniform format
    else:
      total_srccount = 0

    current_score = int(dhsval[1])
    total_score = total_score + current_score
    if current_score > max_score:
      max_score = current_score

    matchings.append((str(dhsval[0]),str(frange[0]),str(frange[1]),str(dhsval[1]),str(dhsval[2])))
  #//end for
  return [total_score, max_score, total_count, total_srccount]
#//end dhs_ops


#----------------------------------------------------------------#
# Function handling for ETS-family transcription factors         #
#----------------------------------------------------------------#
def ets_ops(frange,fmap,total_score,max_score,total_count,total_peaks,matchings):
  etsvallst = fmap[frange] #//<score,float_score,peaks>

  for etsval in etsvallst:
    total_count = total_count + 1
    if etsval[2] != ".":
      total_peaks = total_peaks + int(etsval[2])
    else:
      total_peaks = 0

    current_score = int(etsval[0])
    total_score = total_score + current_score
    if current_score > max_score:
      max_score = current_score

    matchings.append([str(i) for i in (frange[0],frange[1],etsval[0],etsval[1],etsval[2])])
  #//end for
  return [total_score,max_score,total_count,total_peaks]
#//end ets_ops

#-----------------------------------------------------------------------------#
# extract the TFBS (transcription factor binding sites) feature               #
#                                                                             #
# In:                                                                         #
#   chrom - chromosome id                                                     #
#   mutlist - list of mutaion position in that chromosome                     #
#   fmap - feature map (TF binding sites) for that chromosome                 #
# Out:                                                                        #
#   result - list of featured mutations                                       #
#-----------------------------------------------------------------------------#
def labelTFBS(chrom, mutlist, fmap):
  factorSites = fmap.keys()

  # Need to sort the keys even the input file is already sorted by key
  # Key: (start_pos, end_pos)
  factorSites.sort(key=lambda tup:tup[0]) 
  result = []

  i = 0  #the starting search index for each mutation 
  j = 0  #the ending search index for the mutation

  for mut in mutlist:
    mcord = int(mut[0])

    # record the max/average score of those overlapping sites
    maxscore   = 0
    totalscore = 0
    totalcount = 0
    totalSrcCount = 0
    avgscore   = 0
    matchings  = []
    tfMap      = {}

    while i < len(factorSites) and mcord > factorSites[i][1]:
      i = i + 1

    if i == len(factorSites) or (i < len(factorSites) and mcord < factorSites[i][0]):
      # the mutation position out of the sites max range or
      # in-between two ranges
      result.append([0,0,0,0,""])
      continue
    #else:
      # found a match
    
    # memorize the current first matching position 
    j = i
    
    while j < len(factorSites) and mcord >= factorSites[j][0]:
      if mcord <= factorSites[j][1]:
        # found a site overlapping with the mutation position
        tfsites = fmap[factorSites[j]]
        
        for site in tfsites:
          tfname = site[0] 
          tfscore = int(site[1])
          tfsrcount = int(site[2])

          if tfMap.has_key(tfname):
            tfMap[tfname].append((factorSites[j][0], factorSites[j][1], tfscore, tfsrcount))
          else:
            tfMap[tfname] = [(factorSites[j][0], factorSites[j][1], tfscore, tfsrcount)]
      #else:
        #skip this site and continue to search next matching

      j = j + 1
    #//end while

    # We only count once for each transcription factors, if the TF has multiple
    # sites range overlapping with the mutation position, we only need the range
    # with the maximum feature score.
    for tfname in tfMap.keys():
      poslist = tfMap[tfname]
      tf_maxscore = 0
      tf_maxpos = 0
      for pos in poslist:
        if pos[2] > tf_maxscore:
          tf_maxscore = pos[2]
          tf_maxpos = pos
        
      totalscore = totalscore + tf_maxscore
      totalcount = totalcount + 1
      totalSrcCount = totalSrcCount + tf_maxpos[3]

      if tf_maxscore > maxscore:
        maxscore = tf_maxscore

      matchings.append((tfname, str(tf_maxpos[0]), str(tf_maxpos[1]), str(tf_maxpos[2]), str(tf_maxpos[3])))
    #//end for tfname

    if totalcount > 0:
      avgscore = totalscore/totalcount

    #For debugging
    #print chrom,"tfbs",[str(k) for k in mut],[totalcount,totalSrcCount,maxscore,avgscore,matchings]

    result.append([totalcount,totalSrcCount,maxscore,avgscore,matchings])
  #//end for
  return result
#//end labelTFBS
 
      
#---------------------------------------------------------------#
# Extract the DNase I Hypersensitivity Sites (DHS) feature      #
#                                                               #
# IN:                                                           #
#    chrom  - chromosome name                                   #
#    mutlst - list of COSMIC mutations                          #
#    fmap   - the DHS feature for the chromosome                #
# OUT:                                                          #
#    result - list of mutations with appened feature values     #
#---------------------------------------------------------------#
def labelDNase(chrom, mutlist, fmap):
  dhsRanges = fmap.keys() #(start_pos,end_pos)
  
  # Need to sort the keys even the input file is already sorted by key
  # Key: (start_pos, end_pos)
  dhsRanges.sort(key=lambda tup:tup[0])
  result = []

  pi = 0 #the starting searching index for the previous mutation
  pj = 0 #the searching index for the current mutation

  for mut in mutlist:
    pos = int(mut[0])

    maxScore = 0
    avgScore = 0
    totalCount = 0
    totalScore = 0
    totalSrcCount = 0
    matchings = []
    
    # Skip all the ranges without matching
    while pi < len(dhsRanges) and pos > dhsRanges[pi][1]:
      pi = pi + 1 
    
    # No matchings until the last feature range or
    # the position is in-between two feature ranges
    if pi == len(dhsRanges) or (pi < len(dhsRanges) and pos < dhsRanges[pi][0]):
      # Just assign a default feature value for this position
      # and go to next mutation
      result.append([0,0,0,0,""])
      continue
    #else:
      # Found a match

    # Memorize the current first matching range index, when searching
    # the next mutation, we'll start from here
    pj = pi

    # Conintue searching the left ranges 
    while pj < len(dhsRanges) and pos >= dhsRanges[pj][0]:
      if pos <= dhsRanges[pj][1]:
        #Found a match
        dhsVal = fmap[dhsRanges[pj]]

        for dhs in dhsVal:
          totalCount = totalCount + 1
          totalSrcCount = totalSrcCount + dhs[2] 
          totalScore = totalScore + dhs[1]
          if dhs[1] > maxScore:
            maxScore = dhs[1]  

          matchings.append((str(dhs[0]),str(dhsRanges[pj][0]),str(dhsRanges[pj][1]),str(dhs[1]),str(dhs[2])))
      #else:
        #not a match and continue
      pj = pj + 1
    #//end while

    if totalCount > 0:
      avgScore = totalScore / totalCount
    
    #For debugging
    #print chrom,"dhs",[str(k) for k in mut],[totalCount,totalSrcCount,maxScore,avgScore,matchings]

    result.append([totalCount,totalSrcCount,maxScore,avgScore,matchings])
  #//end for
  return result
#//end labelDNase 

def isUniqDirTssLst(tslist):
  strands = [ts[1] for ts in tslist]
  if np.abs(np.sum(strands)) != len(strands):
    return False
  else:
    return True
#//end isValidTssLst

#-------------------------------------------------------------#
# extract the TSS (transcription starting sites) feature      #
#                                                             #
# In:                                                         #
#  chrom - chromosome name                                    #
#  mutlist - list of mutation positions in that chromosome    #
#  fmap - feature map (TS starting sites) for that chrom      #
# Out:                                                        #
#  retult - list of mutations with annotated feature          #
#-------------------------------------------------------------#
def labelTSS(chrom, mutlist, fmap):
  tssPoss = fmap.keys()

  # Need to sort the keys first
  tssPoss.sort()
  result = []

  i = 0 #starting searching index for current mutation
  j = 0 #ending searching index for current mutation

  for mut in mutlist:
    mcord = int(mut[0])
    tsdist = 0
    matchings = []

    while i < len(tssPoss) and mcord > tssPoss[i]:
      i = i + 1
    #//end while

    if i == 0 or i == len(tssPoss):
      # the mutation is NOT winthin the ranges of TSS
      if i == 0:
        j = i
      else:
        j = i - 1

      cur_tslist = fmap[tssPoss[j]]
      if isUniqDirTssLst(cur_tslist):
        cur_dist = (mcord - tssPoss[j]) * cur_tslist[0][1]
      else:
        cur_dist = - np.abs(mcord-tssPoss[j])

      matchings.append((tssPoss[j],cur_tslist))
      result.append([cur_dist,matchings])

      # Go to next mutation point
      continue
    #//end if

    # Now mcord is between tssPoss[i-1] and tssPoss[i]
    # Some exception may exist:
    #   1) different transcripts may use the same transcription
    #      starting site
    #   2) these transcripts may have different transciption direction
    #      which should not happen, if yes, print exception
    pre_tslist = fmap[tssPoss[i-1]]
    cur_tslist = fmap[tssPoss[i]]

    if isUniqDirTssLst(pre_tslist):
      pre_dist = (mcord - tssPoss[i-1]) * pre_tslist[0][1]
    else:
      # In this case, some transcripts have the same TSS coordinate but
      # located on different strands. Here we prefer the upstream - 
      # the mutation position is in upstream of the transcript, so the
      # dist(mut,tss) should be minus value
      pre_dist = -(mcord - tssPoss[i-1])

    if isUniqDirTssLst(cur_tslist):
      cur_dist = (mcord - tssPoss[i]) * cur_tslist[0][1]
    else:
      # We always prefer the upstream
      cur_dist = mcord - tssPoss[i]

    # Get the closest point
    if np.abs(pre_dist) < np.abs(cur_dist):
      tsdist = pre_dist
      matchings.append((tssPoss[i-1],pre_tslist))
    else:
      tsdist = cur_dist
      matchings.append((tssPoss[i],cur_tslist))

    # go to next mutation searching from the current starting index
    i = i - 1

    # Append to result list
    result.append([tsdist,matchings])
  #//end for

  return result
#//end function labelTSS


#-------------------------------------------------------------#
# Extract the GC-Content feature for each mutation            #
#  - calculate the GC-percentage within window of 101 bps     #
#    that is [p-50, p+50] where p is the mutation point       #
# Note: the input reference GC-content data downloaded from   #
# UCSC table browser is based on 5bps. Now we just estimate   #
# the GC-content of 101 bps window centerred at the mutation  #
# position. And also we will calculate the recurrency of the  #
# mutation using the same strategy.                           #
#                                                             #
# IN:                                                         #
#  chrom - chromosome name                                    #
#  mutlist - list of mutation positions in that chromosome    #
#  fmap - feature map (GC-Content/percentage) for that chrom  #
# OUT:                                                        #
#  result - list of mutations with annotated feature          #
#-------------------------------------------------------------#
def labelGCC(fname,mutmap,outhdl):
  chrom = ""
  fmap = {}
  resMap = {}
  fgcc = readFile(fname,'r')
  chrlst = [str(i) for i in range(1,23,1)] + ["X","Y"]

  for line in fgcc:
    fstr = line.strip();
    if "chrom=" in fstr:
      # Now we get to the next chromosome, before moving on to 
      # the next chromosome, just annotate the mutations in current
      # chromosome
      if len(chrom) > 0:
        print "start annotating chromosome:",chrom
        resMap = labelGccByChrom(chrom,mutmap[chrom],fmap,resMap)
        fmap = {}

      # chromosome starting line
      # i.e. (variableStep chrom=chr1 span=5)
      chrom = (fstr.split(" ")[1])[9:]

      if chrom not in chrlst:
        continue
    else:
      frecord = fstr.split("\t")
      pos_start = int(frecord[0])
      gc_percent = float(frecord[1])
      fmap[pos_start] = gc_percent
  #//end for

  # Need to handle the last chromosome
  print "start annotating last chromosome:",chrom
  resMap = labelGccByChrom(chrom,mutmap[chrom],fmap,resMap)
  fmap = {}

  # Write the result to output file
  for ch in chrlst:
    reschrMap = resMap[ch] #<chrom> --> <mut_coordinate> : {(gcc_score,matchings)}
    mutlst = reschrMap.keys()
    mutlst.sort()

    for mut in mutlst:
      gccres = reschrMap[mut]
      restr = ch + "\t" + str(mut)+"\t"+ str(gccres[0])
      restr = restr + "\t" + ';'.join(str(j) for j in gccres[1])
      print restr
      outhdl.write(restr+'\n')
    #//end for
  #//end for
#//end function labelGCC

# Annotate GCC feature by chromosome
def labelGccByChrom(chrom, mutlist, fmap, resMap):
  # read the feature map
  # <chrom_name> -> <start_pos> -> <percentage> 
  # *the percentage is based on 5bps-window, i.e. 
  # GC-Percent(<start_pos>, <start_pos + 4>)
  gccPoss = fmap.keys()
  feaLen = len(gccPoss)

  # Sort the keys first
  gccPoss.sort()
  result = []
  #print gccPoss

  i = 0 #staring searching index for current mutation
  j = 0 #ending searching index for current mutation

  # Traverse each mutation point
  for mut in mutlist:
    gcc = 0
    mcord = int(mut[0])
    w_start = mcord - 50
    w_end = mcord + 50
    matchings = []
    win_size = 101 # 101bps-sliding window

    while i < feaLen and w_start > (gccPoss[i] + 4):
      i = i + 1
    #//end while

    if i == 0 or i == feaLen:
      # Can not find any GC-Content 5bps window 
      # then how to handle this case? Output an exception for now
      #print "exception occured for the mutation:",i,chrom,mut
      #sys.exit()
      matchings.append((0,0,0))
      win_size = 0
      gcc = 0
    else:
      gcc = (gccPoss[i] + 4 - w_start + 1) * fmap[gccPoss[i]]
      j = i + 1

      while j < feaLen and w_end >= (gccPoss[j]+4):
        # calculate the GC-Content count
        gcc = gcc + 5 * fmap[gccPoss[j]]
        j = j + 1
      #//end while

      if j < feaLen:
        if gccPoss[j] > w_end:
          # Shrink the window
          win_size = win_size - (w_end - gccPoss[j-1] - 4)
          matchings.append((gccPoss[i],gccPoss[j-1],win_size))
        else:
          gcc = gcc + (w_end - gccPoss[j] + 1) * fmap[gccPoss[j]]
          matchings.append((gccPoss[i],gccPoss[j],win_size))
      else:
        # Here means the input GC-content data is not complete
        # Or the window [p-50,p+50] is out of range of chromosome length, for this case
        # just calculate the GC-Percent within the max window
        win_size = win_size - (w_end - gccPoss[j-1] - 4)
        matchings.append((gccPoss[i],gccPoss[j-1],win_size))

      gcc = float(gcc) / win_size
    #//end else

    # Build the result map
    if resMap.has_key(chrom):
      (resMap[chrom])[mcord] = (gcc,matchings)
    else:
      resMap[chrom] = {mcord: (gcc,matchings)}
  #//end for
  return resMap
#//end function labelGccByChrom


#-------------------------------------------------------------#
# for each mutation, append features which include:           #
#  - transcription factor binding sites                       #
#  - DNase I Hypersensitivity sites                           #
#  - GC content                                               #
#  - TSS sites                                                #
#                                                             #
# IN:                                                         #
#   mutfl - mutation file handle                              #
#   infl  - feature                                           #
#-------------------------------------------------------------#
def extFeatures(mutfl, infl, outfl, fname, outdetail):
  feaMap = {}
  collst = []
  outhdl = readFile(outfl,'w')

  #specify which columns to be read from the feature input bed file
  if fname == "-tfbs" or fname == "-dhs":
    collst = range(6)
  elif fname == "-dhs-master":
    #<chr,start,end,name,src_count,score>
    collst = [0,1,2,3,4,6]
  elif fname == "-ets":
    #<chr,start,end,score,float_score,#_of_peaks>
    collst = [0,1,2,4,6,9]
  elif fname == "-tss":
    #<chr,ts_name,tss_pos,strand>
    collst = [0,2,1,3]
  elif fname == "-gcc":
    #<chr,cord_start,gc_percent>
    collst = [0,1]
  else:
    print "Error: unkown feature name %s" % (fname)

  #build the features map
  if fname == "-tss":
    feaMap = bldTssFeatureMap(infl, collst)
  elif fname == "-gcc":
    feaMap = {}
  else:
    feaMap = bldBedFeatureMap(infl, collst)

  chrlst = [str(i) for i in range(1,23,1)] + ["X","Y"]
  mutMap = bldMutMap(mutfl)

  # Since gcc input data is too large to fit into the memory, using a different
  # procedure for gcc annotation. Specifically we cannot build the feature map
  # of all chromosomes (>6GB), instead just build the feature map by chromosome
  if fname == "-gcc":
    labelGCC(infl,mutMap,outhdl)
  # For all other features where same procedure can be used
  else:
    # Traverse mutations in each chromosome
    for ch in chrlst:
      mutlst = mutMap[ch]
      chrfea_map = {}

      # Need to sort the keys even the input file is already sorted by key
      mutlst.sort(key=lambda tup: int(tup[0]))    

      # Some features may not have any info on some chromosomes
      if not feaMap.has_key(ch):
        # Nothing to do for all mutations on this chromosome
        # but we still need to set the feature values to 0s
        # or no need to do anything and set to 0s when joining
        # the feature values with the mutation points which will
        # be done separately in a shell script
        continue

      chrfea_map = feaMap[ch]
      res = []

      if fname == "-tfbs":
        #annotate the TFBS feature
        print "---start annotating TFBS feature-----"
        res = labelTFBS(ch,mutlst,chrfea_map)
      elif fname == "-dhs" or fname == "-dhs-master":
        #annotate the DHS feature
        print "---start annotating DHS feature------"
        #dhsres = labelDNase(ch,mutlst,chrdhs_map)
        res = mutSearch(ch,mutlst,chrfea_map,dhs_ops)
      elif fname == "-ets":
        #annotate the ETS TFBS feature
        print "---start annotating ETS feature-----"
        res = mutSearch(ch,mutlst,chrfea_map,ets_ops)
      elif fname == "-tss":
        #annotate the TSS feature
        res = labelTSS(ch,mutlst,chrfea_map)
        print "---start annotating TSS feature"
      #else:
        #other unexpected cases
      
      if len(res) != len(mutlst):
        print "Annotation errors, exit now.."
        sys.exit()

      #format the output string
      i = 0
      for mut in mutlst:
        restr = ch + "\t" + '\t'.join(mut[0:1])+"\t"+'\t'.join(str(j) for j in res[i][0:-1])

        if outdetail == 1:
          restr = restr + "\t" + ';'.join(str(j) for j in res[i][-1])
      
        if outdetail == 0:
          print (restr+"\t"+';'.join(str(j) for j in res[i][-1]))
        else:
          print restr 

        outhdl.write(restr+'\n')
        i = i + 1
      #//end for
    #//end for 
  #//end else
  outhdl.close()  
#//end extFeatures
     
#----------------------------------------------------------
# main function                                                 
#----------------------------------------------------------
if len(sys.argv) < 5:
  printUsage()
  sys.exit()

fea_name = sys.argv[1]
mut_dat = sys.argv[2]
fea_dat = sys.argv[3]
out_dat = sys.argv[4]

with_details = 0

if len(sys.argv) > 5 and sys.argv[5] == "-with_details":
  with_details = 1

extFeatures(mut_dat,fea_dat,out_dat,fea_name, with_details)

#end        