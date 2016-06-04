#!python

########################################################################
# noncoding_mutation_anno.py -                                         #
#  a program to annotate the nocoding mutations with the               #
#  categories below:                                                   #
#   - intergenetic mutation (between genes)                            #
#   - UTR mutation (non-procoding genes)                               #
#   - 5'-utr mutation (mutation in 5'-UTR region)                      #
#   - 3'-utr mutation (*not covered at this time)                      #
#                                                                      # 
# Author: Tanjin Xu, xut@oregonstate.edu                               #
# Org: Ramsey Lab, http://lab.saramsey.org                             #
#      EECS, Oregon State University                                   #
#                                                                      #
# Usage: noncoding_mutation_anno.py                                    #
#            <noncoding_mutation_file>                                 #
#            <gene_structure_data_file>                                #
#                                                                      #
# Change History:                                                      #
# 01/25/2016  Initial version                      tanjin              #
########################################################################

import sys

def printUsage():
  """
  Print usage of the script
  """
  print " Usage: noncoding_mutation_anno.py <mutation_file> <gene_structure_data_file> <output_file>"
  print "        <mutation_file>:  cols(#CHROM, POS, ID, REF, ALT)"
  print "        <gene_data_file>: cols(ts_name,#chrom,strand,ts_start,ts_end,protein_coding,5'_start,5'_end)"
  print "        <output_file>: writing the annotations to the file"

def readFile(fname,mode):
  try:
    infile = open(fname, mode)
  except IndexError:
    print "file read error..."
    sys.exit()
  return infile


#----------------------------------------------------------------------------------------
# build the genen structure  map
# #chrom -> {(ts_start, ts_end) -> {(ts_name,strand, protein_coding,5'_start,5'_end)}} 
#----------------------------------------------------------------------------------------
def bldGeneMap(gene_fhdl):
  geneMap = {}
  tsMap = {}
  for line in gene_fhdl:
    gene_recd = line.strip().split('\t')
    ts_name = gene_recd[0]
    chrom   = gene_recd[1]
    strand  = gene_recd[2]
    ts_start = int(gene_recd[3])
    ts_end   = int(gene_recd[4])
    procoding = gene_recd[5]
    utr5_start = 0
    utr5_end   = 0
    utr3_start = 0
    utr3_end   = 0
    if len(gene_recd) > 6 and len(gene_recd) <= 8:
      utr5_start = int(gene_recd[6])
      utr5_end   = int(gene_recd[7])

    if len(gene_recd) > 8:
      if (len(gene_recd[6]) > 0):
        utr5_start = int(gene_recd[6])
        utr5_end = int(gene_recd[7])
      utr3_start = int(gene_recd[8])
      utr3_end   = int(gene_recd[9])
    
    tscords = (ts_start, ts_end)
    tsval = (ts_name,strand,procoding,utr5_start,utr5_end,utr3_start,utr3_end)

    if geneMap.has_key(chrom):
      if (geneMap[chrom]).has_key(tscords):
        (geneMap[chrom])[tscords].append(tsval)
      else:
        (geneMap[chrom])[tscords] = [tsval]
    else:
      geneMap[chrom] = {tscords: [tsval]}

  return geneMap


#--------------------------------------------------------------------------------------#
# Search the mutation with gene structure data                                         #
#                                                                                      #
# Note: here we only handle noncoding mutations which are categorized as:              #
# - intergenetic mutation (region between genes, i.e. promoters/enhancers/silenters)   #
# - 5'-UTR mutation                                                                    #
# - 3'-UTR mutation                                                                    #
# - intron mutation                                                                    #
# - UTR mutation (overlapping with region of non-protein coding transcript             #
#                                                                                      #
# Since one gene consists of multiple transcripts and each transcript region may       #
# overlap, and also there are non-protein coding transcripts and protein-coding        #
# transcripts for the same gene. We prioritize the search of the regions of the        #
# protein-coding transcripts first and then non-protein-coding transcripts.            #
#                                                                                      #
# Meanwhile, the mutation position could be in different function region of different  #
# transcripts. Therefore, we proritize the searching order:                            #
# - 5'-UTR region                                                                      #
# - 3'-UTR region                                                                      #
# - intron region                                                                      #
#--------------------------------------------------------------------------------------#

def mutSearch(chrom, mcord, geneMap, mutMap):
  if not geneMap.has_key(chrom):
    print "the chromsome %s could not be recognized" % (chrom)
    sys.exit()

  tsmap = geneMap[chrom]

  tscords = tsmap.keys()
  found = 0
  mutlabl = ""
  labelled = 0
  res_tsinfo = ()
  for tsrange in tscords:
    if labelled == 1:
      break

    if mcord >= tsrange[0] and mcord <= tsrange[1]:
      #found a transcript region match
      found = 1

      tsinfolst = tsmap[tsrange]
      for tsinfo in tsinfolst:
        tsname = tsinfo[0]
        strand = tsinfo[1]
        procoding = tsinfo[2]
        utr5_start = int(tsinfo[3])
        utr5_end = int(tsinfo[4])
        utr3_start = int(tsinfo[5])
        utr3_end = int(tsinfo[6])

        if procoding == "p":
          #protein coding transcript
          if utr5_start > 0 and mcord >= utr5_start and mcord <= utr5_end:
            #////////////////////////////////////////////////////////#
            # Note: as long as we found a transcript and its 5'-UTR  #
            # region includes the given mutation position, then we   #
            # annotate this mutation as 5'-utr mutation although it  #
            # is not the case in other transcripts of the same gene. #
            #////////////////////////////////////////////////////////#

            #found a 5'-utr region match
            mutlabl = "5utr"
            labelled = 1
            res_tsinfo = (tsname,tsrange[0],tsrange[1],utr5_start,utr5_end)
            break
          #search 3'-utr region if it exists
          elif utr3_start > 0 and mcord >= utr3_start and mcord <= utr3_end:
            #found a 3'-utr region match
            mutlabl = "3utr"
            labelled = 1
            res_tsinfo = (tsname,tsrange[0],tsrange[1],utr3_start,utr3_end)
            break
          else:
            # Since now no matching with the current UTR regions, we label this mutation as
            # "intron" for now. If in the next iteration, i.e., the mutation position overlapps
            # with a different range, then we search all the transcripts in that range again, 
            # and if found UTR region matching, then this label will be overwritten. 
            mutlabl = "intron"
            res_tsinfo = (tsname,tsrange[0],tsrange[1],utr5_start,utr5_end,utr3_start,utr3_end)
            continue
        else:
          #non-protein coding transcript
          if len(mutlabl) == 0:
            # Similar as intron mutation handling, never overwrite the "intron" label due to
            # the rule - the intron mutation prority is higher than general non-protein coding
            # transcript mutation. And also in next iteration, the current label could be 
            # overwritten by other categories with higher priority:5'-utr,3'-utr,intron
            mutlabl = "utr"
            res_tsinfo = (tsname,tsrange[0],tsrange[1],utr5_start,utr5_end,utr3_start,utr3_end)
    else:
      #no match for this range, then go to next range
      continue

  #no transcript match was found, then label it as intergenetic
  if found == 0:
    mutlabl = "intergenetic"

  # Memorize the annotation by map as it can be reused in query of the mutation
  # in the same position
  mutMap[(chrom,mcord)] = [mutlabl] + list(res_tsinfo)

  #return to the caller
  if len(res_tsinfo) > 0:
    return [chrom,mcord,mutlabl]+list(res_tsinfo)
  else:
    return [chrom,mcord,mutlabl]

  
#---------------------------------------------------------
# annotate the mutations
#---------------------------------------------------------
def annotate(mut_fhdl, out_fhdl, geneMap):
  #skip the header
  next(mut_fhdl)
  
  #memorize the annotated mutaions in previous step
  #(chrom, pos) -> annotation 
  mutMap = {}

  for line in mut_fhdl:
    mutation = line.strip().split('\t')
    chrom = mutation[0]
    pos = int(mutation[1])
    mutid = mutation[2]
    ref = mutation[3]
    alt = mutation[4]

    if chrom == "MT":
      continue

    if mutMap.has_key((chrom,pos)):
      #the mutaion was annotated already
      mutanno = [mutid,ref,alt]+[chrom,pos]+mutMap[(chrom,pos)]
    else:
      #search the mutation
      mutanno = mutSearch(chrom, pos, geneMap, mutMap)
      mutanno = [mutid, ref, alt] + mutanno

    #print mutanno
    out_fhdl.write('\t'.join(str(i) for i in mutanno))
    out_fhdl.write('\n')


#----------------------------------------------------------
# main function
#----------------------------------------------------------
if len(sys.argv) < 4:
  printUsage()
  sys.exit()

mutation_dat = sys.argv[1]
gene_dat = sys.argv[2]
anno_outfile = sys.argv[3]

gene_dathdl = readFile(gene_dat,'r')
mut_outhdl = readFile(anno_outfile,'w+')

genes_map = bldGeneMap(gene_dathdl)

mut_inhdl = readFile(mutation_dat,'r')

#run the annotation function
annotate(mut_inhdl, mut_outhdl, genes_map)

mut_inhdl.close()
mut_outhdl.close()
gene_dathdl.close()

#print genes_map["1"] 
