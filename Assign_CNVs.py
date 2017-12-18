#!/usr/bin/python
#Cameron Grisdale
#Dec 14, 2017

import sys
import re
from collections import Counter

#Compare output of CNV software control-freec against gene annotations

def Read_in_cnv(File):
  '''Read in CNV file, return list'''

  cnvl,cnvd,chrm,start,end,copyn,typea,tml,chrmd=[],{},'','','','','',[],{}

  with open(File, 'r') as f:

    fname=str(File)

    for line in f:
      line=line.strip()
      chrm,start,end,copyn,typea=line.split('\t')

      #Define ranges #ranges end at -1 ie. (5,8) is 5,6,7
      cnvrange=xrange(int(start),int(end))

      tml=[chrm,cnvrange,copyn,typea,fname]

      #Load biglist
      cnvl.append(tml)

      #Use range as key in internal dictionary
      #if cnvd[chrm]:
      #  cnvd[chrm][cnvrange]=tml

      #else: #initialize internal dict
      #  cnvd[chrm]={}
      #  cnvd[chrm][cnvrange]=tml

    #print len(cnvd[fname])

  return cnvl

def Read_in_GTF(gtf):
  tmp,gtfdict,gtfcount,dcount=[],{},0,0

  with open(gtf, 'r') as f:

    for line in f:
      gtfcount+=1
      chrm,start,end,strand,gn,ensid=line.strip().split()
      tmp=[chrm,xrange(int(start),int(end)),strand,gn,ensid]

      try:
        if tmp[1][0]>tmp[1][-1] or tmp[1][0]==tmp[1][-1]:
          print "Exiting: end pos is NOT higher than start pos in GTF",start,end,gn,ensid
          sys.exit(0)
      except IndexError:
        print line,tmp
        sys.exit(0)

      if chrm in gtfdict:
        gtfdict[chrm][tmp[1]]=tmp

      else:#initialize internal dict
        gtfdict[chrm]={}
        gtfdict[chrm][tmp[1]]=tmp

  for x in gtfdict:
    num=len(gtfdict[x])
    dcount+=num

  if gtfcount!=dcount:
    print "dict length doesn't match num of lines:",gtfcount,dcount
    sys.exit(0)

  return gtfdict

def Compare_ranges(cnvlist,gtfd):
  '''Go through list of CNV dict, compare to gene coordinate (dict of dict's), count frequency of recurrence'''
  cnvl,tmp,genes,types,samples=[],[],[],[],[]

  #Go through all CNV lines; list of lists
  for x in cnvlist:
    chrm,cnvr,cpn,etype,filen=x[0:]

    #Go through gtfd for current chromosome
    for i in gtfd[chrm]:
      if range(max(cnvr[0],i[0]),min(cnvr[-1],i[-1])+1): #range matches range key
        #ranges overlap, make note of gene matching cnv region
        genen,ensid=gtfd[chrm][i][3],gtfd[chrm][i][4]
        tmp=[chrm,cnvr,cpn,etype,filen,genen,ensid]
        cnvl.append(tmp)
        genes.append(genen)
        types.append(etype) #gain or loss
        samples.append(filen)

  print Counter(types),'\n'
  print Counter(samples).most_common(10),'\n'
  print Counter(genes).most_common(10)

  return cnvl


if __name__ == "__main__":

  if len(sys.argv)>2:
    GTF=sys.argv[1]
    Files=sys.argv[2:] #

  else:
    sys.exit("Script works for multiple files only; exiting")

  CV,QL=[],[]

  gtfd=Read_in_GTF(GTF)

  for i in range(len(Files)):
    cv=Read_in_cnv(Files[i])
    CV.append(cv)

  #Flatten list of lists of lists down one level to list of lists
  QL=[item for sublist in CV for item in sublist]
  print "Number of CNVs examined: ",len(QL),'\n'

  c=Compare_ranges(QL,gtfd)

  outf=open('CNV.outfile.tsv', 'w')







