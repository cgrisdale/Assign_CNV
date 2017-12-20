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
  cnvl,genes,types,samples={},[],[],[]
  tmp=[0,xrange(1),1,'gain','test','gene1','ENSG00000'] #need initial values to avoid IndexError on first round through loop

  #Go through all CNV lines; list of lists
  for x in range(len(cnvlist)):
    chrm,cnvr,cpn,etype,filen=cnvlist[x][0:]

    #Go through gtfd for current chromosome, key=xrange
    for i in gtfd[chrm]:

      if range(max(cnvr[0],i[0]),min(cnvr[-1],i[-1])+1): #range matches range key
        #ranges overlap, make note of gene matching cnv region
        genen,ensid=gtfd[chrm][i][3],gtfd[chrm][i][4]

        if cnvr!=tmp[1]: #first match for this CNV
          types.append(etype) #gain or loss
          samples.append(filen)

        #currentoverlap=range(max(cnvr[0],i[0]),min(cnvr[-1],i[-1])+1)
        #for k,v in cnvl[chrm].items():
        #  if genen in v[5] and filen==v[2]: #gene already in dict and sample name the same
        #    newoverlap=range(max(cnvr[0],k[0]),min(cnvr[-1],k[-1])+1)
        #    if currentoverlap>newoverlap: #which CNV to GTF overlap is bigger

        tmp=[chrm,cnvr,cpn,etype,filen,[genen],ensid]
        cnvl=Append_genenames(cnvl,cnvr,genen,tmp)
        #cnvl.append(tmp)
        genes.append(genen)

  print '\t'.join([k+': '+str(v) for k,v in Counter(types).items()]),'\n'
  print Counter(samples).most_common(10),'\n'
  print Counter(genes).most_common(10)

  return cnvl

def Append_genenames(cnvd,cnvrange,genename,tmpl):
  #
  if tmpl[0] in cnvd: #if chrom is in dict, check if range is in dict[chrom]
    pass

  else: #initialize dict for this chrom and add current range with value to internal dict
    cnvd[tmpl[0]]={}
    cnvd[tmpl[0]][cnvrange]=tmpl
    return cnvd

  if cnvrange in cnvd[tmpl[0]]: #if CNV range/event already in dict[chrom]
    cnvd[tmpl[0]][cnvrange][5].append(genename) #index 5 is list of gene names

  else: #if cnvrange/key not in dict
    cnvd[tmpl[0]][cnvrange]=tmpl

  return cnvd

#def Count_recurrence():
#Need to find most common overlapping regions, output in descending order of recurrence

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

  outf=open('CNV.outfile1.tsv', 'w')
  y=0
  for l,m in c.items():
    #y+=1
    #print k,v
    #if y>9:
    #  sys.exit(0)
    for k,v in m.items():
      tmpx=[v[0],v[1][0],v[1][-1],v[2],v[3],v[4],v[6]]
      genec=len(v[5])
      genenm=';'.join(str(z) for z in v[5])
      tmpx.append(genec) #append gene count and string of gene names separated by ;
      tmpx.append(genenm)
      myline='\t'.join(str(y) for y in tmpx)
      outf.write(myline)
      outf.write('\n')
  outf.close()






