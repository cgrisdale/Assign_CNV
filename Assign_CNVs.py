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
  cnvl,genes,types,samples,geneslist={},{},[],[],[]
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

        #Make dict for genes with their overlap size so they can be checked when there are duplicates etc.
        #The duplicate genes will be in different ranges! Need to check against other ranges with same genename.
        generange=range(max(cnvr[0],i[0]),min(cnvr[-1],i[-1])+1)
        tmp=[chrm,cnvr,cpn,etype,filen,genen,ensid]
        cnvl=Append_genenames(cnvl,cnvr,genen,generange,tmp)
        #cnvl.append(tmp)

        #Use dict for gene counts
        if genen in genes: #
          if filen in genes[genen]:
            pass #gene already exists with current sample name
          else:
            genes[genen].append(filen)
        else: #First time genename is added to dict, along with current sample name
          genes[genen]=[filen]

  print '\t'.join([k+': '+str(v) for k,v in Counter(types).items()]),'\n'
  print Counter(samples).most_common(10),'\n'
  #for k,v in genes.items(): #Go through dict of genes with samples as values
  #  for i in range(len(v)): #For each sample in value, add genename to list for counting
  #    geneslist.append(k)
  #print Counter(geneslist).most_common(10)

  return cnvl

def Gene_match(cnvdg,currange,genenameg,generangeg,tmplg): #generangeg=overlap
  '''Check for gene name duplicates for same sample/file name, return cnvdict'''
  templg,delcount=[],0
  for l,m in cnvdg[tmplg[0]].items(): #Loop through cnvdict[chrm] to find any matches of gene name in any range/key
    if genenameg in m[5] and tmplg[4]==m[4]: #if genename in any ranges genedict (ie. list[5]) for the same sample/filen
      if generangeg > m[5][genenameg]: #if new overlap range is larger, replace, otherwise continue
        #Assign gene to current range and delete from previous range
        del cnvdg[tmplg[0]][l][5][genenameg]
        delcount+=1
        if delcount > 1: #should not be deleting more than one duplicate
          sys.exit("Trying to delete more than one duplicate in Gene_match(), should not be possible")
        #generd=cnvdg[tmplg[0]][currange][5] #gene:range dict for current gene NOT current l,m
        #generd[genenameg]=generangeg #add current gene
        #templg=[tmplg[0],tmplg[1],tmplg[2],tmplg[3],tmplg[4],generd,tmplg[6]]
        #cnvdg[tmplg[0]][currange]=templg
        #return cnvdg
      else:
        #Previous overlap is at least as long as current, do nothing
        pass #or break because I should only find a match once?
    else:
      #Genename not found in current entry
      pass

  #If no matching genename and sample were found while looping through cnvdict
  if currange in cnvdg[tmplg[0]]: #range key present, then we're just adding a gene to innermost dict
    cnvdg[tmplg[0]][currange][5][genenameg]=generangeg

  else: #add range key for first time
    generd={}
    generd[genenameg]=generangeg
    templg=[tmplg[0],tmplg[1],tmplg[2],tmplg[3],tmplg[4],generd,tmplg[6]]
    cnvdg[tmplg[0]][currange]=templg

  return cnvdg

def Append_genenames(cnvd,cnvrange,genename,generange,tmpl):
  #
  templ=[]

  if tmpl[0] in cnvd: #if chrom is in dict, check if range is in dict[chrom]
    pass

  else: #initialize dict for this chrom and add current range with value to internal dict
    cnvd[tmpl[0]],generd={},{}

    #Create gene,range dict
    generd[genename]=generange
    templ=[tmpl[0],tmpl[1],tmpl[2],tmpl[3],tmpl[4],generd,tmpl[6]]
    cnvd[tmpl[0]][cnvrange]=templ
    return cnvd

#  if cnvrange in cnvd[tmpl[0]]: #if CNV range/event already in dict[chrom]
#    generd=cnvd[tmpl[0]][cnvrange][5]
  cnvd=Gene_match(cnvd,cnvrange,genename,generange,tmpl)
#    if genename in generd:
#      #Compare ranges
#      if generange > generd[genename]: #if new overlap range is larger, replace, otherwise continue
#        generd[genename]=generange
#    else: #generd exists, but genename not in it yet
#      generd[genename]=generange
    #After adding genename to generd, add tmpl list to cnvd
#    templ=[tmpl[0],tmpl[1],tmpl[2],tmpl[3],tmpl[4],generd,tmpl[6]]
#    cnvd[tmpl[0]][cnvrange]=templ
    #cnvd[tmpl[0]][cnvrange][5].append(genename) #index 5 is list of gene names

#  else: #if cnvrange/key not in dict, then generd doesn't exist yet for this entry
#    generd={}
#    generd[genename]=generange
#    templ=[tmpl[0],tmpl[1],tmpl[2],tmpl[3],tmpl[4],generd,tmpl[6]]
#    cnvd[tmpl[0]][cnvrange]=templ

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







