'''
buildBAM.py
rebuilds a BAM file from the chain
Usage: $ python buildBAM.py <chainName> <readLocation> -ml/--multichainLoc <multichainLoc>
Example: python buildBAM.py chain1 
modified by CMB 07/2020
'''

import math
import pysam
import json
import argparse
import unicodedata
import pandas as pd
import os
import re
import zlib, sys, time, base64
import gzip
import numpy as np
#import PrintSequence
import subprocess
from subprocess import Popen, PIPE
#from Bio import SeqIO

#retrieve items from an input stream
def getReadData(chainName, streamName, multichainLoc, count, start, datadir):
	#get the number of items in a stream
        currStreamCommand = multichainLoc+'multichain-cli {} -datadir={} liststreamitems {} false {} {}'.format(chainName, datadir, streamName, count, start)
        with open(os.devnull, 'w') as devnull:
		items = subprocess.check_output(currStreamCommand.split(), stderr= devnull)
		items = json.loads(items, parse_int = int)
                return items
	print("Could not redirect stderr output from getReadData function")    

#get the number of items stored in a stream
def getStreamLen(chainName, streamName, multichainLoc, datadir):
	#get the number of items in a stream
        currStreamCommand = multichainLoc+'multichain-cli {} -datadir={} liststreamkeys {} *'.format(chainName, datadir, streamName)
        with open(os.devnull, 'w') as devnull:
		items = subprocess.check_output(currStreamCommand.split(), stderr= devnull)
		items = json.loads(items, parse_int = int)
		numOfReads=2*len(items)
                return(numOfReads)
        print("Could not redirect stderr output from getStremLen function")

#retrieve all keys for an input stream
def queryKey (chainName, streamName, key, multichainLoc, datadir, count):
	queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems {} {} false {}'.format(chainName, datadir, streamName, key, count)
	with open(os.devnull, 'w') as devnull:
                print(queryCommand.split())
		items = subprocess.check_output(queryCommand.split(), stderr=devnull)
		if("error" in items or "ERROR" in items):
			return -1
		matches = json.loads(items, parse_int= int)
		return matches
	print("Could not redirect stderr output from queryKey function")
	
#
def queryStream(chainName, chrQuery, readLoc, readLength, binLength, numBins, multichainLoc, datadir):
	currStream = readLoc/ binLength + 1
	print(currStream)
	if(currStream > numBins):
		print("Error: "+readLoc+" is a location outside the bounds of the data.")
		return -1
	streamName = "chr"+chrQuery+"stream"+str(currStream)
	currStreamCommand = multichainLoc+'multichain-cli {} -datadir={} liststreamitems {} false 100000000'.format(chainName, datadir, streamName)
	with open(os.devnull, 'w') as devnull:
		items = subprocess.check_output(currStreamCommand.split(), stderr= devnull)
		items = json.loads(items, parse_int = int)

		#determine if query previous stream
		if(currStream > 1 and (readLoc - readLength < (currStream-1) * binLength)):
			streamName = "chr"+chrQuery+"stream"+str(currStream-1)
			prevItems = queryKey(chainName, streamName, "FLANK=1", multichainLoc, 100000000000)
			return items+prevItems
		else:
			return items
	print("Could not redirect stderr output from queryStream function")


def cigparse(cigar):
    l1 = []
    num1 = ""
    for c1 in cigar:
        if c1 in '0123456789':
            num1 = num1 + c1
        else:
            l1.append([int(num1), c1])
            num1 = ""
    return (l1)


def findReads(df, readLoc):
	output = df[(df["POS"]<=readLoc) & (df["endPOS"]>=readLoc)]
        depth = 0
        reads = []
        A=pd.DataFrame(output["cigar"])
        for index, row in output.iterrows():
            if (len(row["cigar"]) == 1):
                reads.append(row["readID"])
            else:
                pos = row["POS"]
                endpos = row["endPOS"]
                loc = -1
                if (row["cigar"][0][1] == "S"):
                    pos = row["POS"] + row["cigar"][0][0]
                if (row["cigar"][len(row["cigar"])-1][1] == "S"):
                    endpos = row["endPOS"] - row["cigar"][len(row["cigar"])-1][0]
                for i in range(0, len(row["cigar"])):
                    if (row["cigar"][i][1] == "M"):
                        pos2 = pos + row["cigar"][i][0]
                    if (row["cigar"][i][1] == "D"):
                        loc = pos2 + row["cigar"][i][0]
                if (pos<=readLoc and endpos>=readLoc and loc != readLoc):
                    reads.append(row["readID"])
	return reads


def parseQuery(query):
	colon_loc = query.find(":")
	dash_loc = query.find("-")
	if colon_loc == -1:
		print("No colon found in query.")
		return -1
	if dash_loc == -1:
		print("No dash found in query.")
		return -1
	if dash_loc>=len(query)-1:
		print("No number after dash.")
		return -1
	if query[0:3] != "chr":
		print("No chromosome abbreviation in query.")
		return -1

	chrQuery = query[3: colon_loc]
	startLoc = int(query[colon_loc+1: dash_loc])
	endLoc = int(query[dash_loc+1: len(query)])
	return (chrQuery, startLoc, endLoc)


def ModifySequence(SOI, loc, mod, rl):
    #print(SOI)
    BB = []
    SOI_index = 0
    j = 0
    for i in xrange(0, len(loc)):
        modtype = loc[i][1]
        modloc = int(loc[i][0])
        if (modtype == "M"):
            BB.append(SOI[SOI_index:SOI_index + modloc])
            SOI_index = (SOI_index + modloc)

        if (modtype == "I"):
            BB.append(mod[j][0])
            j += 1

        if (modtype == "X"):
            BB.append(SOI[SOI_index:SOI_index + modloc])
            BB[i] = mod[j][0]
            SOI_index = SOI_index + len(mod[j][0])
            if (len(mod) > 1):
                j = j + 1

        if (modtype == "D"):  #deletion
            BB.append(' ')
            SOI_index = SOI_index + modloc

        if (modtype == "N"):  #skipped region/intron/same as deletion
            BB.append(' ')
            SOI_index = SOI_index + modloc

        if (modtype == "H"):  #hardclip/ essentially the same as a deletion
            BB.append(' ')
            SOI_index = SOI_index + modloc  #not included in SEQ

        if (modtype == "S"):  #for our purposes, a softclip essentially has the same effect as a mismatch
             BB.append(SOI[SOI_index:SOI_index + modloc])
             BB[i] = mod[j][0]
             SOI_index = SOI_index + len(mod[j][0])
             if (len(mod) > 1):
                 j = j + 1


    while True:
        try:
            BB.remove(' ')
        except ValueError:
            break
    BB = "".join(BB)
    return BB


def col1parser(col1):
    l1 = []
    num1 = ""
    for c1 in col1:
        if c1 in '0123456789':
            num1 = num1 + c1
        else:
            l1.append([int(num1), c1])
            num1 = ""
    return l1


def col2parser(col2):
    l2 = []
    num2 = ""
    for c2 in col2:
        if c2 not in ':-':
            if (c2 in 'XNACTGactgnx'):
                num2 = num2 + c2
            else:
                l2.append([num2, c2])
                num2 = ""
    return l2


def countbps(
        l1
):
    bps = 0
    for i in xrange(0, len(l1)):
        modloc = int(l1[i][0])
        bps = bps + modloc
    return bps


#Function to parse MDZ strings
def parseMDZ(string):
    bigarray = []
    num = ''
    bases = ''
    chunks = string.split(':')
    newstring = chunks[2]
    semifinalstring = re.split('([^0-9]*)', newstring)
    for i in semifinalstring:
        if not re.match(
                '\^', i
        ):  # Do not need to include deletions, because ModifySequence already took care of those
            bigarray.append(i)
    return bigarray


#Modify the sequence that came from ModifySeq based off of parsed MDZ strings
def ModifySeqII(seq, mdz):
    i = 0
    newseq = ""
    for j in range(len(mdz)):
        entry = mdz[j]
        if re.match('[0-9]', entry):
            if j + 1 != len(
                    mdz
            ):  #note: MDZ is based off of counting from base 1, while python counts from base 0. So if we're in the last bit of the mdz column, then we need to add one to the entry. Otherwise, it will cut off the last bp.
                entry = int(entry)
                newseq = newseq + seq[i:(i + entry)]
                i = (i + entry)
            else:
                entry = int(entry)
                newseq = newseq + seq[i:(i + len(seq))]
        else:
            newseq = newseq + entry
            i = i + (len(entry))
    return newseq

def getsequence(SOI, cigar, modcigar, mdz):
    #following is necessary for querying sequences from reference genome 
    l1 = col1parser(cigar)
    l2 = col2parser(modcigar)
    readlength = countbps(l1)
    final = ModifySequence(SOI, l1, l2, readlength)
    if (mdz == 0):
        seq = final.upper()
    else:
        seq = ModifySeqII(final.upper(), mdz)

    return(seq)


def main():
        
	parser = argparse.ArgumentParser()
	parser.add_argument("chainName", help = "the name of the chain storing the data")
	parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
	parser.add_argument("-ref", "--ref",help = "reference genome")
        parser.add_argument("-bam", "--bam",help = "BAM file name")
        parser.add_argument("-dr", "--dir", help= "chain location")
        parser.add_argument("-numc", "--numChromosome", type = int, help = "number of chromosomes in data", default = 24)
        args = parser.parse_args()

	#get info on the existing chain
	chainInfo = queryKey(args.chainName, "metaData", "chainInfo", args.multichainLoc, args.dir, 100000000)
	readLength = chainInfo[0]["data"]["json"]["readLength"]  
	binLength = chainInfo[0]["data"]["json"]["binLength"]
	numBins = chainInfo[0]["data"]["json"]["numBins"]
        chainHeader = queryKey(args.chainName, "metaData", "header", args.multichainLoc, args.dir, 1000000000)
        header=""
        for i in range(0,len(chainHeader)):
            header=header+chainHeader[i]["data"]["json"]["header"]
        
        res = []
        #cigar="30M20S"
        ref=args.ref
        bam=args.bam
        #n=getStreamLen(args.chainName, streamName, args.multichainLoc, args.dir)
        #k=round(float(n)/100)-1;
        #t=(float(float(n)/100)-k)*100
        
        with pysam.AlignmentFile(bam, "wb", text=header) as outf:
            for ii in range(1, args.numChromosome):
                for j in range(1, numBins+1):
                    streamName = "chr"+str(ii)+"stream"+str(j)
                    start=0
                    count=100000000
                    data = getReadData(args.chainName, streamName, args.multichainLoc, count, start, args.dir)
                    n = getStreamLen(args.chainName, streamName, args.multichainLoc, args.dir)
                    for item in data:
                        rdata=item['data'].values()
                        op=rdata[0]['OP']
                        b=op.replace('[','')
                        c=b.replace(']','')
                        d=c.replace('(','')
                        e=d.replace(')','')
                        f=e.replace(' ','')
                        g=f.replace('"','')
                        h=g.replace('\'','')
                        k=h.split(',')
                        l=[str(x) for x in k]
                        opt=''
                        mdz=0
                        a=0
                        for i in range(0, len(l)/2-1):
                            opt=opt+l[a]+':'+l[a+1]+'\t'
                            a=a+2
                        opt=opt+l[len(l)-2]+':'+l[len(l)-1]
                        qname=rdata[0]['QNAME']
                        flag=rdata[0]['FLAG']  
                        rname=rdata[0]['RNAME']
                        pos=rdata[0]['POS']
                        mapq=rdata[0]['MAPQ']
                        cigar=rdata[0]['CIGAR']
                        rnext=rdata[0]['RNEXT']
                        pnext=rdata[0]['PNEXT']
                        tlen=rdata[0]['TLEN']
                        modcigar=rdata[0]['MODCIGAR']

                        if cigar!='None': #and bool(re.search('^chr[1-9]$|^chr[1-2][0-9]$|^chrX$|^chrY$', str(rname))): #ignores None cigars and complex chr names
                            l1 = col1parser(str(cigar))
                            readlength = countbps(l1)
                            l2 = col2parser(str(modcigar))
                            stttr=str(rname)+":"+str(int(pos))+"-"+str(int(pos)+readlength)
                            SOI=pysam.faidx(ref,stttr)
                            record=re.sub('[-:>1234567890\n]','',SOI)
                            if 'chr' in record:
                                record=record.split('chr')[1]
                            seq=getsequence(record,str(cigar),str(modcigar),str(mdz))
                            a = pysam.AlignedSegment() 
                            a.query_name = qname
                            a.query_sequence = seq
                            a.flag = int(flag)
                            #print(rname)
                            if "X" in str(rname): #(rname=="X"):
                                #print(True)
                                rname=23
                            elif "Y" in str(rname):
                                rname=24 
                            
                            if "X" in str(rnext):
                                rnext=23
                            elif "Y" in str(rnext):
                                rnext=24

                            #print(rname)
                            if "chr" not in str(rname):
                                #a.reference_id = int(rname)-1
                                rr = int(rname) -1 
                                #rr = 0
                            else:
                                rr = re.sub('[chr]','',str(rname))
                            
                            #print(rr)
                            a.reference_id = int(rr)-1
                            a.reference_start = int(pos)-1
                            a.mapping_quality = int(mapq)
                            par = col1parser(cigar)
                            cig=()
                            str1=[]
                            str2=[]
                            for kp in range(0, len(par)):
                                if (par[kp][1] == "M"):
                                    par[kp][1] = 0
                                if (par[kp][1] == "I"):
                                    par[kp][1] = 1
                                if (par[kp][1] == "D"):
                                    par[kp][1] = 2
                                if (par[kp][1] == "N"):
                                    par[kp][1] = 3
                                if (par[kp][1] == "S"):
                                    par[kp][1] = 4
                                if (par[kp][1] == "H"):
                                    par[kp][1] = 5
                                if (par[kp][1] == "P"):
                                    par[kp][1] = 6
                                if (par[kp][1] == "="):
                                    par[kp][1] = 7
                                if (par[kp][1] == "X" ):
                                    par[kp][1] = 8
                                if (par[kp][1] == "B" ):
                                    par[kp][1] = 9
                                str1.append([par[kp][1],par[kp][0]])
                            cig = [tuple(kk) for kk in str1]
                            a.cigar = cig
                            a.next_reference_id = int(rnext)
                            a.next_reference_start= int(pnext)
                            a.template_length= int(tlen)
                            qual=''
                            for q in range(0, len(seq)):
                                qual=qual+'<'
                            a.query_qualities = pysam.qualitystring_to_array(qual)
                            b=0
                            for tt in range(0, len(l)/2):
                                str2.append([l[b],l[b+1]])
                                b=b+2
                            newlist=[tuple(kk) for kk in str2]
                            a.tags = [tuple(map(str,eachTuple)) for eachTuple in newlist]
                            outf.write(a)
        import os
        import psutil
        process = psutil.Process(os.getpid())
        print('\n\n Total memory in bytes:\n\n')
        print(process.memory_info().rss)                    
                

if __name__ == "__main__":
    main()
