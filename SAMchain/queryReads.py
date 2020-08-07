'''
queryReads.py
Query number of reads at a specific location from a chain
Usage: $ python queryReads.py <chainName> <readLocation> -ml/--multichainLoc <multichainLoc> -ref <path to ref file>
Example: python queryReads.py chain1 chr1:1005-2005 -ref GRCh38_no_alt_analysis_set_GCA_000001405.15.fa
modified by CMB on 07/2020
'''
from __future__ import print_function
import math
import json
import argparse
import pandas as pd
import subprocess
import re
import pysam
from subprocess import Popen, PIPE

def queryKey (chainName, streamName, key, datadir, multichainLoc):
	queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems {} {}'.format(chainName, datadir, streamName, key)
	items = subprocess.check_output(queryCommand.split())
	if("error" in items or "ERROR" in items):
		return -1
	matches = json.loads(items, parse_int= int)
	return matches
	

def getStreamName(readLoc, lengthBin, lengthChromosome, lengthRead):
    streamName=[]
    num_streams = int(math.ceil(float(lengthChromosome)/float(lengthBin)))
    loc = readLoc.split(":")[1]
    loc = loc.split("-")
    chrName = readLoc.split(":")[0]
    if len(loc)>1:
        for elem in range(int(loc[0]), int(loc[1])):
            stream_num = int(elem)/lengthBin + 1
            toappend = str(chrName)+"stream"+str(stream_num)
            if toappend not in streamName:
                streamName.append(toappend)
                if elem>=(1+stream_num*lengthBin) and elem<=(1+stream_num*lengthBin+lengthRead):
                    streamName.append(str(chrName)+"stream"+str(stream_num-1))
    else:
        stream_num = int(loc[0])/lengthBin + 1
        toappend = str(chrName)+"stream"+str(stream_num)
        streamName.append(toappend)
        if loc[0]>=(1+stream_num*lengthBin) and loc[0]<=(1+stream_num*lengthBin+lengthRead):
            streamName.append(str(chrName)+"stream"+str(stream_num-1))
    
    #get preceeding stream to check for FLANK=1
    first_stream_num = int(loc[0])/lengthBin + 1
    if (first_stream_num != 1):
        streamName.insert(0, str(chrName)+"stream"+str(first_stream_num-1))
    elif (first_stream_num == 1):
        streamName.insert(0, "DEADEND")
    return streamName


def queryStreams(chainName, readLoc, datadir, multichainLoc, lengthChromosome, lengthBin, streamName, numChromosome, lengthRead, reference):
    output = ''
    num_streams = int(math.ceil(float(lengthChromosome)/float(lengthBin)))
    loc = readLoc.split(":")[1].split("-")
    begin = int(loc[0])
    for i in range(0, len(streamName)):
        #first check for reads in the first stream with FLANK=1
        if (i == 0):
            if (streamName[i] != "DEADEND"):
                flankingReads = queryKey(chainName, streamName[i], "FLANK=1", datadir, multichainLoc)
                if (flankingReads != []):
                    for entry in flankingReads:
                        start = entry["data"]["json"]["POS"]
                        end = entry["data"]["json"]["endPOS"]
                        if len(loc) > 1:
                            for x in range(begin, int(loc[1])+1):
                                if x >= int(start) and x <= int(end):
                                    output = output + getItemData(entry, reference)
                                    break
                        elif len(loc)==1:
                            if int(loc[0]) >= int(start) and int(loc[0]) <= int(end):
                                output = output + getItemData(entry, reference)
                            else:
                                output = output
                elif (flankingReads == []):
                    output = output
            elif (streamName[i] == "DEADEND"):
                output = output
        #Now check for all other reads
        elif (i > 0):
            queryCommand = multichainLoc+'multichain-cli {} -datadir={} liststreamitems {} false {}'.format(chainName, datadir, streamName[i], 1000000000)
            items = subprocess.check_output(queryCommand.split())
            if("error" in items or "ERROR" in items):
                return -1
            got_items = json.loads(items, parse_int= int)
            #print(got_items[0]["data"]["json"]["POS"])
            start = ""
            end = ""
            for entry in got_items:
                start = entry["data"]["json"]["POS"]
                end = entry["data"]["json"]["endPOS"]
                if len(loc) > 1:
                    for x in range(begin, int(loc[1])+1):
                        if x >= int(start) and x <= int(end):
                            output = output + getItemData(entry, reference)
                            break
                else:
                    if (int(loc[0]) >= int(start)) and (int(loc[0]) <= int(end)):
                        output = output + getItemData(entry, reference)
    print(output)

def getItemData(streamItem, ref):
    qname=streamItem["data"]["json"]['QNAME']
    flag=streamItem["data"]["json"]['FLAG']
    rname=streamItem["data"]["json"]['RNAME']
    pos=streamItem["data"]["json"]['POS']
    endpos=streamItem["data"]["json"]['endPOS']
    mapq=streamItem["data"]["json"]['MAPQ']
    cigar=streamItem["data"]["json"]['CIGAR']
    rnext=streamItem["data"]["json"]['RNEXT']
    pnext=streamItem["data"]["json"]['PNEXT']
    tlen=streamItem["data"]["json"]['TLEN']
    modcigar=streamItem["data"]["json"]["MODCIGAR"]
    op=streamItem["data"]["json"]["OP"]
    mdz=0
    if cigar!='None' and bool(re.search('^chr[1-9]$|^chr[1-2][0-9]$|^chrX$|^chrY$', str(rname))): #ignores None cigars and complex chr names
        l1 = col1parser(str(cigar))
        readlength = countbps(l1)
        l2 = col2parser(str(modcigar))
        stttr=str(rname)+":"+str(int(pos))+"-"+str(int(pos)+readlength)
        SOI=pysam.faidx(ref, stttr)
        record=re.sub('[-:>1234567890\n]','',SOI)
        if 'chr' in record:
            record=record.split('chr')[1]
        seq=getsequence(record,str(cigar),str(modcigar),str(mdz))
        qual=''
        for q in range(0, len(seq)):
            qual=qual+'<'
        b=op.replace('[','')
        c=b.replace(']','')
        d=c.replace('(','')
        e=d.replace(')','')
        f=e.replace(' ','')
        g=f.replace('"','')
        h=g.replace('\'','')
        k=h.split(',')
        l=[str(x) for x in k]
        output = qname+"\t"+flag+"\t"+rname+"\t"+pos+"\t"+mapq+"\t"+cigar+"\t"+rnext+"\t"+pnext+"\t"+tlen+"\t"+seq+"\t"+qual+"\t" + ''
        counter=0
        for item in l:
            if counter % 2 == 0:
                output = output + item + ":" + ''
            else:
                output = output + item + "\t" + ''
            counter+=1
            if counter== len(l):
                output = output + "\n"
    return(output)


def countbps(
        l1
):  # Count basepairs from diff cigar string rather than pBAM cigar string (otherwise N's are treated as true deletions, rather than introns)
    bps = 0
    for i in xrange(0, len(l1)):
        modloc = int(l1[i][0])
        bps = bps + modloc
    return bps


def getsequence(SOI, cigar, modcigar, mdz):
    #following is necessary for querying sequences from reference genome 
    l1 = col1parser(cigar)
    l2 = col2parser(modcigar)
    readlength = countbps(l1)
    final = ModifySequence(SOI, l1, l2, readlength)
    seq = final.upper()

    return(seq)


def ModifySequence(SOI, loc, mod, rl):
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
            if (c2 in 'NACTGactgn'):
                num2 = num2 + c2
            else:
                l2.append([num2, c2])
                num2 = ""
    return l2



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("chainName", help = "the name of the chain storing the data")
	parser.add_argument("readLoc", help = "the read location to query")
	parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
	parser.add_argument("-dr", "--dir", help= "chain location")
        parser.add_argument("-lb", "--lengthBin", type = int, help = "length of each stream", default = 10000000)
        parser.add_argument("-lc", "--lengthChromosome", type = int, help = "length of the chromosome in data", default = 249000000)
        parser.add_argument("-numc", "--numChromosome", type = int, help = "number of chromosomes in data", default = 24)
        parser.add_argument("-ref", "--ref",help = "reference genome")
        args = parser.parse_args()

        #get info on the existing chain
	#get parameters of this chain from metaData stream, stored under key chainInfo.
        chainInfo = queryKey(args.chainName, 'metaData', "chainInfo", args.dir, args.multichainLoc)	
        readLength = chainInfo[0]["data"]["json"]["readLength"]
	binLength = chainInfo[0]["data"]["json"]["binLength"]
	numBins = chainInfo[0]["data"]["json"]["numBins"]
        streamName=[]
        #determine the relevant streams to check for the position queried. 
        streamName = getStreamName(args.readLoc, binLength, args.lengthChromosome, readLength)
        #The stream names obtained above have a lot of unwanted information. Loop through the data checking the POS field of each read. If the POS field falls within the range queried, get its readID. 
        data = queryStreams(args.chainName, args.readLoc, args.dir, args.multichainLoc, args.lengthChromosome, args.lengthBin, streamName, args.numChromosome, readLength, args.ref)
        #use the readID's obtained above to key into allReadData and get the information to return. 
        #data = queryAllData(args.multichainLoc, args.chainName, args.dir, key, args.ref)
	
        if(data == 0): #read was out of bounds
	    exit(1)


if __name__ == "__main__":
	main()
