'''
insertData.py
Loads tab-separated text file data (SAM or BAM format) onto an existing data stream on an existing multichain
Usage: $ python insertData.py <chainName> <length_chromosome> <length_bin> <length_read> <tab-separated-data>.bam
Example: python insertData.py -cn chain1 -lc 4050 -lb 1000 -lr 100 /home/slw67/fake_data/length_startloc.txt
Example: python insertData.py /home/slw67/nextgen-blockchain/fake_w_chr.txt
Example: python insertData.py splitbam_aa.bam splitbam_ab.bam splitbam_ac.bam splitbam_ad.bam -cn chain1
Example: python insertData.py splitbam_* -rf 10 10 10 10
modified by CMB 07/2020
modified by EN 01/2022
'''

import pysam
import sys
import time
import math
import binascii
import argparse
import subprocess
from subprocess import Popen, PIPE
import os
import psutil
import time
import multiprocessing
import glob
#define global variables
chrType = -1


#Given a chain name and the name of the new stream, subscribe to that stream
def subscribeToStream(chainName, streamName, multichainLoc, datadir):
    subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe {}'.format(chainName,datadir,streamName)

    #subscribe to the stream
    procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    procSubscribe.wait()

#keys is a tuple containing keys to insert
#data is in a form to be inserted into multichain: hex or json
def publishToStream(chainName, datadir, streamName, keys, multichainLoc, data= format(" ".encode('hex'))):
    #make string array of keys
    hasError = 0
    if len(keys) < 1:
        return -1
    #concatenate keys 
    keyString = '["'+keys[0]+'"'
    for key in keys[1:]:
        keyString += ',"'+key+'"'
    keyString += ']'
    #print(keyString)

    publishCommand = [multichainLoc+'multichain-cli', 
    str('{}'.format(chainName)), 
    str('-datadir={}'.format(datadir)),
    'publish',
    str('{}'.format(streamName)), 
    str('{}'.format(keyString)),
    str('{}'.format(data))]
    #supress output of the command
    #print(publishCommand)
    #logf=open("Insert.log", 'w+')
    dummy = subprocess.check_output(publishCommand, stderr=subprocess.STDOUT)

    return hasError

#Making the streams that will be used to store the file
#making the stream bins
#hash for location
def subscribeToStreams(chainName, numChromosome, lengthChromosome, lengthBin, lengthRead, multichainLoc, datadir):
    #data stream by readName
    subscribeToStream(chainName, "metaData", multichainLoc, datadir) #header, binLength, other things for the whole file
    subscribeToStream(chainName, "unmappedANDcontigs", multichainLoc, datadir)
    #hash for location
    num_streams = int(math.ceil(float(lengthChromosome)/float(lengthBin)))
    #range not inclusive, so run num_streams times
    for c in range(1, numChromosome+1):
        for i in range(1, num_streams+1):
            subscribeToStream(chainName, "chr"+str(c)+"stream"+str(i), multichainLoc, datadir)

    return num_streams



def determineChr(input):
    input = str(input)
    chrLookup = {
    "1": 1,"2": 2,"3": 3,"4": 4,"5": 5,
    "6": 6,"7": 7,"8": 8,"9": 9,"10": 10,
    "11": 11,"12": 12,"13": 13,"14": 14,"15": 15,
    "16": 16,"17": 17,"18": 18,"19": 19,"20": 20,
    "21": 21,"22": 22,"X": 23,"Y": 24,
    "chr1": 1,"chr2": 2,"chr3": 3,"chr4": 4,"chr5": 5,
    "chr6": 6,"chr7": 7,"chr8": 8,"chr9": 9,"chr10": 10,
    "chr11": 11,"chr12": 12,"chr13": 13,"chr14": 14,"chr15": 15,
    "chr16": 16,"chr17": 17,"chr18": 18,"chr19": 19,"chr20": 20,
    "chr21": 21,"chr22": 22,"chrX": 23,"chrY": 24
    }
    return chrLookup.get(input, -1)


def isReadWanted (flag):
    wanted = True
    allFlags = list("{0:b}".format(flag))
    #check for FLAG values 4 and 8, which indicate read or read pair is unmapped 
    fourFlag = allFlags[-3]
    eightFlag = allFlags[-4]
    if (fourFlag=='1') or (eightFlag=='1'):
        wanted = False
    else:
        wanted = True
    return(wanted)

def insertHeader(chainName, multichainLoc, data, datadir):
    header = pysam.view("-H", data)
    hed = str(header)
    ss=hed.splitlines(True)

    for i in range(0, len(ss)):
        headerData = '{"json":{"header":"'+str(ss[i])+'"}}'
        publishToStream(chainName, datadir, "metaData", ("header",), multichainLoc, headerData)
######################################################################################################

#for creating modcigar
def cigparse(cigar):
    l1 = []
    num1 = ""
    for c1 in cigar:
        if c1 in '0123456789':
            num1 = num1 + c1
        else:
            l1.append([int(num1), c1])
            num1 = ""
    return (
        l1)


def getfseq(cigar, seq):
    a = cigparse(cigar)
    #b is the sequence
    b = seq
    start = 0
    m = ''
    k = ''
    for i in range(0, len(a)):
        if a[i][1] in 'M':
            start = start + a[i][0]
            k = ''
        if a[i][1] in 'N':
            start = start
            k = ''
        if a[i][1] in 'SIX':
            tup = b[start:start + a[i][0]]
            k = tup + ':' + str(a[i][1]) + '-'
            start = start + a[i][0]
        m = m + k


    return (m[0:len(m) - 1])

#binning and hashing happens here
#we currently bin by location
def classifyReads(chainName, num_streams, lengthBin, lengthRead, multichainLoc, data, datadir, restart_from = 0): 
    hasError = 0

    num_reads = int(pysam.view('-c', data))

    with pysam.AlignmentFile(data, "rb") as bamfile:
        currRead = 0
        for line in bamfile:
            if psutil.disk_usage('/').free < 2**30:
                print('ABORTING: DISK SPACE LOW')
                print('Please restart building chain for {} from read #{}'.format(data, currRead))
                return 

            if currRead < restart_from:
                currRead +=1
                continue

            readName = line.query_name
            flag = line.flag
            chrom = line.reference_name
            startLoc = line.reference_start + 1
            seq = line.query_sequence
            cigar = line.cigarstring
            modcigar=""
            if str(cigar) == "*" or str(cigar) == "" or str(cigar) == "None":
                modcigar == ""
            else:
                modcigar = getfseq(cigar, seq)
            if (modcigar == ""):
                modcigar=0

            lengthFromFile = len(line.query_sequence)
            if lengthFromFile is None:
                lengthFromFile = int(lengthRead)

            #Add to allReadData
            Adata = '{"json":{"QNAME":"'+readName+'"'
            Adata += ', "FLAG":"'+str(flag)+'"'
            Adata += ', "RNAME":"'+str(chrom)+'"'
            Adata += ', "POS":"'+str(startLoc)+'"'
            Adata +=', "MAPQ":"'+str(line.mapping_quality)+'"'
            Adata += ', "CIGAR":"'+str(line.cigarstring)+'"'
            Adata +=', "RNEXT":"'+str(line.next_reference_id)+'"'
            Adata += ', "PNEXT":"'+str(line.next_reference_start+1)+'"'
            Adata += ', "TLEN":"'+str(line.template_length)+'"'
            Adata += ', "MODCIGAR":"'+str(modcigar)+'"'
            Adata += ', "OP":"'+str(line.get_tags())+'"'
            Adata += ', "endPOS":"'+str(startLoc + lengthFromFile+1)+'"'
            Adata += '}}'

            stream_num = startLoc/lengthBin + 1
            insertChrom = determineChr(chrom) #parse different chromosome types
            #check if data is good
            if (isReadWanted(flag) and insertChrom != -1):
                streamName= "chr"+str(insertChrom)+"stream"+str(stream_num)
                if((startLoc + lengthRead) > stream_num*lengthBin): #if the end of the key hangs into the next bin
                    key = "FLANK=1" #it hangs
                else:
                    key = "FLANK=0" #it doesn't hang
                keys = (key,)
                hasError = publishToStream(chainName, datadir, streamName, keys, multichainLoc, Adata)
            else:
                Adata = Adata[:-2]
                Adata += ', "SEQ:":"' + str(line.query_sequence)+'"'
                Adata += '}}'
                streamName = "unmappedANDcontigs"
                key1 = "readID="+readName
                key2 = "FLAG="+str(flag)
                keys=(str(key1), str(key2))
                hasError = publishToStream(chainName, datadir, streamName, keys, multichainLoc, Adata)
            currRead+=1
    bamfile.close()
    return hasError 


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-lc", "--lengthChromosome", type = int, help = "length of the chromosome in data", default = 249000000)
    #parser.add_argument("-numc", "--numChromosome", type = int, help = "number of chromosomes in data", default = 24)
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-lb", "--lengthBin", type = int, help = "length of each stream", default = 10000000)
    parser.add_argument("-lr", "--lengthRead", type = int, help = "max read length in the data", default = 100)
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    parser.add_argument("data", help = "data file(s) to store in chain", nargs='+')
    parser.add_argument("-rf", "--restart-from", type = int, help = "read #(s) to resume building chain from. separate with spaces respectively with data files", nargs='*')

    args = parser.parse_args()

    start = time.time()


    print("--PREPROCESSING--")
    num_streams = subscribeToStreams(args.chainName, 24, args.lengthChromosome, args.lengthBin, args.lengthRead, args.multichainLoc, args.datadir)
    if args.restart_from==None:
        args.restart_from = [0]*len(args.data)
        for data in args.data:
            if data[-6:-4]=='aa' or data[-6:-4]=='01' or len(args.data)==1:
                insertHeader(args.chainName, args.multichainLoc, data, args.datadir)

    print("--INSERTING DATA--")
    processes = []
    for pos, data in enumerate(args.data):
        p = multiprocessing.Process(target = classifyReads, 
            args = (args.chainName, num_streams, args.lengthBin, args.lengthRead, args.multichainLoc, data, args.datadir, args.restart_from[pos]))
        processes.append(p)
        p.start()
    for i, process in enumerate(processes):
        process.join()
        end = time.time()
        e = int(end - start)
        print('Time elapsed for process {}: '.format(i))
        print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

    print("Chain construction complete! Chain name: "+args.chainName)

    end = time.time()
    e = int(end - start)
    print('\n\n Total Time elapsed:\n\n')
    print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

    process = psutil.Process(os.getpid())
    print('\n\n Total memory in bytes:\n\n')
    print(process.memory_info().rss)

if __name__ == "__main__":
    main()
