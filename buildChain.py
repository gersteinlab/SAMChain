'''
buildChain.py
Loads tab-separated text file data onto an existing data stream on an existing multichain
Usage: $ python buildChain.py <chainName> <length_chromosome> <length_bin> <length_read> <tab-separated-data>.bam
Example: python buildChain.py -cn chain1 -lc 4050 -lb 1000 -lr 100 /home/slw67/fake_data/length_startloc.txt
Example: python buildChain.py /home/slw67/nextgen-blockchain/fake_w_chr.txt
Example: python buildChain.py /home/slw67/nextgen-blockchain/test.sorted.bam -t bam
Example: python buildChain.py /home/slw67/nextgen-blockchain/fake_w_chr.txt
modified by SW on 2/22/2019
'''

import pysam
import sys
import time
import math
import binascii
import argparse
import subprocess
from subprocess import Popen, PIPE

#define global variables
chrType = -1 #1 if of type {1, 2, ...X |Y}, 2 if of type {chr1, chr2, ... chrX|chrY}

#Creates a chain given the name of the chain and the location of the multichain program
def createChain(chainName, multichainLoc, datadir):
	hasError = 0
	createCommand=multichainLoc+'multichain-util create {} -datadir={}'.format(chainName, datadir)
	runCommand = multichainLoc+'multichaind {} -datadir={} -daemon'.format(chainName, datadir)

	#make the chain
	procCreate = subprocess.Popen(createCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout_valc, stderr_valc = procCreate.communicate()
	if(("error" in stderr_valc) or ("ERROR" in stderr_valc)):
		print("Could not create chain (probably chain already exists)")
		hasError = -1
	else:
		print("Created chain")

	#start the chain
	daemonOut = subprocess.call(runCommand.split()) #returns 0 whether or not it had an error; subprocess hangs 
	#because output is too long (known subprocess bug); just relying on chain creation to catch bug
	time.sleep(1)
	return hasError

#Given a chain name and the name of the new stream, makes that stream
def makeStream(chainName, streamName, multichainLoc, datadir):
	createStreamCommand=multichainLoc+'multichain-cli {} -datadir={} create stream {} true'.format(chainName,datadir,streamName)
	subscribeStreamCommand=multichainLoc+'multichain-cli {} -datadir={} subscribe {}'.format(chainName,datadir,streamName)

	#make stream of name StreamName
	procCreate = subprocess.Popen(createStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	procCreate.wait()

	#subscribe to the stream
	procSubscribe = subprocess.Popen(subscribeStreamCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	procCreate.wait()

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
#all read data has unmapped reads
def createStreams(chainName, numChromosome, lengthChromosome, lengthBin, lengthRead, multichainLoc, datadir):
	#data stream by readName
	makeStream(chainName, "metaData", multichainLoc, datadir) #header, binLength, other things for the whole file
        makeStream(chainName, "allReadData", multichainLoc, datadir) #stores all the data for each line
        makeStream(chainName, "unmappedANDcontigs", multichainLoc, datadir)
	#hash for location
	num_streams = int(math.ceil(float(lengthChromosome)/float(lengthBin)))
	#range not inclusive, so run num_streams times
	for c in range(1, numChromosome+1):
		for i in range(1, num_streams+1):
			makeStream(chainName, "chr"+str(c)+"stream"+str(i), multichainLoc, datadir)

	#data for chain:
	chainData = '{"json":{"binLength":'+str(lengthBin)+', "readLength":'+str(lengthRead)+', "numBins":'+str(num_streams)+'}}'
	publishToStream(chainName, datadir, "metaData", ("chainInfo",), multichainLoc, chainData)
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
	allFlags = list("{0:b}".format(flag))[::-1]
	dropFlags = [2, 8, 9, 10, 11]
        return True

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
def classifyReads(chainName, num_streams, lengthBin, lengthRead, multichainLoc, data, datadir):
	hasError = 0


	with pysam.AlignmentFile(data, "rb") as bamfile:
		num_dups = 0
                for line in bamfile:
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
			Adata += '}}'
                        #NM_val = 0
                        #for elem in line.get_tags():
                        #    if elem[0] == "NM":
                        #        NM_val = elem[1]

                        #print("NM val is " + str(NM_val))
                        #if NM_val != 0:
                        #    readNameKey = readName + "." + str(num_dups)
                        #    num_dups = num_dups + 1
                        #else:
                        readNameKey = readName
                        keys = (readNameKey,)
                        hasError = publishToStream(chainName, datadir, "allReadData", keys, multichainLoc, Adata)

			#Add to location query stream (format: chr1stream1)
                        #bin the positions
			lengthFromFile = len(line.query_sequence)
			if lengthFromFile is None:
				lengthFromFile = int(lengthRead)
			stream_num = startLoc/lengthBin + 1
			insertChrom = determineChr(chrom) #parse different chromosome types
			#check if data is good
			if (isReadWanted(flag) and insertChrom != -1):
				streamName= "chr"+str(insertChrom)+"stream"+str(stream_num)
				key1 = "readID="+readName
				key2 = "POS="+str(startLoc)
				key4 = "endPOS="+str(startLoc+lengthFromFile+1)
				if((startLoc + lengthRead) > stream_num*lengthBin): #if the end of the key hangs into the next bin
					key3 = "FLANK=1" #it hangs
				else:
					key3 = "FLANK=0" #it doesn't hang
                                key5 = "CIGAR="+str(line.cigarstring)
                                key6 = "MAPQ="+str(line.mapping_quality)
				keys = (str(key1), str(key2), str(key3), str(key4), str(key5), str(key6))
				hasError = publishToStream(chainName, datadir, streamName, keys, multichainLoc)
                        if (insertChrom == -1):
                                streamName = "unmappedANDcontigs"
                                key1 = "readID="+readName
                                key2 = "FLAG="+str(flag)
                                keys=(str(key1), str(key2))
                                hasError = publishToStream(chainName, datadir, streamName, keys, multichainLoc)
	bamfile.close()
	return hasError	

#for test reasons from plaintext
def classifyReadsTest(chainName, num_streams, lengthBin, lengthRead, multichainLoc, data):
	hasError = 0
	try:
		lines = open(data).read().splitlines()
	except IOError:
		print "Error: File does not appear to exist."
		return -1

	for line in lines:
		vals = line.split('\t')
		readName = vals[0]

		#Add to allReadData
		keys = (readName,)
		hasError = publishToStream(chainName, datadir, "allReadData", keys, multichainLoc)

		#Add to location query stream (format: chr1stream1)
		flag = int(vals[1])
		chrom = vals [2]
		startLoc = int(vals[3])
		lengthRead = int(lengthRead)
		stream_num = startLoc/lengthBin + 1
		insertChrom = determineChr(chrom) #parse different chromosome types
		#check if data is good
		if (isReadWanted(flag) and insertChrom != -1):
			streamName= "chr"+str(insertChrom)+"stream"+str(stream_num)
			key1 = "readID="+readName
			key2 = "POS="+str(startLoc)
			if((startLoc + lengthRead) > stream_num*lengthBin): #if the end of the key hangs into the next bin
				key3 = "FLANK=1" #it hangs
			else:
				key3 = "FLANK=0" #it doesn't hang
			keys = (str(key1), str(key2), str(key3))
			#data = '{"json":{"readID":"'+readName+'"}}'
			hasError = publishToStream(chainName, datadir, streamName, keys, multichainLoc)
    
	return hasError

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("data", help = "data file to store in chain")
	parser.add_argument("-lc", "--lengthChromosome", type = int, help = "length of the chromosome in data", default = 249000000)
	parser.add_argument("-numc", "--numChromosome", type = int, help = "number of chromosomes in data", default = 24)
	parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
	parser.add_argument("-lb", "--lengthBin", type = int, help = "length of each stream", default = 10000000)
	parser.add_argument("-lr", "--lengthRead", type = int, help = "max read length in the data", default = 100)
	parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
        parser.add_argument("-dr", "--datadir", help = "path to store the chain")
	args = parser.parse_args()

	#make a chain
	print("--CHAIN CREATION--")
	hasError = createChain(args.chainName, args.multichainLoc, args.datadir)
	if hasError:
		sys.stderr.write("\nERROR: Failed chain creation. Please try again with a different chain name.\n")
		quit()
	
	#make the streams
	print("--STREAM CREATION--")
	num_streams = createStreams(args.chainName, 23, args.lengthChromosome, args.lengthBin, args.lengthRead, args.multichainLoc, args.datadir)
    
	print("--PREPROCESSING--")
	insertHeader(args.chainName, args.multichainLoc, args.data, args.datadir)

	print("--INSERTING DATA--")
	classifyReads(args.chainName, num_streams, args.lengthBin, args.lengthRead, args.multichainLoc, args.data, args.datadir)

	print("Chain construction complete! Chain name: "+args.chainName)

        import os
        import psutil
        process = psutil.Process(os.getpid())
        print('\n\n Total memory in bytes:\n\n')
        print(process.memory_info().rss)

if __name__ == "__main__":
	main()
