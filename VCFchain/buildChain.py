'''
buildChain.py
Loads tab-separated text file data onto an existing data stream on an existing multichain
Usage: $ python buildChain.py <chainName> <length_chromosome> <length_bin> <length_read> <tab-separated-data>.vcf
modified by CMB on 05/28/2020
'''

from pysam import VariantFile
import sys
import re
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
#Different keys stored for SAM vs. VCF
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
        if keys == "header" or keys == "chainInfo":
            keyString=keys

        publishCommand = [multichainLoc+'multichain-cli',
        str('{}'.format(chainName)),
        str('-datadir={}'.format(datadir)),
        'publish',
        str('{}'.format(streamName)),
        str('{}'.format(keyString)),
        str('{}'.format(data))]
        dummy = subprocess.check_output(publishCommand, stderr=subprocess.STDOUT)

        return hasError


#Making the streams that will be used to store the file
def createChrStreams(chainName, numChromosome, multichainLoc, datadir):
        #data stream by readName
        makeStream(chainName, "metaData", multichainLoc, datadir) #stores header, and can store other things for the whole file
        makeStream(chainName, "allVariantData", multichainLoc, datadir) #stores all the data for each line
        #range not inclusive, so run num_streams times
        for c in range(1, numChromosome+1):
            makeStream(chainName, "chr"+str(c), multichainLoc, datadir)


#Same for SAM and VCF
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


def isReadWanted (qual):
        return True #can threshold by qual score here if desired


def insertHeader(chainName, multichainLoc, data, datadir):
	vcf_in = VariantFile(data)
        header = vcf_in.header
        hed = str(header)
        hed_esc = re.sub(r'(["])', r'\"', hed)
        ss=hed_esc.splitlines(True)

        for i in range(0, len(ss)):
	    headerData = '{"json":{"header":"'+str(ss[i])+'"}}'
            publishToStream(chainName, datadir, "metaData", "header", multichainLoc, headerData)


#binning and hashing happens here
#we currently bin by location
def classifyReads(chainName, multichainLoc, data, datadir):
	hasError = 0
	vcf_in = VariantFile(data)	
        num_dups = 0
        for line in vcf_in.fetch():
                chrom = line.chrom
                pos = line.pos
                record = line.id
                if str(record)=="None":
                    record="."
                ref = line.ref
                alt = line.alts
                qual = line.qual
                filt = list(line.filter)
                info = dict(line.info)
                form = list(line.format)
                key = list(line.samples)[0]
                sample = dict(line.samples[key])
                gt = dict(line.samples[key])['GT']
                if str(gt) == "(0, 1)":
                    gt = "1"
                elif str(gt) == "(1, 1)":
                    gt = "2"

                #Add to allReadData
                Adata = '{"json":{"CHROM":"'+str(chrom)+'"'
                Adata += ', "POS":"'+str(pos)+'"'
                Adata += ', "GT":"'+str(gt)+'"'
                Adata += ', "ID":"'+str(record)+'"'
                Adata += ', "REF":"'+str(ref)+'"'
                Adata +=', "ALT":"'+str(alt)+'"'
                Adata += ', "QUAL":"'+str(qual)+'"'
                Adata +=', "FILTER":"'+str(filt)+'"'
                Adata += ', "INFO":"'+str(info)+'"'
                Adata += ', "FORMAT":"'+str(form)+'"'
                Adata += ', "SAMPLEs":"'+str(sample)+'"'
                Adata += '}}'

                posKey = str(pos)
                gtKey = str(gt)
                idKey = str(record)
                keys = (posKey, gtKey, idKey)
                hasError = publishToStream(chainName, datadir, "allVariantData", keys, multichainLoc, Adata)

                #Add to location query stream (format: chr1stream1)
                #bin the positions
                insertChrom = determineChr(chrom) #parse different chromosome types
                #check if data is good
                if (isReadWanted(qual) and insertChrom != -1):
                        streamName = "chr"+str(insertChrom)
                        key1 = "POS="+posKey
                        key2 = "GT="+gtKey
                        key3 = "ID="+idKey
                        keys = (str(key1), str(key2), str(key3))
                        hasError = publishToStream(chainName, datadir, streamName, keys, multichainLoc)

	return hasError	


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("data", help = "data file to store in chain")
	parser.add_argument("-lc", "--lengthChromosome", type = int, help = "length of the chromosome in data", default = 249000000)
	parser.add_argument("-numc", "--numChromosome", type = int, help = "number of chromosomes in data", default = 24)
	parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
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
        createChrStreams(args.chainName, args.numChromosome,  args.multichainLoc, args.datadir)

	print("--PREPROCESSING--")
	insertHeader(args.chainName, args.multichainLoc, args.data, args.datadir)

	print("--INSERTING DATA--")
        classifyReads(args.chainName, args.multichainLoc, args.data, args.datadir)

	print("Chain construction complete! Chain name: "+args.chainName)

        import os
        import psutil
        process = psutil.Process(os.getpid())
        print('\n\n Total memory in bytes:\n\n')
        print(process.memory_info().rss)

if __name__ == "__main__":
	main()
