'''
buildChain.py
Initializes an empty blockchain for storing tab-separated text file data (SAM or BAM format) 
Usage: $ python buildChain.py <chainName> <length_chromosome> <length_bin> <length_read> 
Example: python buildChain.py -cn chain1 -lc 4050 -lb 1000 -lr 100 
Example: python buildChain.py -cn chain1 
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
from tqdm import tqdm
import os
import psutil
import time
#define global variables
chrType = -1


#create a chain given the name of the chain, the location of the multichain program, and the datadir in which to keep the chain
def createChain(chainName, multichainLoc, datadir):
    hasError = 0
    createCommand=multichainLoc+'multichain-util create {} -datadir={}'.format(chainName, datadir)
    runCommand = multichainLoc+'multichaind {} -datadir={} -daemon'.format(chainName, datadir)

    #make the chain
    procCreate = subprocess.Popen(createCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_valc, stderr_valc = procCreate.communicate()
    if(("error" in stderr_valc) or ("ERROR" in stderr_valc)):
        print("Could not create chain (it probably already exists)")
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
def createStreams(chainName, numChromosome, lengthChromosome, lengthBin, lengthRead, multichainLoc, datadir):
    #data stream by readName
    makeStream(chainName, "metaData", multichainLoc, datadir) #header, binLength, other things for the whole file
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



######################################################################################################




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-lc", "--lengthChromosome", type = int, help = "length of the chromosome in data", default = 249000000)
    #parser.add_argument("-numc", "--numChromosome", type = int, help = "number of chromosomes in data", default = 24)
    parser.add_argument("-cn", "--chainName", help = "the name of the chain to store data", default = "chain1")
    parser.add_argument("-lb", "--lengthBin", type = int, help = "length of each stream", default = 10000000)
    parser.add_argument("-lr", "--lengthRead", type = int, help = "max read length in the data", default = 100)
    parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
    parser.add_argument("-dr", "--datadir", help = "path to store the chain")
    args = parser.parse_args()

    start = time.time()
    #make a chain
    print("--CHAIN CREATION--")
    hasError = createChain(args.chainName, args.multichainLoc, args.datadir)
    if hasError:
        sys.stderr.write("\nERROR: Failed chain creation. Please try again with a different chain name.\n")
        quit()

    #make the streams
    print("--STREAM CREATION--")
    num_streams = createStreams(args.chainName, 24, args.lengthChromosome, args.lengthBin, args.lengthRead, args.multichainLoc, args.datadir)
    
    print("Chain construction complete! Chain name: "+args.chainName)

    end = time.time()
    e = int(end - start)
    print('\n\n Time elapsed:\n\n')
    print( '{:02d}:{:02d}:{:02d}'.format(e // 3600, (e % 3600 // 60), e % 60))

    process = psutil.Process(os.getpid())
    print('\n\n Total memory in bytes:\n\n')
    print(process.memory_info().rss)

if __name__ == "__main__":
    main()
