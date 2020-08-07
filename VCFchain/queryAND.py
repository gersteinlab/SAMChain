'''
queryAND.py
Query number of reads at a specific location from a chain
Usage: $ python queryAND.py <chainName> <readLocation> -ml/--multichainLoc <multichainLoc>
Example: python queryAND.py chain1 chr1:1005-2005
modified by CMB on 06/01/2020
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


#get list of all items stored under a given key in a given stream
def queryKey (chainName, streamName, key, datadir, multichainLoc):
	queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeyitems {} {} false 100000000'.format(chainName, datadir, streamName, key)
	items = subprocess.check_output(queryCommand.split())
	if("error" in items or "ERROR" in items):
		return -1
	matches = json.loads(items, parse_int= int)
	return matches
	
#get stream name to search in
def getStreamName(readLoc):
    streamName=str(readLoc.split(":")[0])
    return streamName

#get the list of keys stored in a given stream and keep those that match the query
def getLocKey(chainName, readLoc, rsidInput, gtInput, datadir, multichainLoc, streamName, numChromosome):
    loc = readLoc.split(":")[1].split("-")
    wanted_keys = []
    final_keys = []
    queryCommand=multichainLoc+'multichain-cli {} -datadir={} liststreamkeys {}'.format(chainName, datadir, streamName)
    items = subprocess.check_output(queryCommand.split())
    if("error" in items or "ERROR" in items):
        return -1
    got_keys = json.loads(items, parse_int= int)
    start = ""
    end = ""
    name = ""
    for entry in got_keys:
        key = entry['key']
        pos_key = key.split('=')[1]
        int_pos_key = 0
        try:
            int_pos_key = int(pos_key)
        except ValueError:
            print("Not a valid position")
        if len(loc) > 1:
            if int_pos_key >= int(loc[0]) and int_pos_key <= int(loc[1]):
                wanted_keys.append(str(pos_key))
        elif len(loc)==1:
            #print(int_pos_key, loc[0])
            if int_pos_key== int(loc[0]):
                wanted_keys.append(str(pos_key))
    if (rsidInput=="any" and gtInput=="any"):
        return(wanted_keys)
    

    for elem in wanted_keys:
        ANDgt = False
        ANDrsid = False
        if (rsidInput == "any"):
            ANDrsid = True
        if (gtInput == "any"):
            ANDgt = True

        entry = queryKey(chainName, streamName, "POS="+str(elem), datadir, multichainLoc)
        gt = entry[0]["keys"][1].split("=")[1]
        rsid = entry[0]["keys"][2].split("=")[1]
        if (gt == gtInput):
            ANDgt = True
        if (rsid == rsidInput):
            ANDrsid = True
        if (ANDgt==True and ANDrsid==True):
            final_keys.append(elem)
        print(gt, rsid)
    print(final_keys)
    return final_keys

#query the allReadData stream with the matching keys
def queryAllData(multichainLoc, chainName, dr, key):
    output=''
    for name in key:
        matches = queryKey(chainName, 'allVariantData', name, dr, multichainLoc)
        for result in matches:
            chrom=result["data"]["json"]['CHROM']
            pos=result["data"]["json"]['POS']
            record=result["data"]["json"]['ID']
            ref=result["data"]["json"]['REF']
            alt=result["data"]["json"]['ALT']
            qual=result["data"]["json"]['QUAL']
            filt=result["data"]["json"]['FILTER']
            info=result["data"]["json"]['INFO']
            form=result["data"]["json"]['FORMAT']
            samples=result["data"]["json"]["SAMPLEs"]
            
            if bool(re.search('^chr[1-9]$|^chr[1-2][0-9]$|^chrX$|^chrY$', str(chrom))): #ignores None cigars and complex chr names
                #reformat alt string for output
                a = alt.replace('(','')
                b = a.replace(',','')
                c = b.replace(')','')
                d = c.replace('\'', '')
                #reformat filter string for output
                a1 = filt.replace('(','')
                b1 = a1.replace(')','')
                c1 = b1.replace('\'','')
                d1 = c1.replace('[','')
                e1 = d1.replace(']','')
                #reformat info string for output
                a2 = info.replace('\'','')
                b2 = a2.replace(' ','')
                c2 = b2.replace('{', '')
                d2 = c2.replace('}', '')
                e2 = d2.replace('(','')
                f2 = e2.replace(',)','')
                g2 = f2.replace(':','=')
                h2 = g2.replace(',',';')
                #reformat format string for output
                a3 = form.replace('[','')
                b3 = a3.replace(']','')
                c3 = b3.replace(' ', '')
                d3 = c3.replace(',', ':')
                e3 = d3.replace('\'', '')
                f3 = e3.split(":")
                #reformat samples string for output
                a4 = samples.replace('{','')
                b4 = a4.replace('}','')
                for elem in f3:
                    b4 = b4.replace(elem,'')
                c4 = b4.replace('\'','')
                d4 = c4.replace('(','')
                e4 = d4.replace(')','')
                f4 = e4.replace(', :', ';') 
                g4 = f4.replace(":", "")
                h4 = g4.replace(' ','')
                i4 = h4.replace(";", ":")
                j4 = i4.replace("0,1", "0/1")
                k4 = j4.replace("1,1", "1/1")
                
                output = output+chrom+"\t"+pos+"\t"+record+"\t"+ref+"\t"+d+"\t"+qual+"\t"+e1+"\t"+h2+"\t"+e3+"\t"+k4 + "\n"
    print(output)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("chainName", help = "the name of the chain storing the data")
	parser.add_argument("readLoc", help = "the read location to query")
        parser.add_argument("-rsid", help = "the rsid to query", type = str, default = "any")
        parser.add_argument("-gt", help="the genotype to query", type = str, default = "any")
	parser.add_argument("-ml", "--multichainLoc", help = "path to multichain commands", default = "")
	parser.add_argument("-dr", "--dir", help= "chain location")
        parser.add_argument("-numc", "--numChromosome", type = int, help = "number of chromosomes in data", default = 24)
        args = parser.parse_args()
 
        streamName = getStreamName(args.readLoc) 
        key = getLocKey(args.chainName, args.readLoc, args.rsid, args.gt, args.dir, args.multichainLoc, streamName, args.numChromosome) 
        data = queryAllData(args.multichainLoc, args.chainName, args.dir, key)
	
        if(data == 0): #read was out of bounds
	    exit(1)


if __name__ == "__main__":
	main()
