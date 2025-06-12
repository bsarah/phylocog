import argparse
import pandas as pd
from io import StringIO
from skbio import read
from skbio.tree import TreeNode
import re


##input:
#translation table
#list of patterns from structure tree such that the structural patterns can be plotted onto the sequence tree
#tree
#outputtree


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--intree", help="tree in newick format")
#use translation file AND pattern because the patterns contain the protid to taxid translation but there exist proteins that could not be mapped by HMMER (so no patterns)
parser.add_argument("-j", "--transfile", help="translation file between protein IDs and tax IDs")
parser.add_argument("-p", "--pattern", help="file listing patterns including protein IDs, tax IDs and pattern IDs")
parser.add_argument("-t", "--outtree", help="output tree in newick format")
args = parser.parse_args()

inputtree = ""
if args.intree:
    inputtree = args.intree
else:
    print("no inputtree given!\n")
    exit

patternfile = ""
translationfile= ""
if args.pattern:
    patternfile = args.pattern
else:
    print("No pattern file found, use translation!\n")

if args.transfile:
    translationfile = args.transfile
else:
    print("No translation file found, use pattern!\n")

outputtree = ""
if args.outtree:
    outputtree = args.outtree
else:
    outputtree = inputtree+"-relabeled.newick"


pre2cogid = inputtree.split('/')
precogid = pre2cogid[-1].split('.')
curcog = precogid[0]
print(f'cogid {curcog}')

    
protid2treeid = dict() #protein ids to lines


if(args.transfile):
    with open(translationfile,'r', encoding='UTF-8') as file:
        #use the pattern -PX for all the unmapped patterns
        for line in file:
            line.rstrip()
            dat = line.split("\t")
            if(dat[7] == curcog):
                protid2treeid[dat[3]] = str(dat[3])+"-"+str(dat[0])+"-PX"


if(args.pattern):
    with open(patternfile,'r', encoding='UTF-8') as file:
        for line in file:
            line.rstrip()
            if(line[0] == '>'):
                continue
            else:
                dat = line.split(" ")
 #               print(dat[0])
                #sdat = label without pattern ID
                presdat = dat[0].split('-')
                sdat = '-'.join(presdat[:-1])
 #               print(f'{sdat} vs {dat[0]}')
                #if the label only has the protein id, just take the first part and add the taxid and the patternid
                sdat2 = presdat[0]
                protid2treeid[sdat] = str(dat[0]) #the pattern file already contains as ID the whole ID needed on the tree
                protid2treeid[sdat2] = str(dat[0]) #the pattern file already contains as ID the whole ID needed on the tree

#pattern file format:
#>PatternID 1; #Accessions 5; Structural_annotation NC_0-1=NC_1-1=K=L
#AAM02697.1-a190192-P1 (0,0,START) (4,134,7564.1.1.1) (143,195,187.1.1.9) (196,212,205.1.1.6) (213,226,7564.1.1.1) (241,251,END)
#WP_010876991.1-a187420-P1 (0,0,START) (1,132,7564.1.1.1) (144,192,187.1.1.9) (193,205,205.1.1.6) (206,220,7564.1.1.1) (231,241,END)
            

outfile = open(outputtree,"w")

pretree = read(inputtree, format="newick", into=TreeNode)
pretree.bifurcate()

id = 0 #for postorder traversal
for node in pretree.postorder():
    if(node.is_tip()): #leaves
        protids = node.name.split()
        protid = '_'.join(protids)
        if(protid in protid2treeid):
            #            print(f'{protid} to {protid2treeid[protid]}')
            node.name = str(protid2treeid[protid])
            #there are quotes around the node names because: If your node names contain spaces or special characters, quotes will be added automatically for Newick format validity
        else:
            #            print(f'{protid} not in protid2treeid')
            node.name = node.name+"-PX" #add X for unknown pattern

    else:
        if(node.is_root()):
            node.name = f'r{str(id)}'
        else:
            node.name = f'i{str(id)}i'
    id+=1

#outfile.write(str(pretree))
pretree.write(outfile, format='newick')
