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
#either translation file or pattern because the patterns contain the protid to taxid translation
group = parser.add_argument_group('group')
group.add_argument("-j", "--transfile", help="translation file between protein IDs and tax IDs")
group.add_argument("-p", "--pattern", help="file listing patterns including protein IDs, tax IDs and pattern IDs")
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
elif args.transfile:
    translationfile = args.transfile
else:
    print("A pattern file or translation file has to be given!\n")
    exit

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

if(args.pattern):
    with open(patternfile,'r', encoding='UTF-8') as file:
        for line in file:
            line.rstrip()
            if(line[0] == '>'):
                continue
            else:
                dat = line.split(" ")
                sdat = dat[0].split('-')
                protid2treeid[sdat[0]] = dat[0] #the pattern file already contains as ID the whole ID needed on the tree
elif(args.transfile):
    with open(translationfile,'r', encoding='UTF-8') as file:
        for line in file:
            line.rstrip()
            dat = line.split("\t")
            if(dat[7] == curcog):
                protid2treeid[dat[3]] = str(dat[3])+"-"+str(dat[0])
else:
    print(f'No translation or pattern file specified!')


#pattern file format:
#>PatternID 1; #Accessions 5; Structural_annotation NC_0-1=NC_1-1=K=L
#AAM02697.1-a190192-P1 (0,0,START) (4,134,7564.1.1.1) (143,195,187.1.1.9) (196,212,205.1.1.6) (213,226,7564.1.1.1) (241,251,END)
#WP_010876991.1-a187420-P1 (0,0,START) (1,132,7564.1.1.1) (144,192,187.1.1.9) (193,205,205.1.1.6) (206,220,7564.1.1.1) (231,241,END)
            

outfile = open(outputtree,"w")

pretree = read(inputtree, format="newick", into=TreeNode)
pretree.bifurcate()

for node in pretree.preorder():
    if(node.is_tip()): #leaves
        protids = node.name.split(' ')
        protid = '_'.join(protids)
        #print(f'node name {node.name}')
        #print(f'protid {protid}')
        if(protid in protid2treeid):
            #print(f'{protid2treeid[protid]}')
            node.name = str(protid2treeid[protid])

#outfile.write(str(pretree))
pretree.write(outfile, format='newick')
