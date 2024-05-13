import argparse
import pandas as pd
from io import StringIO
from skbio import read
from skbio.tree import TreeNode
import re


##input:
#translation table
#tree
#outputtree


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--intree")
parser.add_argument("-j", "--transfile")
parser.add_argument("-t", "--outtree")
args = parser.parse_args()

inputtree = ""
if args.intree:
    inputtree = args.intree
else:
    print("no inputtree given!\n")
    exit

translationfile= ""
if args.transfile:
    translationfile = args.transfile
else:
    print("no translation file given!\n")
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

with open(translationfile,'r', encoding='UTF-8') as file:
    for line in file:
        line.rstrip()
        dat = line.split("\t")
        if(dat[7] == curcog):
            protid2treeid[dat[3]] = str(dat[3])+"-"+str(dat[0])



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
