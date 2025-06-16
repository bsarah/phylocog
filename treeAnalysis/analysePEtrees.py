import argparse
#import pandas as pd
#from io import StringIO
from skbio import read
from skbio.tree import TreeNode
#import subprocess
#import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--intree", help="tree in newick format with node labels, sankoff output")
parser.add_argument("-c", "--cogid", help="COGid for the output")
parser.add_argument("-f", "--fcat", help="functional category")
parser.add_argument("-a", "--asciitree", help="print ascii tree")
#parser.add_argument("-o", "--outtable", help="table with calculates tree scores")


args = parser.parse_args()

inputtree = ""
if args.intree:
    inputtree = args.intree
else:
    print("no inputtree given!\n")
    exit

cogid = "COGXXXX"
fcat = "XX"

if(args.cogid):
    cogid = args.cogid

if(args.fcat):
    fcat = args.fcat

asciitree = False
if(args.asciitree):
    asciitree = True
    
#outputfile = ""
#if args.outtable:
#    outputfile = args.outtable
#else:
#    print("no output table given!\n")
#    exit




#output file header:
#COG numa numb numx numtotal numpattern funcCat numsplits S_root delta_S_splits delta_S_nonsplit



#read tree in postorder
#label nodes bottom-up (A/B/W/S)
#->modified Fitch algorithm

#top down approch:
#label nodes with (A/B/S)
#collect root score
#collect scores for split branches and non-split branches
#count splits

pretree = read(inputtree, format="newick", into=TreeNode)
pretree.bifurcate()
#print(pretree)
#print(pretree.ascii_art())
nodenum = pretree.count()
node2labels = dict() #nodename to set of labels (a,b,e,x)
patterns= set()
numa = 0
numb = 0
nume = 0
numx = 0
#leaf-label: protid-[a|b|e|x]num-Pnum
#internal label: nodeid-score-Pnum

for node in pretree.postorder():
    #print(f'curid {curid}')
    nns = node.name.split('-')
    pattern = nns[-1]
    patterns.add(pattern)
    if(node.is_tip()):
        if("-a" in node.name):
            node2labels[node.name] = {'a'}
            numa += 1
        elif("-b" in node.name):
            node2labels[node.name] = {'b'}
            numb += 1
        elif("-e" in node.name):
            node2labels[node.name] = {'e'}
            nume += 1
        else:
            node2labels[node.name] = {'x'}
            numx += 1
    else:
        child1 = node.children[0]
        c1name = child1.name
        c1labels = node2labels[c1name]
        c2labels = set()
        if(len(node.children) > 1):
            child2 = node.children[1]
            c2name = child2.name
            c2labels = node2labels[c2name]
        #take intersection if non-empty
        #otherwise take union
        c1c2ints = c1labels.intersection(c2labels)
        if(len(c1c2ints) > 0):
            node2labels[node.name] = c1labels.intersection(c2labels)
        else:
            node2labels[node.name] = c1labels.union(c2labels)


#go in preorder and label nodes definitively, collect score differences
node2deflabel = dict() #node name to one letter (a,b,e,s)
numsplits = 0 #number of split branches
numbranches = 0 #total number of branches
sumsplitscores = 0
sumnonsplitscores = 0
rootscore = 0

#leaf-label: protid-[a|b|e|x]num-Pnum
#internal label: nodeid-score-Pnum
for node in pretree.preorder():
    nns = node.name.split('-')
    if(node.is_root()):
        rootscore = float(nns[1])
        rootls = node2labels[node.name]
        #take the first label of the set
        node2deflabel[node.name] = list(rootls)[0]
    else: #internal or leaf
        numbranches+=1
        #check if the parent node has the same label
        ownlabs = node2labels[node.name]
        parentname = node.parent.name
        ps = parentname.split('-')
        parentscore = float(ps[1])
        parentlab = node2deflabel[parentname]
        if(node.is_tip()):
            ownscore = 0
        else:
            ownscore = float(nns[1])
        deltas = parentscore-ownscore
        #check if parentlab in ownlabs
        #if yes, take the parent label-->non-split branch
        #otherwise, take your own label-->split branch
        if(parentlab in ownlabs):
            #no split nodes
            sumnonsplitscores += deltas
            node2deflabel[node.name] = parentlab
        else:
            #split branch
            sumsplitscores += deltas
            numsplits+=1
            node2deflabel[node.name] = list(ownlabs)[0]

if(asciitree):
    print(pretree.ascii_art())
            
#output file header:
#COG numa numb numx numtotal numpattern funcCat numsplits S_root delta_S_splits delta_S_nonsplit
avdeltasplits = 0
if(numsplits > 0):
    avdeltasplits = round(sumsplitscores/numsplits,2)
nonsplitbranches = numbranches-numsplits
avdeltanonsplits = round(sumnonsplitscores/nonsplitbranches,2)
#print(f'COGid\tnuma\tnumb\tnume\tnumx\tnumtotal\tnumpattern\tfcat\tnumsplits\ts_root\tavdelta_splits\tavdelta_nonsplits')
print(f'{cogid}\t{numa}\t{numb}\t{nume}\t{numx}\t{numa+numb+nume+numx}\t{len(patterns)}\t{fcat}\t{numsplits}\t{rootscore}\t{avdeltasplits}\t{avdeltanonsplits}')


#TODO:
#output file
#add cogid and fcat

#check else cases for -P and -(a|b|e|x) labels, sankoff will report an error if it doesn't find any -P label (sometimes, the - is missing, only PX
