import argparse
#import pandas as pd
#from io import StringIO
from skbio import read
from skbio.tree import TreeNode
import numpy as np
#import subprocess
#import os



def calcABsplits(inputtree,cogid, fcat,asciitree):
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

    #print(node2labels)
    #go in preorder and label nodes definitively, collect score differences
    node2deflabel = dict() #node name to one letter (a,b,e,s)
    numsplits = 0 #number of split branches/outside clade
    numbranches = 0 #total number of branches
    sumsplitscores = 0#outside clade score diff
    sumnonsplitscores = 0#inside clade score diff
    rootscore = 0
    
    scorelist = [] #sankoff scores in preorder
    cnumlist = []#num different clades in preorder
    
    #leaf-label: protid-[a|b|e|x]num-Pnum
    #internal label: nodeid-score-Pnum
    for node in pretree.preorder():
        nns = node.name.split('-')
        if(node.is_root()):
            rootscore = float(nns[1])
            rootls = node2labels[node.name]
            #take the first label of the set
            node2deflabel[node.name] = list(rootls)[0]
            #print(node2union[node.name])
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
            #only append if cnum changed between parent and current?
            #deltas = parentscore-ownscore
            #check if parentlab in ownlabs
            #if yes, take the parent label-->non-split branch
            #otherwise, take your own label-->split branch
            if(parentlab in ownlabs):
                #no split nodes
                sumnonsplitscores += ownscore
                node2deflabel[node.name] = parentlab
            else:
                #split branch
                sumsplitscores += ownscore
                numsplits+=1
                node2deflabel[node.name] = list(ownlabs)[0]

            #print(node2deflabel)
    if(asciitree):
        print(pretree.ascii_art())

    #output file header:
    #COG numa numb numx numtotal numpattern funcCat numsplits S_root delta_S_splits delta_S_nonsplit
    avdeltasplits = 0
    if(numsplits > 0):
        avdeltasplits = round(sumsplitscores/numsplits,5)
    nonsplitbranches = numbranches-numsplits
    avdeltanonsplits = round(sumnonsplitscores/nonsplitbranches,5)
    #print(f'COGid\tnuma\tnumb\tnume\tnumx\tnumtotal\tnumpattern\tfcat\tnumsplits\ts_root\tavdelta_splits\tavdelta_nonsplits')
    return (f'{cogid}\t{numa}\t{numb}\t{nume}\t{numx}\t{numa+numb+nume+numx}\t{len(patterns)}\t{fcat}\t{numsplits}\t{rootscore}\t{avdeltasplits}\t{avdeltanonsplits}')




def calcCladesplits(inputtree, cogid, fcat,taxid2clade,asciitree,both=False):
    #read tree in postorder
    #node labels = clades as in taxid2clade
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
    node2union = dict() #nodename to set of all clades in subtree
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
            pretaxid = nns[-2]
            taxid = pretaxid[1:] #remove a,b,e,x
            if(taxid in taxid2clade):
                node2labels[node.name] = {taxid2clade[taxid]}
                node2union[node.name] = {taxid2clade[taxid]}
                # print(taxid2clade[taxid])
            else:
                print(f'taxid {taxid} of {node.name} not existing\n')
                node2labels[node.name] = {"OTHER"}
                node2union[node.name] = {"OTHER"}
            if("-a" in node.name):
                numa+=1
            elif("-b" in node.name):
                numb+=1
            elif("-e" in node.name):
                nume+=1
            else:
                numx+=1
        else: #inner nodes
            child1 = node.children[0]
            c1name = child1.name
            c1labels = node2labels[c1name]
            c1union = node2union[c1name]
            c2labels = set()
            if(len(node.children) > 1):
                child2 = node.children[1]
                c2name = child2.name
                c2union = node2union[c2name]
                c2labels = node2labels[c2name]
                #take intersection if non-empty
                #otherwise take union
            c1c2ints = c1labels.intersection(c2labels)
            node2union[node.name] = c1union.union(c2union)
            if(len(c1c2ints) > 0):
                node2labels[node.name] = c1labels.intersection(c2labels)
            else:
                node2labels[node.name] = c1labels.union(c2labels)

    #print(node2labels)
    #go in preorder and label nodes definitively, collect score differences
    node2deflabel = dict() #node name to one letter (a,b,e,s)
    numsplits = 0 #number of split branches/outside clade
    numinclade = 0#number of non-split/inside clade
    numbranches = 0 #total number of branches
    sumsplitscores = 0#outside clade score diff
    sumnonsplitscores = 0#inside clade score diff
    rootscore = 0
    rootclades = 0

    #scorelist = [] #sankoff scores in preorder
    #cnumlist = []#num different clades in preorder

    #leaf-label: protid-[a|b|e|x]num-Pnum
    #internal label: nodeid-score-Pnum
    for node in pretree.preorder():
        nns = node.name.split('-')
        if(node.is_root()):
            rootscore = float(nns[1])
            rootls = node2labels[node.name]
            #take the first label of the set
            node2deflabel[node.name] = list(rootls)[0]
            rootclades = len(node2union[node.name])
            #print(node2union[node.name])
            #scorelist.append(rootscore)
            #cnumlist.append(rootclades)
        else: #internal or leaf
            numbranches+=1
            #check if the parent node has the same label
            ownlabs = node2labels[node.name]
            parentname = node.parent.name
            ps = parentname.split('-')
            parentscore = float(ps[1])
            parentlab = node2deflabel[parentname]
            pcnum = len(node2union[parentname])
            if(node.is_tip()):
                ownscore = 0
                piesize = 1
            else:
                ownscore = float(nns[1])
                piesize = len(node2union[node.name])
                #only append if cnum changed between parent and current?
                #what is the baseline to compare to?
            deltas = parentscore-ownscore
            if(pcnum > piesize): #outside clade
                sumsplitscores += ownscore
                numsplits+=1
                #scorelist.append(ownscore/rootscore)
                #cnumlist.append(len(node2union[node.name])/rootclades)
            else: #inside clade
                sumnonsplitscores+=ownscore
                numinclade+=1
            
            #check if parentlab in ownlabs
            #if yes, take the parent label-->non-split branch
            #otherwise, take your own label-->split branch
            if(parentlab in ownlabs):
                #no split nodes
                #sumnonsplitscores += ownscore
                node2deflabel[node.name] = parentlab
            else:
                #split branch
                #sumsplitscores += ownscore
                #numsplits+=1
                node2deflabel[node.name] = list(ownlabs)[0]

    #print(node2deflabel)
    if(asciitree):
        print(pretree.ascii_art())

    #correlation = np.corrcoef(scorelist, cnumlist)[0, 1]

    #output file header:
    #COG numa numb numx numtotal numpattern funcCat numsplits S_root delta_S_splits delta_S_nonsplit
    avdeltasplits = 0
    if(numsplits > 0):
        avdeltasplits = round(sumsplitscores/numsplits,5)
    avdeltanonsplits = round(sumnonsplitscores/numinclade,5)
    #print(f'COGid\tnuma\tnumb\tnume\tnumx\tnumtotal\tnumpattern\tfcat\tnumsplits\ts_root\tavdelta_splits\tavdelta_nonsplits')
    if(both):
        return (f'{avdeltasplits}\t{avdeltanonsplits}\t{rootclades}')
    return (f'{cogid}\t{numa}\t{numb}\t{nume}\t{numx}\t{numa+numb+nume+numx}\t{len(patterns)}\t{fcat}\t{numsplits}\t{rootscore}\t{avdeltasplits}\t{avdeltanonsplits}\t{rootclades}')







def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--intree", help="tree in newick format with node labels, sankoff output")
    parser.add_argument("-c", "--cogid", help="COGid for the output")
    parser.add_argument("-f", "--fcat", help="functional category")
    parser.add_argument("-s", "--clades", help="translate tax ids to subclades and visualize them")
    parser.add_argument("-b", "--both", help="get clade and domain scores",action='store_true')
    parser.add_argument("-a", "--asciitree", help="print ascii tree",action='store_true')
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


    clades = ""
    if(args.clades):
        clades = args.clades
    
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
    if(args.clades or args.both):
        taxid2clade = dict()
        with open(clades,'r') as file:
            for line in file:
                l2 = line.strip()
                ls = l2.split(',')
                taxid = ls[-2]
                cl = ls[-1]
                taxid2clade[taxid] = cl
        #call function to calculate clade score
        outstr = calcCladesplits(inputtree, cogid, fcat,taxid2clade,asciitree,args.both)
        if(args.both):
            tmpout = calcABsplits(inputtree, cogid, fcat,asciitree)
            realout = f'{tmpout}\t{outstr}'
            outstr = realout
    else:
        outstr = calcABsplits(inputtree, cogid, fcat,asciitree)
    #print(f'COGid\tnuma\tnumb\tnume\tnumx\tnumtotal\tnumpattern\tfcat\tnumsplits\ts_root\tavscore_outsideClade\tavscore_insideClade\tnumclades')
    print(outstr)

    

if __name__ == "__main__":
    main()
