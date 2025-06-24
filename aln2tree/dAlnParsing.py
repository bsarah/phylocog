import argparse

#import numpy as np
#import pandas as pd
import math
#import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors

from skbio import DistanceMatrix
from skbio.tree import nj
from skbio.tree import TreeNode
#from skbio.tree import read
import sys

import re

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


#argparse input

parser = argparse.ArgumentParser(description="Parsing pairwise domain alignment output, calculation of pairwise distances and plotting")
parser.add_argument("infile", help = "dNWA output file")
parser.add_argument("-f", "--translationfile", help="optional, file to translate between protein IDs and tax IDs to distinguish between archaea and bacteria")
parser.add_argument("-d", "--distmat", help = "optional, distance matrix output file, if given, the program will calculate the distance matrix")
parser.add_argument("-o","--outfile", help = "output file for modified/simplified dNWA alignment")
#parser.add_argument("-t","--treefiletaxid", help = "output file for newick tree with tax IDs as labels")
parser.add_argument("-p","--treefile", help = "output file for newick tree with IDs as labels")
parser.add_argument('-v','--verbose', action='store_true')

#add some parameters to change?
#alpha/beta for the Feng-Doolittle algorithm?



args = parser.parse_args()

inputfile = ""
if args.infile:
    inputfile = args.infile
else:
    print("no inputfile given!\n")
    exit

distmatfile = ""
if args.distmat:
    distmatfile = args.distmat
else:
    print("no file for distance matrix given!\n")
    exit

transfile = ""
if args.translationfile:
    transfile = args.translationfile
else:
    print("no translation file given, assume already translated IDs!\n")
    exit

outputname = inputfile+".out"
if args.outfile:
    outputname = args.outfile

#in current versions of the previous steps, taxid and protid are merged
treeprotname = inputfile+"prot-tax.newick"
if args.treefile:
    treeprotname = args.treefile

#treetaxname = inputfile+"taxid.newick"
#if args.treefiletaxid:
#    treetaxname = args.treefiletaxid

    

def formatAln(seq1,seq2):
    #go through seq1 and seq2 to get the same number of tuples and matching indices
    i = 0
    j = 0
    res1 = []
    res2 = []
    while(i < len(seq1) and j < len(seq2)):
        (a,b,e) = seq1[i]
        (c,d,f) = seq2[j]
        if(a==c and b==d):
            res1.append((a,b,e))
            res2.append((c,d,f))
            i+=1
            j+=1
        elif(a==c and b!=d):
            #take min of b and d
            if(b<d):
                res1.append((a,b,e))
                res2.append((c,b,f))
                i+=1
            else:
                res1.append((a,d,e))
                res2.append((c,d,f))
                j+=1
        elif(a!=c and b==d):
            #take max of a and c
            if(a>c):
                res1.append((a,b,e))
                res2.append((a,d,f))
            else:
                res1.append((c,b,e))
                res2.append((c,d,f))
            i+=1
            j+=1
        else: #both different
            #take max(a,c) and min(b,d)
            newstart = max(a,c)
            if(b<d):
                res1.append((newstart,b,e))
                res2.append((newstart,b,f))
                i+=1
            else:
                res1.append((newstart,d,e))
                res2.append((newstart,d,f))
                j+=1
    return((res1,res2))

#example
#(0,20), (21,42)
#00000000------
#(0,15), (16,32),(33,42)
#00000AAAAA----BBBBB
#(0,15,0:0),(16,20,0:A),(21,32,-:A),(33,42,-:B) 
#i=0,j=0,    i=0,j=1     i=1,j=1    i=1,j=2


def writeFile(id2aln,outname):
    outf = open(outname,"a")
    for(k,val) in id2aln.items():
        outf.write(f'{str(k)}\n')
        (v,w)= val
        vstr = "(0,0,START)"
        for (a,b,c) in v:
            vstr+=f' ({a},{b},{c})'
        (x,y,z) = v[-1]
        vstr += f' ({y+1},{y+1},END)'
        outf.write(f'{vstr}\n')
        wstr = "(0,0,START)"
        for (a,b,c) in w:
            wstr+=f' ({a},{b},{c})'
        (x1,y1,z1) = w[-1]
        wstr += f' ({y1+1},{y1+1},END)'
        outf.write(f'{wstr}\n')

        
#read and parse input file
def readInput(filename):
    id2aln = dict()
    id2selfaln = dict()
    with open(filename) as f:
        curid = ""
        seq1 = []
        seq2 = []
        for pline in f:
            if(pline[0] == ">"):
                line = pline.strip()
                lc = line.split(',')
                if(len(lc) < 3):
                    continue
                #check for self alignment
                if(lc[1] == lc[3]):
                    id2selfaln[lc[1]] = (lc[4],lc[5])
                if(curid != ""):
                    #add current entry to dict
                    if(curid not in id2aln):
                        #go through seq1 and seq2 to get the same number of tuples and matching indices
                        id2aln[curid] = formatAln(seq1,seq2)
                curid = line
                seq1 = []
                seq2 = []
            else:
                line = pline.strip()
                lls = line.split(',')
                tmpl = []
                curelem = lls[0]
                curstart = 0
                for i in range(1,len(lls)):
                    if(lls[i] != curelem):
                        newpair = (curstart,i-1, curelem)
                        tmpl.append(newpair)
                        curelem = lls[i]
                        curstart = i
                newpair = (curstart,len(lls)-1, curelem)
                tmpl.append(newpair)
                if(seq1 == []):
                    seq1 = tmpl
                else:
                    seq2 = tmpl
        if(curid not in id2aln):
            #go through seq1 and seq2 to get the same number of tuples and matching indices
            id2aln[curid] = formatAln(seq1,seq2)
    return((id2aln,id2selfaln))


def scoringFun(i,j):
    """Scoring function, assign the values for match calculation and return the value.
    This is the scoring function as given in the dNWA algorithm.
    TODO: externalize such that it can be changed more easily.
    """

    # variables
    score = 0

    #scoring function
    # GAP values
    PENALTY_GAP_OPENING = -1
    PENALTY_GAP_EXTENSION = -1

    # Scoring Matrix
    WEIGHT_MATCH_F_GROUP = 8
    WEIGHT_MATCH_T_GROUP = 4
    WEIGHT_MATCH_H_GROUP = 2
    WEIGHT_MATCH_X_GROUP = 1
    WEIGHT_MATCH_NO_FID = 0
    WEIGHT_MISMATCH_FIDS = -5
    WEIGHT_MISMATCH_NOFID_FID = -1

    # split label to compare the different groups
    list_with_splitted_groups_i = i.split('.')
    list_with_splitted_groups_j = j.split('.')
    if len(list_with_splitted_groups_i)<4:
        list_with_splitted_groups_i.append('None')
    if len(list_with_splitted_groups_j)<4:
        list_with_splitted_groups_j.append('None')

    # scoring logic
    # 0 indicates "no labeled proteindomain" 
    if i == '0' and j == '0':
        score = WEIGHT_MATCH_NO_FID
    # match
    elif i.strip() == j.strip():
        score = WEIGHT_MATCH_F_GROUP
    # one sequence with NOFID, but not both
    elif i == '0' or j == '0':
        score = WEIGHT_MISMATCH_NOFID_FID
    #compare X-group, then H-group, ...
    elif list_with_splitted_groups_i[0] == list_with_splitted_groups_j[0]:
        if list_with_splitted_groups_i[1] == list_with_splitted_groups_j[1]:
            if list_with_splitted_groups_i[2] == list_with_splitted_groups_j[2]:
                if list_with_splitted_groups_i[3] == list_with_splitted_groups_j[3] and list_with_splitted_groups_i[3] != None:
                    score = WEIGHT_MATCH_F_GROUP
                else:
                    score = WEIGHT_MATCH_T_GROUP
            else:
                score = WEIGHT_MATCH_H_GROUP
        else:
            score = WEIGHT_MATCH_X_GROUP
    else:
        score = WEIGHT_MISMATCH_FIDS

    return score


def calcDistMat(id2aln, id2selfaln, distmatfile, treeprotname):
    #create distance matrix for pairs of IDs in id2selfaln
    #calculation based on Feng-Doolittle algorithm
    #https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Feng-Doolittle

    #sequences a,b, pairwise score S(a,b)
    #distance: D(a,b) = -ln (Seff(a,b)*f)  --with f =1
    #Seff(a,b) = (S(a,b)-Srand(a,b))/(Smax(a,b)-Srand(a,b))
    #Smax(a,b) = (S(a,a)+S(b,b))/2 --self alignments
    #Srand(a,b) = 1/len(aln(a,b)) * (sum_x sum_y s(x,y) * Na(x)*Nb(y)) + beta*Ne + alpha *No
    #Na(x), Nb(y) number of characters x,y in sequences a,b
    #Ne,No number of gap openings and gap extensions
    #alpha/beta: parameters from Gotoh (-3,-1) or from our alignment algorithm, we only use -1, no gap extensions, thus use -(number of gaps)
    print("calculating distances")
    if(distmatfile):
        outd = open(distmatfile,"a")
    idlist = list(id2selfaln.keys())
    numids = len(idlist)
    data = [[0 for x in range(numids)] for y in range(numids)] 
    id2dist = dict()
    for (k,v) in id2aln.items():
        ks = k.split(',') #[>clus1,id1,>clus2,id2,score,len]
        (sa,sb) = v #two lists of triples (start,end,character)
        saa = 0
        if(ks[1] in id2selfaln):
            (saa,lensa) = id2selfaln[ks[1]]
        else:
            print("no self aln for seq A!")
        sbb = 0
        if(ks[3] in id2selfaln):
            (sbb,lensb) = id2selfaln[ks[3]]
        else:
            print("no self aln for seq B!")
        Smax = (float(saa)+float(sbb))/2
        Sab = float(ks[4])
        lenab = int(ks[5])
        numg = 0
        gapsa = 0
        gapsb = 0
        sacounts = dict()
        #sacountsum = 0
        for (x,y,z) in sa:
            if(z == "START" or z == "END"):# or z =='0' or z =='-'):
                continue
            if(z in sacounts):
                sacounts[z] += (y-x)
            else:
                sacounts[z] = (y-x)
            if(z == "-"):
                gapsa += (y-x)
        sbcounts = dict()
        #sbcountsum = 0
        for (x,y,z) in sb:
            if(z == "START" or z == "END"):# or z == '0' or z =='-'):
                continue
            if(z in sbcounts):
                sbcounts[z] += (y-x)
            else:
                sbcounts[z] = (y-x)
            if(z == "-"):
                gapsb += (y-x)
        #calculate this sum:(sum_x sum_y s(x,y) * Na(x)*Nb(y))
        totsum = 0

        for (ka,va) in sacounts.items():
            for (kb,vb) in sbcounts.items():
                totsum+= scoringFun(ka,kb) * va * vb
        alphabeta = -1 #this should be two parameters of the feng-doolittle algorithm, for gap opening/gap extension, in the freiburg teaching page, they are set to -2 or -3: https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Feng-Doolittle
        Srand = 1/lenab * totsum + (alphabeta*(gapsa+gapsb)) 
        #distance: D(a,b) = -ln (Seff(a,b)*f)  --with f =1
        #Seff(a,b) = (S(a,b)-Srand(a,b))/(Smax(a,b)-Srand(a,b))
        #what do we do if Sab < 0?
        if(Sab<0 and args.verbose):
            print(f'Sab {Sab} < 0 for k:{k}')
            Sab = 0.00001
        Seff = (Sab-Srand)/(Smax-Srand)
        #print(f'{Seff}')
        if(Seff < 0 and args.verbose):
            print(f'neg seff {Seff} for {sa} and {sb} ({k}): srand: {Srand}, lenab:{lenab}, smax: {Smax}, sab: {Sab}, saa: {saa}, sbb: {sbb}, lensa:{lensa}, lensb:{lensb}, gapsab: {gapsa},{gapsb}, totsum:{totsum}')
        if(Seff == 0):
            Seff = 0.001
        Sdist = -1 * math.log(Seff)
        if(Sdist < 0 and args.verbose):
            print(f'Distance < 0: {Sdist} for k:{k}')
#        print(f'eff score and distance for {ks[1]} and {ks[3]}: {Seff}, {Sdist}')
        #test tree
        kid1 = idlist.index(ks[1])
        kid2 = idlist.index(ks[3])
        data[kid1][kid2] = Sdist
        data[kid2][kid1] = Sdist
        if(distmatfile):
            outstr1 = f'{ks[1]}\t{ks[3]}\t{Sdist}\n'
            outstr2 = f'{ks[3]}\t{ks[1]}\t{Sdist}\n'
            outd.write(outstr1)
            outd.write(outstr2)


    if(len(data) > 2):
        #skip creating the cluster tree if there are not enough nodes,
        #this occurs when creating the averaged patterns and the tree of averaged patterns
        #print(f'ERROR create distance matrix for {distmatfile}', file=sys.stderr)

        dm = DistanceMatrix(data, idlist)
        tree = nj(dm)
        curnodeid = 0
        #go through tree and add internal node names
        for node in tree.postorder():
            if(node.is_root()):
                tmpname = "r"+str(curnodeid)
                node.name = tmpname
            elif(node.has_children()):
                tmpname = "i"+str(curnodeid)+"i"
                node.name = tmpname
            else:
                tmpname = ""
                #we don't want to overwrite the leave names!
            curnodeid+=1
        #print(tree.ascii_art())
        tree.write(treeprotname)
        #newick_str = nj(dm, result_constructor=str)
        #    print(newick_str[:55], "...")
        #return newick_str
        #print(newick_str)
        #in a different script?
        #use neighbor joining with midpoint rooting to construct the tree?
        #http://scikit-bio.org/docs/0.2.1/generated/skbio.tree.nj.html


#check for duplications
#-calculate score at the same time?
#-sort input and find self-alignments?
#-create distance matrix?

def reformat(protid2geneid):
    print(f'read from file {treeprotname}')
    pretree = TreeNode.read(treeprotname, format="newick")
    pretree.bifurcate()

    id = 0
    redcol = "#ff0000"
    bluecol = "#0000ff"
    blackcol = "#000000"
    lab = "label"
    fontspec = "normal"
    fontsize = "1.5"
    numarc = 0
    numbac = 0
    numx = 0
    
    for node in pretree.preorder():
        if(node.is_root()):
            tmpname = "r"+str(id)
        elif(node.has_children()):
            #internal node
            #node.name = support
            #node.length = length/distance to parent
            #print("internal node")
            #print(node.name)
            #print(node)
            #print(node.length)
            tmpname = "i"+str(id)+"i"
            node.support = 0.0
            if(node.name != None):
                node.support = float(node.name)
        else:
            protid = node.name
            protid = re.sub(r"\s+", '_', protid)
            #print(protid)
            #node.length = length/distance to parent
            #node.name = identifier, protein id
            #exchange protein ID by taxid
            taxid = "xxxx"
            curcol = blackcol
            curdat = ["x"]
            if(protid in protid2geneid):
                taxid = protid2geneid[protid]
                #curdat = curline.split("\t")
                #taxid = curdat[0]
                curcol = blackcol
                if(taxid[0] == 'a'):
                    numarc = numarc+1
                    curcol = redcol
                elif(taxid[0] == 'b'):
                    numbac = numbac+1
                    curcol = bluecol
            else:
                numx+=1
            #leaves
            tmpname = str(id)+'-'+taxid
            curdat[0] = tmpname
            #outline = "\t".join(curdat)
            #tabfile.write(outline)
            #colfile.write(tmpname+"\t"+lab+"\t"+curcol+"\t"+fontspec+"\t"+fontsize+"\n")
            #print(node.name)
            #print(tmpname)
            #protid = translation[node.name]
            #print(protid)
        node.name = tmpname
        id=id+1
    
    #print("Number of archaea", numarc)
    #print("Number of bacteria", numbac)
    #print("Number of unknown", numx)
    #print(pretree.ascii_art())

    tree = pretree
    tree2 = str(pretree)
    #print(tree2)
    ##go again through the string and change format

    newstr = "" #this will be the new tree string with names replaced
    curword = "" #this will be the new complete name
    cursupport = "" #this will be the support value that is added to the name at the end
    dowrite = 0 #if 1, add chars to curword
    doadd = 0 #if 1 add supp
    for c in tree2:
        if(c == '\'' and dowrite == 0):
            dowrite = 1
        elif(c == '\'' and dowrite == 1):
            dowrite = 0
            doadd = 1
            #change name and add support value
            vals = curword.split(':')
            cursupport = vals[0]
            curword = vals[1]
        elif(dowrite == 1):
            curword = curword+c
        elif(doadd == 1 and (c == ',' or c == ')' or c == '(')):
            #stop adding to curword
            doadd = 0
            newstr = newstr + curword + '[' + cursupport + ']' + c
            curword = ""
            cursupport = ""
        elif(doadd==1):
            curword = curword+c
        else:
            newstr = newstr + c

    #print(newstr) #output tree
    #tabfile.close()
    #colfile.close()
    
    #use own string to write to file such that iTOL can read it
    treefile = open(treetaxname,"w")
    treefile.write(newstr)
    treefile.close()




def main():
    (id2aln,id2selfaln) = readInput(inputfile)
    writeFile(id2aln,outputname)

    #translation file
    protID2geneid = dict()
    protID2line = dict()

    if(args.translationfile):
        file1 = open(transfile, 'r')
        inlines = file1.readlines()
    
        for line in inlines:
            line2 = line.rstrip()
            cols = line2.split("\t")
            protID2geneid[cols[3]] = cols[0]
            protID2line[cols[3]] = line2
        
        file1.close()


    calcDistMat(id2aln, id2selfaln,distmatfile,treeprotname)

    #reformatting is not necessary here
    #reformat(protID2geneid)
    

main()
