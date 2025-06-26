import argparse
import pandas as pd
from io import StringIO
from skbio import read
from skbio.tree import TreeNode
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--intree", help="tree in newick format")
parser.add_argument("-p", "--inpat", help="list of patterns")
parser.add_argument("-c", "--cogid", help="COGid for the output")
parser.add_argument("-f", "--fcat", help="functional category")
parser.add_argument("-o", "--outfolder", help="folder for output files")
#parser.add_argument("-o", "--outpat", help="list of simplified patterns")
#parser.add_argument("-t", "--outtree", help="output tree in newick format")
args = parser.parse_args()

inputtree = ""
if args.intree:
    inputtree = args.intree
else:
    print("no inputtree given!\n")
    exit

patternfile = ""
if args.inpat:
    patternfile = args.inpat
else:
    print("no patternfile given!\n")
    exit

cogid = "COGXXXX"
fcat = "XX"

if(args.cogid):
    cogid = args.cogid

if(args.fcat):
    fcat = args.fcat


outputfolder = "."
if(args.outfolder):
    outputfolder = args.outfolder
#outpatternfile = ""
#if args.outpat:
#    outpatternfile = args.outpat
#else:
#    print("no output patternfile given!\n")
#    exit


#outputtree = ""
#if args.outtree:
#    outputtree = args.outtree
#else:
#    outputtree = inputtree+"-pescores.newick"



##set path to domain alignment if needed
dnwa_path = "/home/sarah/projects/cogupdate/phylocog/pairwise-dNWA"
alnparse_path = "/home/sarah/projects/cogupdate/phylocog/aln2tree"

##calculate the average pattern in case a pattern has more than 1 sequence
def avPattern(patternfile, outpatternfile):
    outf = open(outpatternfile,"a")
    clusid2label = dict()
    clus2ranges = dict() #clusid to [([a],[b],domid)] -> lists to calculate average later
    clus2maxend = dict() #clusid to max end cordinate
    curclusid = -1
    curacc = 0
    curann = ""
    curcolid = 0 #to set unique colors
    with open(patternfile) as f:
        for line in f:
            curline = line.strip()
            if(curline[0] == '>'): #id line
                cs = curline.split(';')
                #clusnum
                cs0 = cs[0].split(" ")
                curclusid = int(cs0[-1])
#                print(f'current pattern: {curclusid}')
                #accessions
                cs1 = cs[1].split(" ")
                curacc = cs1[-1]
                #annotation
                cs2 = cs[2].split(" ")
                curann = cs2[-1]
                clusid2label[curclusid] = curline
            else: #line describing one sequence
                #>avgPattern-P0 (0,0,START) (92,160,109.54.1.1) (199,317,266.1.1.1) (330,340,END)
                dseq = curline.split(" ")[2:] #this is the domain sequence including end but without start
                #print(dseq)
                prestr = dseq[-1][1:-1] #END tag
                #print(prestr)
                ps = prestr.split(',')
                a = int(ps[0])
                b = int(ps[1])
                ef = ps[2]
                if curclusid in clus2maxend:
                    (e,f) = clus2maxend[curclusid]
                    if(b > f):
                        clus2maxend[curclusid] = (a,b)
                    else:
                        clus2maxend[curclusid] = (a,b)
                else:
                    clus2maxend[curclusid] = (a,b)
                #print(len(clus2maxend.items()))
                #remaining sequence
                if curclusid in clus2ranges:
                    dictlist = clus2ranges[curclusid]
                    for j in range(len(dictlist)):
                        (xs,ys,z) = dictlist[j]
                        #(ac,bc,cc) = dseq[j]
                        predseq = dseq[j][1:-1]
                        psd = predseq.split(',')
                        ac = int(psd[0])
                        bc = int(psd[1])
                        cc = psd[2]
                        #z and cc should be the same
                        if(cc != z):
                            print(f'WARNING {cc} unequals {z}!\n')
                        xs.append(ac)
                        ys.append(bc)

                else:
                    newdlist = []
                    for k in range(len(dseq)-1):
                        predseq = dseq[k][1:-1]
                        psd = predseq.split(',')
                        ac = int(psd[0])
                        bc = int(psd[1])
                        cc = psd[2]
                        newdlist.append(([ac],[bc],cc))
 #                   print(f'append {newdlist}')
                    clus2ranges[curclusid] = newdlist




    #analyse input before creating the plot
    #x-limits [0, max END]
    #y-limits [number of clusters * 10]
    #x-ticks
    #y-ticks
    #different colors
    #length and height of the bars

    ends_sorted = sorted(clus2maxend.items(),key=lambda x: x[1])
    maxlim = ends_sorted[-1][1][0]
    #print(maxlim)
    numclus = len(ends_sorted)
    #print(numclus)

    clus2avranges = dict() #take the average in case of 
    for (clusid,v) in clus2ranges.items():
        newv = []
        #use realv in order to check that the next start coordinate is larger than current end coordinate, otherwise, dNWA will throw an error
        realv = []
        realv.append((0,0,"START"))
        for (xs,ys,z) in v:
            #print(f'averaging: xs {xs}, ys {ys}')
            avxs = int(sum(xs)/len(xs))
            avys = int(sum(ys)/len(ys))
            (pa,pb,pc) = realv[-1]
            if (avys < pa):
                avys = pa+1
            newv.append(f'({avxs},{avys},{z})')
            realv.append((avxs,avys,z))
        #print(clusid)
        #print(realv)
        (a,b) = clus2maxend[clusid]
        #print(f'END {a} {b}')
        (pa,pb,pc) = realv[-1]
        if(a<pb):
            a=pb+1
        newv.append(f'({a},{b},END)')
        #print(f'result {newv}')
        #print(clusid2label[clusid])
        clus2avranges[clusid] = newv
        outstr = f'avgPattern-P{clusid} (0,0,START)'
        for xs in newv:
            outstr = outstr+" "+xs
        outf.write(f'{clusid2label[clusid]}\n')
        outf.write(f'{outstr}\n')
        
    #print(len(clus2avranges.items()))


def read_dNWA_aln(alnfile):
    aln2score = dict() #aln between two patterns -> normalized score
    xlab = 'X' #include 'no-sequence-score'
    maxval = 100000
    with open(alnfile) as f:
        for line in f:
            curline = line.strip()
            ls = curline.split('\t')
            preid1 = ls[0]
            ps1 = preid1.split('-')
            id1 = ps1[-1][1:]
            preid2 = ls[1]
            ps2 = preid2.split('-')
            id2 = ps2[-1][1:]
            nscore = round(float(ls[2]),3)
            alnid1 = (id1,id2)
            alnid2 = (id2,id1)
            aln2score[alnid1] = nscore
            aln2score[alnid2] = nscore
            alnidx2 = (xlab,id2)
            alnid2x = (id2,xlab)
            alnidx1 = (xlab,id1)
            alnid1x = (id1,xlab)
            aln2score[alnid1x] = 0
            aln2score[alnidx1] = 0
            aln2score[alnid2x] = 0
            aln2score[alnidx2] = 0
    #print(aln2score)
    alnidxx = (xlab,xlab)
    aln2score[alnidxx] = maxval
    return aln2score



#def calculateChildOptimum(parentlabel, aln2score, unique_pids, dpm):
    



def sankoff(aln2score,inputtree,outputtree):
    #possible labels can be found in the dictionary
    xlab = 'X' #include 'no-sequence-score', this is not allowed at the internal nodes! no propagation of X towards the root
    ks = list(aln2score.keys())
    kks = []
    for (a,b) in ks:
        kks.append(f'{a}')
        kks.append(f'{b}')

    unique_pids = list(set(kks))
    #print(unique_pids)
    #we need a DP matrix node IDs vs labels
    #postorder node IDs!
    #outfile = open(outputtree,"w")
    nodenames = []

    #use idx function to get id of pid or node id
    
    pretree = read(inputtree, format="newick", into=TreeNode)
    pretree.bifurcate()
#    print(pretree)
    #print(pretree.ascii_art())
    nodenum = pretree.count()
    curid = 0 #for postorder traversal
    nodename2id = dict()
    for node in pretree.postorder():
#        if(node.is_tip()):
#            nns = node.name.split('-')
#            nodelabel = nns[-1]
        nodenames.append(node.name)
        nodename2id[node.name] = curid
        curid+=1
        
    #create matrix and init    
    maxval = 100000
    w, h = len(unique_pids), len(nodenames) #h = rows
    #access matrix with Matrix[h][w]
    Matrix = [[maxval for x in range(w)] for y in range(h)] 

    #forward step
    for node in pretree.postorder():
        nid = nodenames.index(node.name)
        if(node.is_tip()):
            #print(f'enter leaves: {nodelabel} {curid}')
            nns = node.name.split('-')
            nodelabel = nns[-1][1:]
            pid = unique_pids.index(nodelabel)
            Matrix[nid][pid] = 0
        else:
            #calculate score from children nodes for all possible labels
            #fill dpm[(up,curid)] for current node
            child1 = node.children[0]
            c1name = child1.name
            c1id = nodenames.index(c1name)
            #print(f'child1 {c1name} {c1id}')
            if(len(node.children) > 1):
                child2 = node.children[1]
                c2name = child2.name
                c2id = nodenames.index(c2name)
                #print(f'child2 {c2name} {c2id}')
            for up in unique_pids:
                #assume current node has label up, calculate possible scores for children and take mininmum
                uid = unique_pids.index(up)
                if(up == xlab):
                    #skip xlab for internal nodes
                    continue
                labnscore1 = [] #store labels with scores for child1
                labnscore2 = [] #store labels with scores for child2
                for vp in unique_pids:
                    #print(vp)
                    if(vp == xlab):
                        #skip xlab for internal nodes
                        continue
                    vid = unique_pids.index(vp)
                    c1score = Matrix[c1id][vid]+aln2score[(up,vp)]
                    labnscore1.append((vp,c1score))
                    if(len(node.children) > 1):
                        c2score = Matrix[c2id][vid]+aln2score[(up,vp)]
                        labnscore2.append((vp,c2score))
                #get min from both lists
                labnscore1.sort(key=lambda x:x[1])
                min1 = labnscore1[0][1]
                min2 = 0
                if(len(node.children) > 1):
                    labnscore2.sort(key=lambda x:x[1])
                    min2 = labnscore2[0][1]
                #print(f'up {up} at {curid} min1 {min1} min2 {min2}')
                curscore = round(min1+min2,2)
                Matrix[nid][uid] = curscore
    #DP matrix dpm is filled, do backtracking to find the optimal solution for the complete tree
#    print("forward step done\n")
#    upstr = ' '.join(unique_pids)
#    print(f'{upstr}\n')
#    for v in range(nodenum):
#        linestr = f'{v} '
#        for pp in range(len(unique_pids)):
#            linestr+=f'{Matrix[v][pp]} '
#        print(f'{linestr}\n')
    ####BACKTRACK
    nodename2bestpid = dict() #nodename -> pid
    nodename2update = dict()
    rootscore = 0
    rootlabel = ""
    for node in pretree.preorder():
        nid = nodenames.index(node.name)
        if(node.is_root()):
            minpid = len(unique_pids)+2
            minscore = maxval
            for t in range(len(Matrix[nid])):
                cursc = Matrix[nid][t]
                if(cursc < minscore):
                    minpid = t
                    minscore = cursc 
            #mindpid is the label chosen for the current node
            nodename2bestpid[node.name] = minpid
            nodename2update[node.name] = f'{node.name}-{minscore}-P{unique_pids[minpid]}'
            rootscore = minscore
            rootlabel = unique_pids[minpid]
        else:
            parentname = node.parent.name
            pnid = nodenames.index(parentname)
            ppid = nodename2bestpid[parentname]
            pp = unique_pids[ppid]
            #calculate the own best score + label based on the parent label set
            labnscore1 = [] #store labels with scores for this node
            for vp in unique_pids:
                if(vp == xlab):
                    #skip xlab for internal nodes
                    continue
                vid = unique_pids.index(vp)
                c1score = Matrix[nid][vid]+aln2score[(pp,vp)]
                labnscore1.append((vp,c1score))
            #get min from both lists
            labnscore1.sort(key=lambda x:x[1])
            min1 = labnscore1[0][1]
            minp = labnscore1[0][0]
            if(not node.is_tip()):
                #only for inner nodes!
                #mindpid is the label chosen for the current node
                nodename2bestpid[node.name] = unique_pids.index(minp)
                nodename2update[node.name] = f'{node.name}-{min1}-P{minp}'
                
    #    #recreate the tree with internal node labels corresponding to the best fitting labels
    for node in pretree.postorder():
        if(node.name in nodename2update):
            node.name = nodename2update[node.name]
    #print(pretree)
    #print(pretree.ascii_art())
    pretree.write(outputtree,format='newick')
    return (rootscore,rootlabel)
        
#calculate average patterns
#cwd = os.getcwd()
base2=os.path.basename(patternfile)
outpatternfile = f'{outputfolder}/averaged_{base2}'
outpatname = f'averaged_{base2}'

print(f'create {outpatternfile}')
avPattern(patternfile, outpatternfile)

#run dNWA
#python ../pairwise-dNWA/dNWA.py outpattern3.txt
print(f'create alignments_{outpatternfile}')
cmd1 = f'python {dnwa_path}/dNWA.py -o {outputfolder} {outpatternfile}'
subprocess.call(cmd1,shell=True)

#outputfile will be called: alignments_outpattern3.txt
alnfile = f'{outputfolder}/alignments_{outpatname}'

alnfile_format = f'{outputfolder}/alignments_formatted_{outpatname}'

alncluster = f'{outputfolder}/alignments_formatted_{outpatname}.newick'

alndistmat = f'{outputfolder}/distancematrix_formatted_{outpatname}.tsv'
#run parsing of output
#parse dNWA output, get formatted version and cluster tree for av patterns
#python /home/sarah/projects/cogupdate/phylocog/aln2tree/dAlnParsing.py -v -o alignments_outpattern3_formatted.txt -p outpattern3_formatted.newick alignments_outpattern3.txt
#format: >pid1,seqid1,>pid2,seqid2,totalscore,alnlength
#>0,avgPattern-P0,>1,avgPattern-P1,1951,326
print(f'create {alnfile_format} {alncluster} {alndistmat}')
cmd2 = f'python {alnparse_path}/dAlnParsing.py -o {alnfile_format} -p {alncluster} -d {alndistmat} {alnfile}'#we need the distance matrix as Neya's program calculates similarity scores!
subprocess.call(cmd2,shell=True)

#read formatted alignment file and store in a dictionary, normalize scores!
aln2score = read_dNWA_aln(alndistmat) #this includes self alignments

base=os.path.basename(inputtree)
#its = inputtree.split('/')
#reipt = its[-1] #real input tree without outfolder

outputtree = f'{outputfolder}/sankoff_{base}'
print(f'create {outputtree}')

#call sankoff algorithm to calculate score for 
finalscore = sankoff(aln2score,inputtree,outputtree)
print(f'{cogid} {finalscore}')

#calculate split nodes?



