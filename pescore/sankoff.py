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



def sankoff(aln2score,inputtree,outputtree):
    #possible labels can be found in the dictionary
    xlab = 'X' #include 'no-sequence-score', this is not allowed at the internal nodes! no propagation of X towards the root
    ks = list(aln2score.keys())
    kks = []
    for (a,b) in ks:
        kks.append(a)
        kks.append(b)

    unique_pids = list(set(kks))
    #print(unique_pids)
    #we need a DP matrix node IDs vs labels
    #postorder node IDs!
    #outfile = open(outputtree,"w")

    pretree = read(inputtree, format="newick", into=TreeNode)
    pretree.bifurcate()
#    print(pretree)
    #print(pretree.ascii_art())
    nodenum = pretree.count()

    #create matrix
    maxval = 100000
    dpm = dict()
    for up in unique_pids:
        for i in range(nodenum):
            dpm[(up,i)] = maxval

    nodeid2bestmatch = dict()
    nodename2id = dict()
    nodename2update = dict()
    curid = 0 #for postorder traversal
    for node in pretree.postorder():
        #print(f'curid {curid}')
        if(node.is_tip()):
            nns = node.name.split('-')
            prelab = nns[-1]
            nodelabel = prelab[1:]
            #print(f'enter leaves: {nodelabel} {curid}')
            dpm[(nodelabel,curid)] = 0
            nodeid2bestmatch[curid] = (nodelabel,0)
            nodename2id[node.name] = curid         
        else:
            #calculate score from children nodes for all possible labels
            #fill dpm[(up,curid)] for current node
            curlab2min = []
            nodename2id[node.name] = curid
            child1 = node.children[0]
            c1name = child1.name
            c1id = nodename2id[c1name]
            #print(f'child1 {c1name} {c1id}')
            if(len(node.children) > 1):
                child2 = node.children[1]
                c2name = child2.name
                c2id = nodename2id[c2name]
                #print(f'child2 {c2name} {c2id}')
            for up in unique_pids:
                if(up == xlab):
                    #skip xlab for internal nodes
                    continue
                #print(f'cur up {up}')
                labnscore1 = [] #store labels with scores for child1
                labnscore2 = [] #store labels with scores for child2
                for vp in unique_pids:
                    #print(f'cur vp {vp}')
                    prevsc1 = dpm[(vp,c1id)]
                    dist1 = aln2score[(up,vp)]
                    #print(f'prevsc1 {prevsc1} dist1 {dist1}')
                    sum1 = prevsc1+dist1
                    labnscore1.append((vp,sum1))
                    if(len(node.children) > 1):
                        prevsc2 = dpm[(vp,c2id)]
                        dist2 = aln2score[(up,vp)]
                        #print(f'prevsc2 {prevsc2} dist2 {dist2}')
                        sum2 = prevsc2+dist2
                        labnscore2.append((vp,sum2))
                #get min from both lists
                labnscore1.sort(key=lambda x:x[1])
                min1 = labnscore1[0][1]
                min2 = 0
                if(len(node.children) > 1):
                    labnscore2.sort(key=lambda x:x[1])
                    min2 = labnscore2[0][1]
                #print(f'min1 {min1} min2 {min2}')
                curscore = round(min1+min2,2)
                dpm[(up,curid)] = curscore
                curlab2min.append((up,curscore))
            curlab2min.sort(key=lambda x:x[1])
            nodeid2bestmatch[curid] = (curlab2min[0][0],curlab2min[0][1])
            nodename2update[node.name] = f'{node.name}-{curlab2min[0][1]}-P{curlab2min[0][0]}'
        curid+=1
    #print(nodeid2bestmatch)

    #recreate the tree with internal node labels corresponding to the best fitting labels
    for node in pretree.postorder():
        if(node.name in nodename2update):
            node.name = nodename2update[node.name]
    #print(pretree)
    #print(pretree.ascii_art())
    pretree.write(outputtree,format='newick')
    return nodeid2bestmatch[nodenum-1] #return results for root
        
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
#print(f'{cogid} {finalscore}')

#calculate split nodes?



