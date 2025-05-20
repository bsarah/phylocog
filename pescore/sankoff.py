import argparse
import pandas as pd
from io import StringIO
from skbio import read
from skbio.tree import TreeNode


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--intree", help="tree in newick format")
parser.add_argument("-p", "--inpat", help="list of patterns")
parser.add_argument("-o", "--outpat", help="list of simplified patterns")
parser.add_argument("-t", "--outtree", help="output tree in newick format")
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

outpatternfile = ""
if args.outpat:
    outpatternfile = args.outpat
else:
    print("no output patternfile given!\n")
    exit


outputtree = ""
if args.outtree:
    outputtree = args.outtree
else:
    outputtree = inputtree+"-pescores.newick"



##set path to domain alignment if needed
dnwa_path = "../pairwise-dNWA"

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
        for (xs,ys,z) in v:
            print(f'averaging: xs {xs}, ys {ys}')
            avxs = int(sum(xs)/len(xs))
            avys = int(sum(ys)/len(ys))
            newv.append(f'({avxs},{avys},{z})')
        (a,b) = clus2maxend[clusid]
        newv.append(f'({a},{b},END)')
        print(f'result {newv}')
        print(clusid2label[clusid])
        clus2avranges[clusid] = newv
        outstr = f'avgPattern-P{clusid} (0,0,START)'
        for xs in newv:
            outstr = outstr+" "+xs
        outf.write(f'{clusid2label[clusid]}\n')
        outf.write(f'{outstr}\n')
        
    print(len(clus2avranges.items()))




avPattern(patternfile, outpatternfile)
