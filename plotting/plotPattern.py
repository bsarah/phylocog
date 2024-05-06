import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import argparse

#input: clustered structure domain sequences (no alignment!)
#this plot is only to visualize the structural domain composition of sequences
#inside a COG

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile")#clustered structural domain sequences
parser.add_argument("-o", "--outfile") #name for plot
args = parser.parse_args()


inputfile = ""
if args.infile:
    inputfile = args.infile
else:
    print("no inputfile given!\n")
    exit


outputfile = ""
if args.outfile:
    outputfile = args.outfile
else:
    print("no outputfile given!\n")
    exit




colordict = mcolors.CSS4_COLORS
colordict.pop('black', None)
colordict.pop('white', None)
colordict.pop('whitesmoke', None)
colordict.pop('seashell', None)
colordict.pop('mistyrose', None)
colordict.pop('linen', None)
colordict.pop('cornsilk', None)
colordict.pop('ivory', None)
colordict.pop('beige', None)
colordict.pop('honeydew', None)
colordict.pop('lightcyan', None)
colordict.pop('lavenderblush', None)
colordict.pop('lemonchiffon', None)
colordict.pop('ghostwhite', None)
colordict.pop('floralwhite', None)
colordict.pop('oldlace', None)
colordict.pop('azure', None)
colordict.pop('aliceblue', None)
colordict.pop('mintcream', None)
colordict.pop('snow', None)
colordict.pop('antiquewhite', None)
colordict.pop('lightyellow', None)
colordict.pop('lightgoldenrodyellow', None)

colordict.pop('aqua', None)
colordict.pop('aquamarine', None)
colordict.pop('bisque', None)
colordict.pop('blanchedalmond', None)
    
colors = list(colordict.keys())




clusid2label = dict()
#clusdom2ranges = dict() #(clusnum,domID) to list of ranges of this domain in the cluster
clus2ranges = dict() #clusid to [([a],[b],domid)] -> lists to calculate average later
clus2maxend = dict() #clusid to max end cordinate
domid2col = dict()
domid2col["END"] = "gray"
curclusid = -1
curacc = 0
curann = ""
curcolid = 0 #to set unique colors
with open(inputfile) as f:
    for line in f:
        curline = line.strip()
        if(curline[0] == '>'): #id line
            cs = curline.split(';')
            #clusnum
            cs0 = cs[0].split(" ")
            curclusid = int(cs0[-1])
            #accessions
            cs1 = cs[1].split(" ")
            curacc = cs1[-1]
            #annotation
            cs2 = cs[2].split(" ")
            curann = cs2[-1]
            clusid2label[curclusid] = f'#{curclusid}:{curacc}x:{curann}'
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
                    if cc not in domid2col:
                        domid2col[cc] = colors[curcolid]
                        curcolid+=1
            else:
                newdlist = []
                for k in range(len(dseq)-1):
                    #(ac,bc,cc) = dseq[k]
                    predseq = dseq[k][1:-1]
                    psd = predseq.split(',')
                    ac = int(psd[0])
                    bc = int(psd[1])
                    cc = psd[2]
                    newdlist.append(([ac],[bc],cc))
                
                    if cc not in domid2col:
                        domid2col[cc] = colors[curcolid]
                        curcolid+=1          
                clus2ranges[curclusid] = newdlist

#print(domid2col)




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
        avxs = int(sum(xs)/len(xs))
        avys = int(sum(ys)/len(ys))
        newv.append((avxs,avys,z))
    clus2avranges[clusid] = newv

#print(len(clus2avranges.items()))


# Declaring a figure "gnt"
fig, gnt = plt.subplots()
 
# Setting Y-axis limits
gnt.set_ylim(0, numclus*10+5)
 
# Setting X-axis limits
gnt.set_xlim(0, maxlim+30)
 
# Setting labels for x-axis and y-axis
gnt.set_xlabel('Sequence length')
gnt.set_ylabel('Pattern')

ytickpos = []
yticklabs = []
for i in range (0,numclus):
    curpos = i*10+5
    ytickpos.append(curpos)
    yticklabs.append(clusid2label[i])

#print(ytickpos)
#print(yticklabs)

# Setting ticks on y-axis
gnt.set_yticks(ytickpos)
# Labelling tickes of y-axis
gnt.set_yticklabels(yticklabs, fontsize=4)

# Setting graph attribute
gnt.grid(True)
gnt.set_axisbelow(True)
#gnt.yaxis.grid(color='gray', linestyle='dashed')

lineoffset = 0

#plot END markers
curcol = domid2col["END"]

for (k,v) in clus2maxend.items():
    cury = k*10+lineoffset
    curx = v[0]
    #gnt.broken_barh([(xcoord, length fo bar)], (ycoord, height of bar), facecolors =color)
    gnt.broken_barh([(curx, 20)], (cury+1, 8), facecolors =curcol)
    gnt.text(s= "END", x=v[0]+6, y=cury+4.5,color='white', weight='bold',fontsize=4)

#plot 
for (k,v) in clus2avranges.items():
    cury = k*10+lineoffset
    offset = 3
    count = 0
    #sort v
    sortedV = sorted(v , key=lambda x: x[0])
    for (a,b,c) in sortedV:
        count+=1
#        print(f'{a},{b},{c}')
        if c in domid2col:
            curcol = domid2col[str(c)]
        else:
            print(f'{c} not in domid2col')
        #print(curcol)
        curx = a
        curlen = b-a
        gnt.broken_barh([(curx, curlen)], (cury+1, 8), facecolors =curcol)
        if(count % 2 == 0):
            offset = 6
        else:
            offset = 3
        gnt.text(s= c, x=curx+1, y=cury+offset,color='black', weight='bold',fontsize=4)

#plt.show()
plt.savefig(outputfile,bbox_inches = "tight")
