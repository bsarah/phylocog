import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import argparse
import seaborn as sns

#input: input = patternstats file
#output: output folder for pie charts for single COGs
#and barplot file to be plotted at the end showing number of COGs per Xgroup

#format: xgroup\tcog

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile")#
parser.add_argument("-c", "--cog", help="cog ID")
parser.add_argument("-o", "--outfolder") #
parser.add_argument("-b", "--barplot") #name for barplot text file to be plotted at the end
args = parser.parse_args()


inputfile = ""
if args.infile:
    inputfile = args.infile
else:
    print("no inputfile given!\n")
    exit

incog = ""
if args.cog:
    incog = args.cog
else:
    print("no inputcog given!\n")
    exit

outputfolder = ""
if args.outfolder:
    outputfolder = args.outfolder
else:
    print("no outputfolder given!\n")
    exit



outputbarplot = ""
if args.barplot:
    outputbarplot = args.barplot
else:
    print("no output barplot given!\n")
    exit



#plan:
#- keep patternstats file per cog and write another script analysing all the files
#- create pie charts with percentages of xgroups per COG
#- and barplot with xaxis=xgroups and yaxis = counts of COGs, which xgroups occurs in how many cogs


#read patternstats file
pats = pd.read_csv(inputfile,sep='\t')

cog2xlist = dict()
x2coglist = dict()
x2numseq = dict()
for (index,row) in pats.iterrows():
    cogid = row[1]
    numseqs = int(row[2])+int(row[3])
    ps = row[0].split('_')
    for p in ps:
        pps = p.split('.')
        curxgroup = pps[0]
        if curxgroup in x2coglist:
            x2coglist[curxgroup].append(cogid)
            x2numseq[curxgroup] += numseqs
        else:
            x2coglist[curxgroup] = [cogid]
            x2numseq[curxgroup] = numseqs
        if cogid in cog2xlist:
            cog2xlist[cogid].append(curxgroup)
        else:
            cog2xlist[cogid] = [curxgroup]


#open barplottextfile appending text
#transform x2coglist into lines
#xgroup\tcogid
outf = open(outputbarplot,"a")
#prepare data for pie chart
#print(x2coglist)
#plot barplot with bar height=number of same cogid
#problem: no insight into type of structural domain, e.g. L, NC etc
xgs = []
numpats = []
for (k,v) in x2numseq.items():
    if(k=="ENDEND"):
        continue
    xgs.append(str(k))
    numpats.append(v)

    
#df = pd.DataFrame()
#df['xgroup'] = xgs
#df['counts'] = numpats
#dataset = pd.DataFrame({'xgroup': xgs, 'numcogs': numcogs}, columns=['xgroup', 'numcogs'])
#print(df)
#df[:5]

#calculate percentages for xgroups
sumcogs = sum(numpats)
percogs = []

xgs_filtered = []
percogs_filtered = []
for n in range(len(numpats)): #numpats=numpatterns
    #print(f'calculation: {numpats[n]}/{sumcogs}')
    perc = round(numpats[n]/sumcogs,2)
    percogs.append(perc)
    if(perc > 0):
        xgs_filtered.append(xgs[n])
        percogs_filtered.append(perc)
        outf.write(f'{xgs[n]}\t{incog}\t{perc}\n')
        
rem = 1-sum(percogs_filtered)
xgs_filtered.append("rest")
percogs_filtered.append(rem)
#print(percogs)
#print(percogs_filtered)


            
##pie chart:
piefile = f'{outputfolder}/{incog}_patternpie.pdf'

#define Seaborn color palette to use
colors = sns.color_palette('pastel')#[0:5]

#create pie chart
plt.figure(figsize=(10, 10))
plt.pie(percogs_filtered, 
        labels = xgs_filtered, 
        colors = colors, 
        #colors=sns.color_palette('Set2'),
        autopct='%.0f%%', 
        textprops={'fontsize':8})

plt.title(
    label=f'{incog} ECOG X groups', 
    fontdict={"fontsize":10},
    pad=20
)

plt.savefig(piefile)


