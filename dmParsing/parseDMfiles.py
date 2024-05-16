import argparse
import subprocess

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


parser = argparse.ArgumentParser(description="Cluster domain sequences and plot one representative per cluster for all COGs in the input list")
parser.add_argument("-i", "--infile", help = "list of DM output files")
parser.add_argument("-f" ,"--infolder", help = "location of DM output files")
parser.add_argument("-s", "--seqfile", help = "optional, list of fasta files with corresponding protein sequences")
parser.add_argument("-g" ,"--fastafolder", help = "optional, location of fasta files")
parser.add_argument("-t", "--transfile", help = "optional, translation file to translate between protein IDs and tax IDs as in the trees.")
parser.add_argument("-o","--outname", help = "name for produced files as output") #we can have the plot as well as the translation between clusters and contained sequences/accessions
parser.add_argument("-p","--outfolder", help = "location for output files") #we can have the plot as well as the translation between clusters and contained sequences/accessions

args = parser.parse_args()

inputfile = ""
if args.infile:
    inputfile = args.infile
else:
    print("no inputfile given!\n")
    exit

inputfolder = ""
if args.infolder:
    inputfolder = args.infolder
else:
    print("no inputfolder given!\n")
    exit


    
doseqlen = False
sequencefile = ""
if args.seqfile:
    sequencefile = args.seqfile
    doseqlen= True

fastafolder = ""
if args.fastafolder:
    fastafolder = args.fastafolder

dotranslate = False
translationfile = ""
if args.transfile:
    translationfile = args.transfile
    dotranslate = True


outputname = inputfile+".out"
if args.outname:
    outputname = args.outname
outputfolder = "."
if args.outfolder:
    outputfolder = args.outfolder





cog2fasta = dict()
if(doseqlen):
    with open(sequencefile) as f:
        for line in f:
            curline = line.strip()
            ls = curline.split(".")
            curcog = ls[0]
            curfasta = f'{fastafolder}/{curline}'
            cog2fasta[curcog] = curfasta


print("COG\tnProteins\tnDomains\tnL(%)\tnNC(%)\tnCP(%)\tnIS(%)\tnClus\tnSDoms\tnSeqs")


with open(inputfile) as l:
    for line in l:
        curline = line.strip()
        ls = curline.split("_")
        curcog = ls[0]
        #print(f'{curcog}\n')
        curfile = f'{inputfolder}/{curline}'
        cmd1 = ""
        outfile = f'{outputfolder}/{curcog}_{outputname}'
        if(doseqlen and dotranslate):
            fafile = ""
            if(curcog in cog2fasta):
                fafile = cog2fasta[curcog]
            else:
                print(f'no fasta file for {curcog}!')

            cmd1 = f'python3 /home/sarah/projects/cogupdate/phylocog/dmParsing/groupPatterns.py {curfile} -s {fafile} -o {outfile} -t {translationfile}'
        elif(doseqlen):
            cmd1 = f'python3 /home/sarah/projects/cogupdate/phylocog/dmParsing/groupPatterns.py {curfile} -s {fafile} -o {outfile}'
        elif(dotranslate):
            cmd1 = f'python3 /home/sarah/projects/cogupdate/phylocog/dmParsing/groupPatterns.py {curfile} -o {outfile} -t {translationfile}'
        else:
            cmd1 = f'python3 /home/sarah/projects/cogupdate/phylocog/dmParsing/groupPatterns.py {curfile} -o {outfile}'

        #print(cmd1)
        p1 = subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True).stdout.splitlines()
        #print(p1)
        res1 = p1[0].decode("utf-8")
        print(res1)
