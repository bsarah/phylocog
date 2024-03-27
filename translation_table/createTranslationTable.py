import argparse


#this program combines two cog-20 files and outputs a translation table
#to translate between protein ID and tax ID and domain (a/b)
#to get the domain letter (a/b) check in input 3 using the taxid if it is ARCHAEA or BACTERIA

#input cog-20.cog.csv
#format: 
#Columns:
#
#1.	Gene ID (GenBank or ad hoc)
#
#2.	NCBI Assembly ID
#
#3.	Protein ID (GenBank if conforms to [A-Za-z0-9_]+\.[0-9]+ regex; ad hoc otherwise)
#
#4.	Protein length
#
#5.	COG footprint coordinates on the protein. "201-400" means "from position 201 to position 400"; "1-100=201-300" indicates a segmented footprint, 1-100 AND 201-300
#
#6.	Length of the COG footprint on the proteins
#
#7.	COG ID
#
#8.	reserved
#
#9.	COG membership class (0: footprint covers most of the protein and most of the COG profile; 1: footprint covers most of the COG profile and part of the protein; 2: footprint covers most of the protein and part of the COG profile; 3: partial match on both protein and COG profile)
#
#10.	PSI-BLAST bit score for the match between the protein and COG profile
#
#11.	PSI-BLAST e-value for the match between the protein and COG profile
#
#12.	COG profile length
#
#13.	Protein footprint coordinates on the COG profile


#input2: cog-20.org.csv
#format: 
#Columns:
#
#1.	NCBI Assembly ID
#
#2.	Organism (genome) name
#
#3.	NCBI Tax ID of the assembly
#
#4.	Taxonomic category used in COGs


#input3: cog-20.tax.csv
#Columns:
#
#1.	Taxonomic category in COGs
#
#2.	Parent taxonomic category (self, if top of the hierarchy)
#
#3.	NCBI Tax ID of the assembly


#output: combination of the two tables, just add matching rows of tab2 to tab1
#additionally, add as a first column the tax id with a or b on first position
#eg a01234, b56789

#output: cog-20.translation.csv
#Columns:
#
#0.     DomainLetter+Taxid
#
#1.	Gene ID (GenBank or ad hoc)
#
#2.	NCBI Assembly ID
#
#3.	Protein ID (GenBank if conforms to [A-Za-z0-9_]+\.[0-9]+ regex; ad hoc otherwise)
#
#4.	Protein length
#
#5.	COG footprint coordinates on the protein. "201-400" means "from position 201 to position 400"; "1-100=201-300" indicates a segmented footprint, 1-100 AND 201-300
#
#6.	Length of the COG footprint on the proteins
#
#7.	COG ID
#
#8.	reserved
#
#9.	COG membership class (0: footprint covers most of the protein and most of the COG profile; 1: footprint covers most of the COG profile and part of the protein; 2: footprint covers most of the protein and part of the COG profile; 3: partial match on both protein and COG profile)
#
#10.	PSI-BLAST bit score for the match between the protein and COG profile
#
#11.	PSI-BLAST e-value for the match between the protein and COG profile
#
#12.	COG profile length
#
#13.	Protein footprint coordinates on the COG profile
#
#14.	Organism (genome) name
#
#15.	NCBI Tax ID of the assembly
#
#16.	Taxonomic category used in COGs


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cog20cog")#input1
parser.add_argument("-o", "--cog20org")#input2
parser.add_argument("-t", "--cog20tax")#input3
parser.add_argument("-s", "--cog20translation") #this is the output file
args = parser.parse_args()

input1 = ""
if args.cog20cog:
    input1 = args.cog20cog
else:
    print("no cog20cog given!\n")
    exit

input2 = ""
if args.cog20org:
    input2 = args.cog20org
else:
    print("no cog20org given!\n")
    exit

input3 = ""
if args.cog20tax:
    input3 = args.cog20tax
else:
    print("no cog20tax given!\n")
    exit

outputfile = ""
if args.cog20translation:
    outputfile = args.cog20translation
else:
    outputfile = "cog-20.translation.csv"



taxid2letter = dict()
spec2letter = dict()
file3 = open(input3, 'r',encoding="utf8", errors='ignore')
taxlines = file3.readlines()

for line in taxlines:
    line2 = line.strip()
    cols = line2.split(",")
    mylet = "x"
    myspec = cols[0]
    if(cols[1] == "ARCHAEA"):
        mylet = "a"
    elif(cols[1] == "BACTERIA" or cols[1] == "PROTEOBACTERIA" or cols[1] == "FIRMICUTES"):
        mylet = "b"
    else:
        print(f'no clear letter for {cols[1]}\n')
    taxid2letter[cols[2]] = mylet
    spec2letter[myspec] = mylet
file3.close()

protid2line = dict()
protid2taxletter = dict()
file2 = open(input2, 'r',encoding="utf8", errors='ignore')
orglines = file2.readlines()

for gline in orglines:
    gline2 = gline.strip()
    #print(gline2)
    gls = gline2.split(",")
    if(gls[0] in protid2line):
        print(f'{gls[0]} exists twice!\n')
    else:
        protid2line[gls[0]] = f'{gls[1]}\t{gls[2]}\t{gls[3]}'
    taxid = gls[2]
    spec = gls[-1]
    if(spec in spec2letter):
        curlet = spec2letter[spec]
        curtaxid = f'{curlet}{taxid}'
        protid2taxletter[gls[0]] = curtaxid
    else:
        print(f'{spec} not in tax file or dict!\n')
        curtaxid = f'x{taxid}'
        protid2taxletter[gls[0]] = curtaxid

file2.close()


outf = open(outputfile, 'a')

file1 = open(input1, 'r',encoding="utf8", errors='ignore')
for cline in file1:
    cline2 = cline.strip()
    cls = cline2.split(",")
    protid = cls[1]
    clstr = "\t".join(cls)
    if(protid in protid2line):
        preline = protid2line[protid]
        taxlet = protid2taxletter[protid]
        outline = f'{taxlet}\t{clstr}\t{preline}\n'
        outf.write(outline)
