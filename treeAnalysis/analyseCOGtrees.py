from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace, add_face_to_node
import re
import pandas as pd
import argparse
from io import StringIO

#input:
#list of trees to work on
#infolder (containing the trees)
#translation file for all COGs
#outfolder (for relabeled trees and plots)
#outfilein
#outfilenameout (table listing information for each COG, extend already existing table= outfilein)
#format
#COGid	numA	numB	numTotal	funcCat	proteinName	NumMemClass0	NumMemClass1	NumMemClass2	NumMemClass3	ProfileLength
#+ numsplits avsupport avdistAA avdistBB avdistAB Doverline Donly


#output
#updated translation file
#outtable with information for each COG in the input list 
#+ files for each COG: relabeled tree and plot

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--intrees")#list of trees
parser.add_argument("-j", "--infolder") #location of trees
parser.add_argument("-f", "--infile") #translation file
parser.add_argument("-o", "--outfolder") #location for output files
parser.add_argument("-s", "--stdout",default=1) #only plotting, write further output to stdout, thus next two are not needed. val = 1 or 0
parser.add_argument("-g", "--outfile") #updated translation file
parser.add_argument("-d", "--dvalues") #table including split and dvalues for each COG-tree in the input list
args = parser.parse_args()

treelist = ""
if args.intrees:
    treelist = args.intrees
else:
    print("no inputtree given!\n")
    exit
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


outputfolder = ""
if args.outfolder:
    outputfolder = args.outfolder
else:
    print("no outputfolder given!\n")
    exit


write2stdout = 0
if args.stdout:
    write2stdout = int(args.stdout)
    

outputfile = ""
if args.outfile and write2stdout==0:
    outputfile = args.outfile
elif write2stdout==1:
    print("output to stdout\n")
else:
    print("no outputfile given!\n")
    exit


outputtable = ""
if args.dvalues and write2stdout==0:
    outputtable = args.dvalues
elif write2stdout==1:
    print("output to stdout\n")
else:
    print("no name for the outputtable given!\n")
    exit




def labcount(ls,nodename):
    numc = len(ls)
    labcounts = [0,0,0,0,0] #a,b,s,w, other
    for l in ls:
        if(l == 'a' or "-a" in l):
            labcounts[0]+=1
        elif(l=='b' or "-b" in l):
            labcounts[1]+=1
        elif(l=='s'):
            labcounts[2]+=1
        elif(l=='w'):
            labcounts[3]+=1
        else:
            labcounts[4]+=1
#    print("Labcounts: "+str(labcounts)+" at "+nodename)
    if(sum(labcounts) == labcounts[0]):
        return 'a'
    elif(sum(labcounts) == labcounts[1]):
        return 'b'
    elif(labcounts[3] == sum(labcounts)):
        return 'w'
    elif(labcounts[0] > 0 and labcounts[1] > 0):
        #a and b
        #print("Splitnode1: "+nodename)
        return 's'
    elif(labcounts[0] == 0 and labcounts[1] > 0 and labcounts[2] > 0):
        return 'b'
    elif(labcounts[0] > 0 and labcounts[1] == 0 and labcounts[2] > 0):
        return 'a'
    elif(labcounts[0] == 0 and labcounts[1] > 0 and labcounts[3] > 0):
        return 'b'
    elif(labcounts[0] > 0 and labcounts[1] == 0 and labcounts[3] > 0):
        return 'a'
    else:
        return 'w'


def layout(node):  
    if node.is_leaf():
        # Add node name to leaf nodes
        parts = node.name.split('-')
        if(parts[1][0] == 'a' or "-a" in node.name):
            F = TextFace(node.name, tight_text=True, fsize=60, ftype="Arial", fgcolor="red",bold=True)
            #N = AttrFace("name", fsize=30,fgcolor="red")
        elif(parts[1][0] == 'b' or "-b" in node.name):
            #N = AttrFace("name", fsize=30, fgcolor="blue")
            F = TextFace(node.name, tight_text=True, fsize=60, ftype="Arial", fgcolor="blue",bold=True)
        else:
            F = AttrFace("name", fsize=30, fgcolor="black")
        #faces.add_face_to_node(N, node, 0)
        faces.add_face_to_node(F, node, column=0, position="branch-right")
        #node.set_style(nst0)
    else:
        F = TextFace(node.name, tight_text=True, fsize=60, ftype="Arial", fgcolor="black",bold=True)
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"
        C = CircleFace(radius=node.support*75, color="DarkGreen", style="sphere")
        # Let's make the sphere transparent
        C.opacity = 0.7
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")
        faces.add_face_to_node(F, node, column=0, position="branch-right")
    curmark = nmarks[node.name]
    if(curmark == 'a'):
        node.set_style(nst2)
    elif(curmark == 'b'):
        node.set_style(nst1)
    else:
        node.set_style(nst0)



transupdate = None
cogupdate = None
header = "COGid\tnumA\tnumB\tnumTotal\tNumSplits\tavSupportSplits\tavSupportTree\tavDistAA\tavDistBB\tavDistAB\tDbar\tD\n"
if write2stdout==0:
    #translation file
    transupdate = open(outputfile, 'a')
    #cogfile
    cogupdate = open(outputtable,'a')
    cogupdate.write(header)
else:
    print(header)

#translation file
protID2geneid = dict()
protID2line = dict()

file1 = open(inputfile, 'r')
inlines = file1.readlines()

for line in inlines:
    line2 = line.rstrip()
    cols = line2.split("\t")
    protID2geneid[cols[3]] = cols[0]
    protID2line[cols[3]] = line2

file1.close()


curcog = ""
file3 = open(treelist, 'r')
#trees = file3.readlines()

for line in file3:
    line2 = line.rstrip('\n')
    tfile = inputfolder+"/"+line2
    ts  = line2.split('.')
    curcog = ts[0][0:7]
    #print("COG:")
    print(curcog)
    print(tfile)

    #read tree from file
    t = Tree(tfile, format=1)

    #bifurcate the tree
    root = t.get_tree_root()
    root.resolve_polytomy(recursive=False)




    #postorder traversal
    id = 0
        ##calculate splits
    nmarks = dict() #node identifier to marks a, b, s or w
    numsplits = 0
    splitsnodes = []
    mixednodes = []
    sumsupport = 0
    numinnodes = 0
    anodes = [] #names of archaeal nodes
    bnodes = [] #names of bacterial nodes
    totnodes = 0

    for node in t.traverse("postorder"):
        nname = node.name
        totnodes+=1
        if(node.is_leaf()):
            protid = re.sub(r"\s+", '_', nname)
            taxid = protID2geneid[protid]
            #find taxid
            tmpname = str(id)+'-'+taxid
            node.name = tmpname
            #update translation file
            trline = protID2line[protid].split("\t")
            trline[0] = tmpname
            outline = "\t".join(trline)
            if write2stdout==0:
                transupdate.write(outline+"\n")
            else:
                print(outline,"\n")
            parts = tmpname.split('-')
            #        print(parts[1][0])

            if(parts[1][0] == 'b' or "-b" in tmpname):
                bnodes.append(node.name)
                nmarks[tmpname] = "b"
            else:
                anodes.append(node.name)
                nmarks[tmpname] = "a"
        else:
            if(node.is_root()):
                tmpname = "r"+str(id)
            else:
                tmpname = "i"+str(id)+"i"
            node.name = tmpname
            node.support = 0.0
            if(nname != None and nname != ''):
                node.support = float(nname)
            sumsupport+= node.support
            numinnodes+= 1
            cs = [] #labels at children nodes
            for c in node.children:
                cs.append(nmarks[c.name])
            res = labcount(cs,tmpname) 
            nmarks[tmpname] = res
            if(res == 's'):
                numsplits+=1
                splitsnodes.append(tmpname)
            if(res == 'w'):
                #print(node.name)
                mixednodes.append(tmpname)
        id=id+1


    print(f'anodes: {len(anodes)}, bnodes: {len(bnodes)}')

    #iterate through splitsnodes and check which of the children will be cut, as the children have two distinct labels, we chose the one to keep by having the same label as the brother, thus we cut the other child
    cutnodes = []
    sumsplitsupport = 0
    for sn in splitsnodes:
        #sn = name of current splitsnode
        wsn = t&sn
        cs = wsn.children
        if(wsn.is_root()):
            cutnodes.append(cs[0].name)
            sumsplitsupport+= cs[0].support
        else:
            curpar = wsn.up
            sibs = curpar.children
            #node has a sibling:
            if(len(sibs) == 2):
                siblab = ''
                if(sibs[0].name == sn):
                    siblab = nmarks[sibs[1].name]
                else:
                    siblab = nmarks[sibs[0].name]
                    #there should be two children!
                    #we take the one that has another label!
                if(nmarks[cs[0].name] == siblab):
                    cutnodes.append(cs[1].name)
                    sumsplitsupport+= cs[1].support
                else:
                    cutnodes.append(cs[0].name)
                    sumsplitsupport+= cs[0].support
            else:#no sibling
                #take any node?
                cutnodes.append(cs[0].name)
                sumsplitsupport+= cs[0].support
    #print("Number of splits: "+ str(numsplits))
    #print(splitsnodes)
    #print(cutnodes)
    avsupport = 0
    if(numsplits>0):
        avsupport = sumsplitsupport/numsplits
    #print(avsupport)
    avallsupport = 0
    if(numinnodes>0):
        avallsupport = sumsupport/numinnodes
    #print(avallsupport)


    ##plot tree
    #define more node styles
    nst0 = NodeStyle()
    nst0["hz_line_width"]=5
    nst0["vt_line_width"]=5
    nst1 = NodeStyle()
    nst1["bgcolor"] = "LightSteelBlue"
    nst1["hz_line_width"]=5
    nst1["vt_line_width"]=5
    nst2 = NodeStyle()
    nst2["bgcolor"] = "LightSalmon"
    nst2["hz_line_width"]=5
    nst2["vt_line_width"]=5
    #nst3 = NodeStyle()
    #nst3["bgcolor"] = "DarkSeaGreen"
    #nst4 = NodeStyle()
    #nst4["bgcolor"] = "Khaki"





    #only create plots when number of nodes is below 150 nodes
    #check with tree data structure:
    
    if(totnodes <= 10250):
        print("plotting..")
    
        #plotting
        plotfile = outputfolder+"/"+curcog+"_colors.pdf"

        ts = TreeStyle()
        #ts.legend.add_face(CircleFace(10, "red"), column=1)
        ts.show_leaf_name = False
        ts.show_branch_length = True
        ts.show_branch_support = True
        ts.mode = "c"
        ts.layout_fn = layout
        ts.arc_start = 0 # 0 degrees = 3 o'clock
        ts.arc_span = 360


        # Add legend
        legend={}
        ts.legend.add_face(CircleFace(100, "red"), column=0)
        ts.legend.add_face(TextFace("  Archaea", fsize=200, ftype="Arial"), column=1)
        
        ts.legend.add_face(CircleFace(100, "blue"), column=0)
        ts.legend.add_face(TextFace("  Bacteria", fsize=200, ftype="Arial"), column=1)
        
        ts.legend.add_face(CircleFace(100, "green"), column=0)
        ts.legend.add_face(TextFace("  support", fsize=200, ftype="Arial"), column=1)
        
        # Position legend in lower right corner
        ts.legend_position=4

        # Add title
        ts.title.add_face(TextFace(curcog, fsize=300, ftype="Arial", fstyle="italic"), column=0)

        # Additional info
        #ts.title.add_face(TextFace("Total spent by referred customers: £" + str(int(data["Amount spent (£)"].sum())), fsize=20, ftype="Arial", fstyle="italic"), 
        #  column=0)
        #ts.scale = 1000
        ts.optimal_scale_level = 'full'
        ts.allow_face_overlap = True

        ts.show_scale = True
        #t.show(tree_style=ts)
        ts.scale = 700
        #t.render("%%inline", tree_style=ts, w=800)
        t.render(plotfile, w=800, tree_style=ts)



    #calculate distances
    print("calculate pairwise distances..")
    #average pairwise A-A distance
    distsumaa = 0
    distnumaa = 0
    for i in range(0,len(anodes)):
        an = anodes[i]
        va = t&an
        #va = t.get_leaves_by_name(an)
        for j in range(i+1,len(anodes)):
            am = anodes[j]
            wa = t&am
            if(am != an):
                da = t.get_distance(va,wa)
                distsumaa = distsumaa + da
                distnumaa += 1
    avdistaa = 0
    if(distnumaa>0):
        avdistaa = distsumaa/distnumaa
    #average pairwise B-B distance
    distsumbb = 0
    distnumbb = 0
    for i in range(0,len(bnodes)):
        bn = bnodes[i]
        vb = t&bn
        for j in range(i+1,len(bnodes)):
            bm = bnodes[j]
            wb = t&bm
            if(bm != bn):
                db = t.get_distance(vb,wb)
                distsumbb = distsumbb + db
                distnumbb += 1
    avdistbb = 0
    if(distnumbb>0):
        avdistbb = distsumbb/distnumbb
    #average pairwise A-B distance
    distsumab = 0
    distnumab = 0
    for aa in anodes:
        va = t&aa
        for bb in bnodes:
            vb = t&bb
            dab = t.get_distance(va,vb)
            distsumab = distsumab + dab
            distnumab += 1
    avdistab = 0
    if(distnumab>0):
        avdistab = distsumab/distnumab
    #print("AA distances: " + str(avdistaa) + ", " + str(distsumaa) + ", "+ str(distnumaa))
    #print("BB distances: " + str(avdistbb) + ", " + str(distsumbb) + ", "+ str(distnumbb))
    #print("AB distances: " + str(avdistab) + ", " + str(distsumab) + ", "+ str(distnumab))

    #calculate D values
    dol=0
    don=0
    if(avdistab>0):
        #doverline = 0.5*(avdaa+avdbb)/avdab
        dol = 0.5*(avdistaa+avdistbb)/avdistab
        #d = min(avdaa,avdbb)/avdab
        don = (min(avdistaa,avdistbb))/avdistab

    #print(dol)
    #print(don)


    #add infos: numsplits avdistAA avdistBB avdistAB Doverline Donly
    curline = f'{curcog}\t{len(anodes)}\t{len(bnodes)}\t{len(anodes)+len(bnodes)}'
    cogoutline = curline+"\t"+str(numsplits)+"\t"+str(avsupport)+"\t"+str(avallsupport)+"\t"+str(avdistaa)+"\t"+str(avdistbb)+"\t"+str(avdistab)+"\t"+str(dol)+"\t"+str(don)+"\n"
    if write2stdout==0:
        cogupdate.write(cogoutline)
    else:
        print(cogoutline)


if write2stdout==0:
    cogupdate.close()
    transupdate.close()
