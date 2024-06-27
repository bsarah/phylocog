from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace, add_face_to_node, SVG_COLORS, random_color
import re
import pandas as pd
import argparse
from io import StringIO


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--intree")#input tree in newick format
parser.add_argument("-o", "--outfile") #output file (plot)
args = parser.parse_args()



treefile = ""
if args.intree:
    treefile = args.intree
else:
    print("no inputtree given!\n")
    exit

#dotranslate = False
#translationfile = ""
#if args.transfile:
#    translationfile = args.transfile
#    dotranslate = True


outputfile = ""
if args.outfile:
    outputfile = args.outfile
else:
    print("no outputfile given!\n")
    exit





colorlist = list(SVG_COLORS)



def labcount(ls,nodename):
    numc = len(ls)
    labcounts = [0,0,0,0,0] #a,b,s,w, other
    for l in ls:
        if(l == 'A' or "-A" in l):
            labcounts[0]+=1
        elif(l=='B' or "-B" in l):
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
        #parts = node.name.split('-')
        if("-A" in node.name):
            F = TextFace(node.name, tight_text=True, fsize=60, ftype="Arial", fgcolor="red",bold=True)
            #N = AttrFace("name", fsize=30,fgcolor="red")
        elif("-B" in node.name):
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



def calcSplits(tree, nmarks):
    #postorder traversal
    id = 0
    ##calculate splits

    numsplits = 0
    splitsnodes = []
    mixednodes = []
    sumsupport = 0
    numinnodes = 0
    anodes = [] #names of archaeal nodes
    bnodes = [] #names of bacterial nodes
    totnodes = 0
    numaBnodes =0
    numbAnodes = 0
    others = 0

    for node in t.traverse("postorder"):
        nname = node.name
        totnodes+=1
        if(node.is_leaf()):
            sn = nname.split('-')
            tmpname = str(id)+'-'+nname
            node.name = nname
            
            if("-b" in tmpname):
                bnodes.append(node.name)
                nmarks[node.name] = "b"
                if(tmpname[-2:] == "-A"):
                    numbAnodes+=1
            elif("-a" in tmpname):
                anodes.append(node.name)
                nmarks[node.name] = "a"
                if(tmpname[-2:] == "-B"):
                    numaBnodes+=1
            else:
                others += 1
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
    curline = f'{curcog}\t{len(anodes)}\t{len(bnodes)}\t{len(anodes)+len(bnodes)}\t{numaBnodes}\t{numbAnodes}\t{numsplits}\t{others}\n'
    return curline

        
curcog = ""

ts = treefile.split('/')
tts = ts[-1].split('.')
curcog = tts[0][0:7]

#read tree from file
t = Tree(treefile, format=1)

#bifurcate the tree
root = t.get_tree_root()
root.resolve_polytomy(recursive=False)


nmarks = dict() #node identifier to marks a, b, s or w
outline = calcSplits(t, nmarks)

print(outline)

lso = outline.split("\t")
totnodes = int(lso[-2])

if(totnodes < 150):
    

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
    
    
    
    #check with tree data structure:
    
    print("plotting..")
    
    #plotting
    
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
    ts.legend.add_face(CircleFace(100, "lightsalmon"), column=0)
    ts.legend.add_face(TextFace("  Archaea", fsize=200, ftype="Arial"), column=1)
    
    ts.legend.add_face(CircleFace(100, "lightsteelblue"), column=0)
    ts.legend.add_face(TextFace("  Bacteria", fsize=200, ftype="Arial"), column=1)
    
    #ts.legend.add_face(CircleFace(100, "green"), column=0)
    #ts.legend.add_face(TextFace("  support", fsize=200, ftype="Arial"), column=1)
    
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
    t.render(outputfile, w=800, tree_style=ts)

