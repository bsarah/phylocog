from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace, add_face_to_node, SVG_COLORS, random_color
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
parser.add_argument("-i", "--intree")#input tree in newick format
parser.add_argument("-o", "--outfile") #output file (plot)
args = parser.parse_args()

treefile = ""
if args.intree:
    treefile = args.intree
else:
    print("no inputtree given!\n")
    exit


outputfile = ""
if args.outfile:
    outputfile = args.outfile
else:
    print("no outputfile given!\n")
    exit




#create dict from Pid to color based on input tree file
#then layout(node) can use the dict to create the layout
colorlist = list(SVG_COLORS)

colorlist.remove('black')
colorlist.remove('white')
colorlist.remove('whitesmoke')
colorlist.remove('seashell')
colorlist.remove('mistyrose')
colorlist.remove('linen')
colorlist.remove('cornsilk')
colorlist.remove('ivory')
colorlist.remove('beige')
colorlist.remove('honeydew')
colorlist.remove('lightcyan')
colorlist.remove('lavenderblush')
colorlist.remove('lemonchiffon')
colorlist.remove('ghostwhite')
colorlist.remove('floralwhite')
colorlist.remove('oldlace')
colorlist.remove('azure')
colorlist.remove('aliceblue')
colorlist.remove('mintcream')
colorlist.remove('snow')
colorlist.remove('antiquewhite')
colorlist.remove('lightyellow')
colorlist.remove('lightgoldenrodyellow')
colorlist.remove('aqua')
colorlist.remove('aquamarine')
colorlist.remove('bisque')
colorlist.remove('blanchedalmond')
colorlist.remove('lightpink')
#colorlist.remove('LightSteelBlue')
#colorlist.remove('LightSalmon')

pid2col = dict()


def layout(node):
    if node.is_leaf():
        # Add node name to leaf nodes
        if("-P" in node.name):
            sn = node.name.split('-')
            pid = ""
            for snelem in sn:
                if len(snelem) > 1 and snelem[0] == 'P':
                    pid = snelem
            curcol = ""
            if(pid in pid2col):
                curcol = pid2col[pid]
            else:
                print(f'WARNING: no color for pid {pid} and node name {node.name}!')
            print(f'color {curcol}')
            F = TextFace(node.name, tight_text=True, fsize=50, ftype="Arial", fgcolor=curcol,bold=True)
            faces.add_face_to_node(F, node, column=0, position="branch-right")
        else:
            print(f'WARNING: no patternID for node {node.name}')
    else:
        F = TextFace(node.name, tight_text=True, fsize=50, ftype="Arial", fgcolor="black",bold=True)
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"
        #C = CircleFace(radius=node.support*75, color="DarkGreen", style="sphere")
        # Let's make the sphere transparent
        #C.opacity = 0.7
        # And place as a float face over the tree
        #faces.add_face_to_node(C, node, 0, position="float")
        faces.add_face_to_node(F, node, column=0, position="branch-right")

    if("-a" in node.name):
        node.set_style(nst2)
    elif("-b" in node.name):
        node.set_style(nst1)
    else:
        node.set_style(nst0)






curcog = ""

ts = treefile.split('/')
tts = ts[-1].split('.')
curcog = tts[0][0:7]

#read tree from file
t = Tree(treefile, format=1)

#bifurcate the tree
root = t.get_tree_root()
root.resolve_polytomy(recursive=False)

curcolnum = 0
#fill dict
for node in t.traverse("postorder"):
    nname = node.name
    if(node.is_leaf()):
        #check Pid
        if("-P" in nname):
            sn = nname.split('-')
            pid = ""
            for snelem in sn:
                if snelem[0] == 'P':
                    pid = snelem
            if pid not in pid2col:
                pid2col[pid] = colorlist[curcolnum]
                curcolnum+=3
    



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
