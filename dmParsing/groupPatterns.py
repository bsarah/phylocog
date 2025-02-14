import argparse

#import numpy as np
import pandas as pd

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


parser = argparse.ArgumentParser(description="Separate domain sequences into groups of patterns")
parser.add_argument("infile", help = "DM output file")
parser.add_argument("-s", "--seqfile", help = "optional, fasta file with corresponding protein sequences")
parser.add_argument("-t", "--transfile", help = "optional, translation file to translate between protein IDs and tax IDs as in the trees.")
parser.add_argument("-o","--outname", help = "name for produced files as output")
parser.add_argument("-p","--patternlist", help = "optional, output file listing all patterns with counts for archaeal/bacterial proteins")

args = parser.parse_args()

inputfile = ""
if args.infile:
    inputfile = args.infile
else:
    print("no inputfile given!\n")
    exit

doseqlen = False
sequencefile = ""
if args.seqfile:
    sequencefile = args.seqfile
    doseqlen= True
#else:
#    print("no fasta file given!\n")
#    exit


dotranslate = False
translationfile = ""
if args.transfile:
    translationfile = args.transfile
    dotranslate = True


outputname = inputfile+".out"
if args.outname:
    outputname = args.outname


dodetails = False
patternfile = ""
if args.patternlist:
    patternfile = args.patternlist
    dodetails = True


class SDomain:
    """a structural domain with corresponding information"""

    def __init__(self, ann, acc, low, up, fid, evalue, subdoms, numdoms, prop):
        self.ann = ann
        self.acc = acc #accession number
        self.low = low #in case of several subdomains, this is the lowest coordinate in all of them
        self.up = up #analoguous to low
        self.fid = fid
        self.evalue = evalue
        self.subdoms = subdoms # this is a list with coordinate pairs for the subdomains, thus (start,end)
        self.numdoms = numdoms #number of subdomains
        self.prop= prop

    def __eq__(self,sdom2):
        if (self.ann == sdom2.ann and
            self.acc == sdom2.acc and
            self.low == sdom2.low and
            self.up == sdom2.up and
            self.fid == sdom2.fid and
            self.evalue == sdom2.evalue and
            self.subdoms == sdom2.subdoms and
            self.numdoms == sdom2.numdoms):
            return True
        else:
            return False

    def overlap11(self,sdom2):
        #same number of subdomains that overlap with each other
        #only makes sense when they have the same fid
        if(self.numdoms != sdom2.numdoms):
            #print(f'different numdoms: {self.numdoms}, {sdom2.numdoms}')
            return False
        if(self.fid != sdom2.fid):
            #print("different fids")
            return False
        ocount = 0
        for i in range(len(self.subdoms)):
            (s1, e1) = self.subdoms[i]
            (s2, e2) = sdom2.subdoms[i]
            if((s1 <= e2) and (e1 >= s2)):
                ocount+=1
            #else:
                #return False
        if(ocount == 0):
            return False
        return True
        

    def overlap1n(self,sdom2):
        #different numbers of subdomains, we only want all of the smaller subdomains to all have an overlap
        #only makes sense when they have the same fid
        if(self.fid != sdom2.fid):
            return False
        
        mindoms = []
        maxdoms = []
        if(self.numdoms < sdom2.numdoms):
            mindoms = self.subdoms
            maxdoms = sdom2.subdoms
        else:
            maxdoms = self.subdoms
            mindoms = sdom2.subdoms
            #now, are we doing alignments?
            #TODO: finish!
            #not yet, at the moment, we want to have the same annotation

##################end class SDomain




##########################start class SSequence#############
#this class will join all the SDomains of the same accession such that they all stay together
class SSequence:
    """A group of SDomains belonging to the same accession number (= protein)"""
    """joining two sseqs will result in a SCat object"""

    def __init__(self, sdomlist, ssid):
        self.sdomlist = sorted(sdomlist, key = lambda x: x.low, reverse = False)
        self.ssid = ssid #artifical id for identification purposes
        self.ann = self.createAnnStr()
        self.acc = sdomlist[0].acc #should be the same for all SDomains
        self.low = sdomlist[0].low
        self.up = sdomlist[-1].up

    def __eq__(self,sseq2):
        if(self.sdomlist == sseq2.sdomlist):
            return True
        else:
            return False

    def numSdoms(self): #==num contained sdomains
        return len(self.sdomlist)


    def createAnnStr(self):
        if(len(self.sdomlist) == 0):
            return ""
        elif(len(self.sdomlist) == 1):
            return self.sdomlist[0].ann
        else: #merge different annotations by sorting them accordinly to their ranges
            fid2low = dict()
            fid2ann = dict()
            for c in self.sdomlist:
                fid2low[c.fid] = c.low
                fid2ann[c.fid] = c.ann
            sorted_fid_by_low = sorted(fid2low.items(), key=lambda x:x[1])
            curann = ""
            for (fid,low) in sorted_fid_by_low:
                if(fid2ann[fid] != "END"):
                    if(curann == ""):
                        curann += f'{fid2ann[fid]}'
                    else:
                        curann += f'={fid2ann[fid]}'
            return curann

    def getFIDs(self):
        fids = []
        for s in self.sdomlist:
            fids.append(s.fid)
        return list(set(fids))

    def createDIDStr(self):
        if(len(self.sdomlist) == 0):
            return ""
        elif(len(self.sdomlist) == 1):
            return self.sdomlist[0].fid
        else: #put all fids and annotations in the order as a string together
            curstr = ""
            for c in self.sdomlist:
                if(c.ann != "END"):
                    curstr += c.fid+c.ann+'_'
                else:
                    curstr += c.fid+c.ann
            return curstr



    def createDIDStrLength(self):
        if(len(self.sdomlist) == 0):
            return ""
        elif(len(self.sdomlist) == 1):
            return self.sdomlist[0].fid
        else: #put all fids and annotations in the order as a string together
            curstr = ""
            totlen = 0 #check the summed up length of the sdomains to sort into shorter and longer ones
            for c in self.sdomlist:
                totlen += abs(c.up-c.low)
                if(c.ann != "END"):
                    curstr += c.fid+c.ann+'_'
                else:
                    curstr += c.fid+c.ann
            nstr = str(totlen)
            nstrl = len(nstr)
            #take 9 positions
            missing = 9-nstrl
            zerolist = ['0']*missing
            zlstr = ''.join(zerolist)
            curstr += '_'+zlstr+nstr
            return curstr

##########################end class SSequence#############




class SClus:
    """This class will summarize SSeq objects that can be joined (= have exactly the same sequence of
    SDomains. In this way, we always keep all the SDomains from the same accession together and only
    join exactly the same patterns, except for the lengths of the single Sdomains.
    attributes: artificial cid, min low, max up, list sseqs, list av_subranges list of fids, list of accs, ann
    list of avsubranges and fids refer to the sequences of SDomains in the SClus that should be the
    same for all the accessions
    sdomains in one sseq are sorted increasingly by low value
    """
    
    def __init__(self, low, up, sseqs, uid, ann):
        self.ann = ann
        self.low = low #min values of all
        self.up = up #max value of all
        self.sseqs = sseqs #list of contained ssequences
        self.uid = uid
        self.accs = sseqs[0].getFIDs()    #getAccList()
#        self.num11map = 0 #domains that are added and map 1:1
#        self.num1nmap = 0 #domains that are added and map 1:n
        #summarize ranges for subdomains into average ranges to be plotted
        #all domains with the same pattern should have the same number of subdomains with overlapping ranges
        self.avsubranges = []
        self.fids = []
        numsubdoms = len(self.sseqs[0].sdomlist)
        for s in range(numsubdoms):
            suma = 0
            sumb = 0
            curfid = self.sseqs[0].sdomlist[s].fid
            for i in range(len(sseqs)):
                (a,b) = (self.sseqs[i].sdomlist[s].low,self.sseqs[i].sdomlist[s].up)
                suma += a
                sumb += b
                #print(f'subsubdoms {a},{b}')
            ava = int(suma/numsubdoms)
            avb = int(sumb/numsubdoms)
            #print(f'{fid} subdomain ({ava},{avb})')
            self.avsubranges.append((ava,avb))
            self.fids.append(curfid)
        
    def __eq__(self,clus2):
        #uid doesn't have to be the same, in this way, we can check if the same cluster has two different uids
        if(self.ann == clus2.ann and
           self.low == clus2.low and
           self.up == clus2.up and
           self.fids == clus2.fids):
            return True
        else:
            return False

    def getAccList(self):
        res = []
        for sd in self.sseqs:
            res.append(sd.acc)
        res2 = list(set(res))
        return res2



    def addSSeq(self,sseq):   
        if(sseq.acc in self.accs):
#            print(f'addSSeq 0: {sseq.acc}')
            return False
        #check conditions to attachment first
        checkann = (sseq.ann == self.ann)
        checkfids = (set(self.fids) == set(sseq.getFIDs()))
        #check lengths of sdoms and also check if up and low are very far from the current average

        #what is very far? 50?
        checklens = True
        distthreshold = 75
        avlens = [0]*len(sseq.sdomlist)
        
        for sq in self.sseqs:
            for iss in range(len(sq.sdomlist)):
                sd = sq.sdomlist[iss]
                lsd = sd.up - sd.low
                avlens[iss] += lsd

        avlens = [a/len(self.sseqs) for a in avlens]

        for sis in range(len(sseq.sdomlist)):
            ssd = sseq.sdomlist[sis]
            ssdlen = ssd.up-ssd.low
            delta = abs(ssdlen-avlens[sis])
            if(delta > distthreshold):
               checklens = False 
                
                    
        if(checkann and checkfids and checklens):
            #check if it will change up and low (max, min)
            if(sseq.low < self.low):
                self.low = sseq.low
            if(sseq.up > self.up):
                self.up = sseq.up
            #add sseq to current SClus
            #take new average values for low and up for each sdomain
            self.sseqs.append(sseq)
            newavsubranges = []
            numsubdoms = len(self.sseqs[0].sdomlist)
            for s in range(numsubdoms):
                suma = 0
                sumb = 0
                for i in range(len(self.sseqs)):
                    (a,b) = (self.sseqs[i].sdomlist[s].low,self.sseqs[i].sdomlist[s].up)
                    suma += a
                    sumb += b
                ava = int(suma/numsubdoms)
                avb = int(sumb/numsubdoms)
                #print(f'{fid} subdomain ({ava},{avb})')
                newavsubranges.append((ava,avb))
            self.avsubranges = newavsubranges
            return True
        else:
#            print(f'anns {sseq.ann} and {self.ann}; fids: {sseq.getFIDs()}  {self.fids} ')
            return False


####################end class SClus#########


class Cluster:
    """A cluster of structural domain annotations for a COG, different accessions but always same fid"""

    #def __repr__(self):
    #    return f"Ship('{self.name}', {self.positions})"

    #def __str__(self):
    #    return f'{repr(self)} with hits {self.hits}'

    def __init__(self, fid, low, up, sdomains, uid, ann):
        self.fid = fid #all sdomains have the same fid inside a cluster
        self.ann = ann
        self.low = low #min values of all
        self.up = up #max value of all
        self.sdomains = sdomains #list of contained sdomains
        self.uid = uid
        self.num11map = 0 #domains that are added and map 1:1
        self.num1nmap = 0 #domains that are added and map 1:n
        #summarize ranges for subdomains into average ranges to be plotted
        #all domains in the cluster should have the same number of subdomains with overlapping ranges
        self.avsubranges = []
        numsubdoms = len(self.sdomains)
        numsdomainparts = self.sdomains[0].numdoms
        #print(f'numsubdoms {numsubdoms}')
        for s in range(numsdomainparts):
            suma = 0
            sumb = 0
            for i in range(numsubdoms):
                (a,b) = self.sdomains[i].subdoms[s]
                suma += a
                sumb += b
                #print(f'subsubdoms {a},{b}')
            ava = int(suma/numsubdoms)
            avb = int(sumb/numsubdoms)
            #print(f'{fid} subdomain ({ava},{avb})')
            self.avsubranges.append((ava,avb))
        
    def __eq__(self,clus2):
        #uid doesn't have to be the same, in this way, we can check if the same cluster has two different uids
        if(self.fid == clus2.fid and
           self.low == clus2.low and
           self.up == clus2.up and
           self.ann == clus2.ann and
           self.sdomains == clus2.sdomains):
            return True
        else:
            return False

    def getAccList(self):
        res = []
        for sd in self.sdomains:
            res.append(sd.acc)
        res2 = list(set(res))
        return res2


    def addSdom11(self,sdomain):
        numol = 0
        for sd in self.sdomains:
            if(sd.overlap11(sdomain)):
                numol+=1
        if(numol == len(self.sdomains)):
            self.sdomains.append(sdomain)
            self.num11map+=1
            if(sdomain.low < self.low):
                self.low = sdomain.low
            if(sdomain.up > self.up):
                self.up = sdomain.up
            return True
        return False

    def canbejoined(self,clus2):
        #print("cluster can be joined")
        if(self.fid != clus2.fid):
            return False
        if(len(self.sdomains) != len(clus2.sdomains)):
            return False
        for curs in clus2.sdomains:
            numol = 0
            for sd in self.sdomains:
                if(sd.overlap11(curs)):
                    numol+=1
            if(numol != len(self.sdomains)):
                return False
        return True

    def join(self,clus2):
        #for all sdomains in clus2, they will be joined
        for curs in clus2.sdomains:
            if(self.addSdom11(curs) == False):
                return False
        return True

    

    def addSdom1n(self,sdomain):
        numol = 0
        for sd in self.sdomains:
            if(sd.overlap1n(sdomain)):
                numol+=1
        if(numol == len(self.sdomains)):
            self.sdomains.append(sdomain)
            self.num1nmap+=1
            if(sdomain.low < self.low):
                self.low = sdomain.low
            if(sdomain.up > self.up):
                self.up = sdomain.up
            return True
        return False

    #join two clusters: rather add one sequence at a time?
    #can we join without additional checks when assuming that both clusters map 1:1 and the low and up of the clusters are similar?


####################end class Cluster

class SCat:
    #needed?
    
    def __init__(self, sseqlist, scid):
        self.sseqlist = sseqlist
        self.scid = scid #artifical id for identification purposes
        tmpacc = []
        for c in self.sseqlist:
            tmpacc.append(c.acc)
        self.accs = list(set(tmpacc)) #different accession but same sseq pattern
        self.ann = sseqlist[0].ann #should be the same for all the ssequences here
        self.numsseq = len(accs)

    def __eq__(self,scat2):
        if(set(self.accs) == set(cat2.accs) and self.ann == scat2.ann):
            return True
        else:
            return False

##needed?

#######################end class SCat ###########



class Category:
    """this class is able two join several clusters with different Fids but
    with the same set of accessions, thus we can plot the categories and
    summarize accessions with the same sdomain sequence including different
    fids"""

    def __init__(self, cluslist, cid):
        self.cluslist = cluslist
        self.cid = cid #artifical id for identification purposes
        tmpacc = []
        for c in self.cluslist:
            curacc = c.getAccList()
            tmpacc += curacc
        self.accs = list(set(tmpacc))
        self.ann = self.createAnnStr()

    def __eq__(self,cat2):
        if(self.cluslist == cat2.cluslist):
            return True
        else:
            return False

    def add(self, clus): #we want to add another cluster to the category
        #the cluster has to have the same set of accessions than the category
        clusaccs = clus.getAccList()
        if(set(clusaccs) == set(self.accs) or len(self.accs)==0):
            joined = False
            for c in self.cluslist:
                if(c.canbejoined(clus)):
                    res = c.join(clus)
                    joined = True
                    break
            if(joined == False):
                #add as a separate cluster
                self.cluslist.append(clus)
            self.ann = self.createAnnStr()
            return True
        return False

    def canbejoined(self,cat2): #can only be joined if we can successfully join all clusters
        if(len(self.cluslist) != len(cat2.cluslist)):
            return False
        foundpartner = 0
        for c in self.cluslist:
            for c2 in cat2.cluslist:
                if(c.canbejoined(c2)):
                    foundpartner+=1
        if(foundpartner >= len(self.cluslist)):
            return True
        return False

    def join(self,cat2): #what happens to the cid if two clusters will be joined?
        #return a new cluster here that is the join of the initial ones?
        if(len(self.cluslist) != cat2.cluslist):
            return False
        foundpartner = 0
        for c in self.cluslist:
            for c2 in cat2.cluslist:
                if(c.join(c2)):
                    foundpartner+=1
        if(foundpartner >= len(self.cluslist)):
            return True
        return False

    def numFid(self): #==num contained clusters
        return len(self.cluslist)

    def createAnnStr(self):
        if(len(self.cluslist) == 0):
            return ""
        elif(len(self.cluslist) == 1):
            return self.cluslist[0].ann
        else: #merge different annotations by sorting them accordinly to their ranges
            uid2low = dict()
            uid2ann = dict()
            for c in self.cluslist:
                uid2low[c.uid] = c.low
                uid2ann[c.uid] = c.ann
            sorted_uid_by_low = sorted(uid2low.items(), key=lambda x:x[1])
            curann = ""
            for (uid,low) in sorted_uid_by_low:
                if(uid2ann[uid] != "END"):
                    curann += uid2ann[uid]
            return curann


############end class Category


def writetofileverbose(catlist, outputname, pid2length):
    outputfile = outputname
    outf = open(outputfile,"a")
    #verbose because we don't write the summarized/average start/end position but all the details of every single sdomain

    #NEXT

    #write mapping to accessions in a separate file that is easily parseable
    #format:

    #>artificial catID and number of sequences and structural domain annotation as in plot
    ##list of sequences contained:
    #Accession [(0,0,START), (start,end,Fid)...(n,n,END)] :n = len(prot)
    ##...

    #>artificial catID and number of sequences and structural domain annotation as in plot
    #Accession [(0,0,START), (start,end,Fid)...(n,n,END)] :n = len(prot)
    ##...
    starttuple = (0,0,"START")
    for cc in range(len(catlist)):
        c = catlist[cc] 
        acc2ll = dict()
        idline = f'>PatternID {c.cid}; #Accessions {len(c.accs)}; Structural_annotation {c.ann}\n'
        if(c.ann == "END"):
            continue
        for s in c.cluslist: #s=cluster
            for t in s.sdomains: #t=sdomain
                if(t.ann == "END"):
                    continue
                sacc = t.acc
                endcoord = pid2length[sacc]
                if(sacc not in acc2ll):
                    acc2ll[sacc] = []
                    acc2ll[sacc].append(starttuple)
                    acc2ll[sacc].append((endcoord,endcoord+10,"END"))
                for (u,v) in t.subdoms:
                    acc2ll[sacc].append((u,v,t.fid))
        #sort ll by first tuple entry
        outf.write(idline)
        #print(idline)
        for (k,v) in acc2ll.items():
            v.sort(key=lambda a: a[0])
            vstr = str(k)
            for (a,b,c) in v:
                vstr+=f' ({a},{b},{c})'
            outf.write(f'{vstr}\n')
            #print(f'{vstr}\n')


#======================================end writetofile



def writeSClustofileverbose(scluslist, outputname, pid2length, protID2geneid):
    outputfile = outputname
    outf = open(outputfile,"a")
    #write mapping to accessions in a separate file that is easily parseable
    #format:

    #>artificial catID and number of sequences and structural domain annotation as in plot
    ##list of sequences contained:
    #Accession [(0,0,START), (start,end,Fid)...(n,n+10,END)] :n = len(prot)
    ##...

    #>artificial catID and number of sequences and structural domain annotation as in plot
    #Accession [(0,0,START), (start,end,Fid)...(n,n+10,END)] :n = len(prot)
    ##...
    #starttuple = (0,0,"START")
    for cc in range(len(scluslist)):
        c = scluslist[cc] #current pattern
        if(c.ann == "END"):
            continue
        annlist = c.ann.split('=')
        correctannlist = []
        curidx = 0
        anndone = 0 #only needed for first sequence in pattern
        numa = 0
        numb =0   
        for s in c.sseqs: #s=sseq
            #we have to write down all the SSequences per cluster
            vstr = str(s.acc)
            if(s.acc in protID2geneid):
                vstr += f'-{protID2geneid[s.acc]}-P{c.uid}' #add patternID: pid = f'P{c.uid}'
                #print(vstr)
            #check if vstr contains -a or -b and count
            if("-a" in vstr):
                numa += 1
            if("-b" in vstr):
                numb += 1
            vstr+=f' (0,0,START)'
            flatsdomlist = []
            for t in s.sdomlist:
                for (u,v) in t.subdoms:
                    flatsdomlist.append((u,v,t.fid))
                    #print(t.fid)
                    if not anndone:
                        if(t.fid != "END" and curidx < len(annlist)):
                            correctannlist.append((u,annlist[curidx]))
                            #print(f'u {u} curidx {curidx} ann {annlist[curidx]}')
                            curidx+=1
            sdomlistsorted = sorted(flatsdomlist, key = lambda x: x[0], reverse = False)
            if not anndone:
                correctannlistsorted = sorted(correctannlist, key = lambda x: x[0], reverse = False)
                #print(f'corrected sorted annlist {correctannlistsorted}')
                #collect corrected annotation
                newann = ""
                for (spos,a) in correctannlistsorted:
                    if newann == "":
                        newann = newann+str(a)
                    else:
                        newann = newann+"="+str(a)
                idline = f'>PatternID {c.uid}; #Accessions {len(c.sseqs)}; Structural_annotation {newann}\n'
                #print(idline)
                outf.write(idline)
                anndone=1
            preup = 0 #previous up-value (right border)
            for (u,v,fid) in sdomlistsorted:
                curleft = u
                if(curleft < preup):
                    curleft = preup + 1
                preup = v
                vstr+=f' ({curleft},{v},{fid})'
            outf.write(f'{vstr}\n')

        #print(idline)
        #print(f'{vstr}\n')
    return ((numa,numb))


#======================================end writetofileverboseclus



def readFasta(fastafile):
    #reads the faste files and returns a dict with accession -> sequence length
    pid2length = dict()
    with open(fastafile) as ff:
        curid = ""
        curlen = 0
        for ll in ff:
            if(ll[0] == ">"):
                ls = ll.split(" ")
                if(curid != ""):
                    cs = curid.split('|')
                    tmpcurid = cs[0]
                    pid2length[tmpcurid] = curlen
                    #print(f'PID2LEN:{curid}:{tmpcurid}:{curlen}')
                curid = ls[0][1:]
                curlen=0
            else:
                lls = ll.strip()
                curlen+=len(lls)
        cs = curid.split('|')
        tmpcurid = cs[0]        
        pid2length[tmpcurid] = curlen
        #print(f'PID2LEN:{curid}:{tmpcurid}:{curlen}')

    return pid2length


def sortRanges(propset, rangestr):
    #get the set of properties and the ranges
    #return the annotation, up, low and the list of subdomain ranges
    sansp = propset
    doms = rangestr.split(',')
    if(len(doms) == 1):
        ft = doms[0].split('-')
        low = int(ft[0])
        up = int(ft[1])
        subranges = [(low,up)]
        if(sansp == {''} or sansp == {" "} or sansp == {}): #take single letters, otherwise they won't be disinguishable, e.g LLL and LLLIS
            cann = "L" #Property field is empty
        elif(sansp == {"IS"}):
            cann = "K"#"LIS"
        elif(sansp == {"CP"}):
            cann = "M"#"LCP"
        elif(sansp == {"CP", "IS"}):
            cann = "A"#"LCI"
        else:
            print("Should not appear:")
            print(sansp)
    else: #len(doms) > 1
        ofnum = len(doms)-1
        subranges = []
        annlist = []
        for i in range(len(doms)):
            ft = doms[i].split('-')
            subranges.append((int(ft[0]),int(ft[1])))
            if(sansp == {"NC"}):
                tmpann = "NC"+str(i)+"-"+str(ofnum)
            elif(sansp == {"NC","IS"}):
                tmpann = "NI"+str(i)+"-"+str(ofnum)
            elif(sansp == {"NC","CP"} ):
                tmpann = "NP"+str(i)+"-"+str(ofnum)
            elif(sansp == {"NC","CP", "IS"} ):
                tmpann = "NB"+str(i)+"-"+str(ofnum) #NB instead of NPI
            else:
                print(f'unknown structure: {sansp}')
            annlist.append(tmpann)

        low = subranges[0][0]
        up = subranges[-1][1]
        cann = "=".join(annlist)

    return (cann, low, up, subranges)


#==========================end sortRanges




#===================================main======================================
#input format:
# Accession	E-Value	Residue Range	Property	Architecture	X-group	T-group	F-group	F-id

pid2length = dict()
if(doseqlen):
    pid2length = readFasta(sequencefile)
#create END SDomains for each category which show the average sequence length
#set fid = END and map a certain color (black?) to the end sdomain


protID2geneid = dict()
#protID2line = dict()

if(dotranslate):
    file1 = open(translationfile, 'r')
    inlines = file1.readlines()

    for line in inlines:
        line2 = line.rstrip()
        cols = line2.split("\t")
        protID2geneid[cols[3]] = cols[0]
        #protID2line[cols[3]] = line2

    file1.close()





sdomlist = [] #number of sdomains, however accessions can appear several times
cluslist = [] #number of clusters, however accessions can appear several times
catlist = [] #list of categories regarding clusters of domain sequences
countsdomains = 0

acc2fids = dict()
acc2sdoms = dict()

#print additional outline for COG statistics
#outline:
#number of clusters: 527
#number of SDomains: 3366
#number of categories: 96
#Total Proteins:   1410
#Total Domains:    3555
# NC : 105 (2.95%)
# CP :   1 (0.03%)
# IS :  68 (1.91%)


#format: COGid numprot numdoms L L% NC NC% CP CP% IS IS% sdoms clus cats



pref = inputfile.split("/")
ps = pref[-1].split("_")
cogid = ps[0]
nump = 0
numd = 0
ltot = 0
lperc = ""
nctot = 0
ncperc = ""
cptot = 0
cpperc = ""
istot = 0
isperc = ""



with open(inputfile) as f:
    for line in f:
        curline = line.strip()
        if(curline[0] == '#'):
            nline = curline[1:].strip()
            if(nline.startswith("Total")):
                ns = nline.split()
                nump = int(ns[2])
                numd = int(ns[5])
            if(nline.startswith("NC :")):
                ns = nline.split()
                nctot = int(ns[2])
                ncperc = ns[3]
            if(nline.startswith("CP :")):
                ns = nline.split()
                cptot = int(ns[2])
                cpperc = ns[3]
            if(nline.startswith("IS :")):
                ns = nline.split()
                istot = int(ns[2])
                isperc = ns[3]
            continue
        ls = curline.split('\t')
        #one line = one sdomain

        #FID, we want to use the complete FID and not only parts of it
        fid = ls[-1]
        #prefid = ls[-1]
        #fid = "N.A"
        #if(prefid != "N/A"):
        #    ps = prefid.split('.')
        #    fid = str(ps[0])+"."+str(ps[1])

        #ANNOTATION
        #print(ls)
        ansp = ls[3].split(" ")
        sansp = set(ansp)
        (cann, low, up, subranges) = sortRanges(sansp, ls[2]) #there can be several intervals of domains in one line

        #CREATE SDomain
        acc = ls[0]
        curs = SDomain(cann, acc, low, up, fid, ls[1], subranges, len(subranges), ls[3])
        sdomlist.append(curs)
        countsdomains += 1
        #print(f'curfid {fid} and acc {acc}\n')

        if(acc in acc2fids):
            if(fid in acc2fids[acc]):
                acc2sdoms[acc].append(curs) #several sdomains per accession
            else:
                acc2fids[acc].append(fid)
                acc2sdoms[acc].append(curs)
        else:
            acc2fids[acc] = [fid]
            acc2sdoms[acc] = [curs]
        #first, sort by accession, and then check, if we can create clusters (same fids) or categories from it


#UPDATE (03/2024): create categories from every accession (one protein per category)
#then join two categories if they have the same pattern of sDomains, join into a cluster
# we need to make sure that the sdomains of one sequence stay together

#print(f'number of sdomains: {countsdomains}')

sseqlist = []
maxlen = 800

ssid = 0
for (acc,sdlist) in acc2sdoms.items():
    #print(f'acc and len sdoms: {acc} with {len(sdlist)}\n')
    if(doseqlen):
        lena = maxlen
        tmpas = acc.split('|')
        tmpacc=tmpas[0]
        if(tmpacc in pid2length):
            lena = pid2length[tmpacc]
        else:
            print(f'WARNING: neither {acc} nor {tmpacc} in pid2length, is it in fasta?\n')
        endomain = SDomain("END", acc, lena, lena+10, "END", 0, [(lena,lena+10)], 1, "END")
        sdlist.append(endomain)
    sseq = SSequence(sdlist,ssid)
    ssid+=1
    sseqlist.append(sseq)
    lenfs = len(sseq.getFIDs())
    lfs = sseq.getFIDs()
    #print(f'num fids cursseq: {lenfs} and lfs {lfs}, {sseq.acc}, {sseq.ssid}, {sseq.low}\n')


#create clusters of SSequence if they share the same pattern
scluslist = []
uid=0
#sseqlist_lensorted = sorted(sseqlist, key = lamda x: len(x.sdomlist), reverse = False)
#sseqlist_annsorted = sorted(sseqlist, key = lambda x: x.ann, reverse = False)
sseqlist_fidsorted = sorted(sseqlist, key = lambda x: x.createDIDStrLength(), reverse=False)




#for s in sseqlist_fidsorted:
#    print(s.createDIDStr())

if(len(sseqlist_fidsorted)==0):
    #print(f'WARNING: sseqlist sorted is empty!\n')
    resratio = 0
    if(nump >0):
        resratio = round(len(sseqlist)/nump,2)
    print(f'{cogid}\t{nump}\t{numd}\t{ltot}({lperc}%)\t{nctot}{ncperc}\t{cptot}{cpperc}\t{istot}{isperc}\t0\t{countsdomains}\t{len(sseqlist)}\t{resratio}')
    exit
else:
    i = 0
    cursseq = sseqlist_fidsorted[i]
    #print(len(sseqlist_annsorted))
    cursclus = SClus(cursseq.sdomlist[0].low, cursseq.sdomlist[-1].up, [cursseq], uid, cursseq.ann)
    uid+=1
    #create extra statistics regarding the patterns
    if(dodetails):
        outpat = open(patternfile,"a") #patternstats
        outpat.write("Pattern\tCOGid\tnumA\tnumB\n")
        numa = 0
        numb = 0
        #directly count for cursclus if it is a or b
        #TODO
        if(cursseq.acc in protID2geneid):
           # print(f'accession {cursseq.acc}')
            #print(protID2geneid[cursseq.acc])
            if("a" in protID2geneid[cursseq.acc]):
                numa+=1
            elif("b" in protID2geneid[cursseq.acc]):
                numb+=1
            else:
                print(f'WARNING counts1: {protID2geneid[cursseq.acc]} of {cursseq.acc}')
        curfidstr = cursseq.createDIDStr()
        #print(f'curfid str {curfidstr} ')
    while(i < len(sseqlist_fidsorted)-1):
        cursseq = sseqlist_fidsorted[i]
        for j in range(i+1,len(sseqlist_fidsorted)):
            jsseq = sseqlist_fidsorted[j]
            added = False
#            print(f'i:{i} idstr {cursseq.createDIDStr()} j:{j} idstr {jsseq.createDIDStr()}')
            if(cursseq.createDIDStr() == jsseq.createDIDStr()):
                #check if they can be in the same Sclus
                if(cursclus.addSSeq(jsseq)):
                    added = True
                    if(dodetails):
                        if(jsseq.acc in protID2geneid):
                            #print(f'accession jsseq {jsseq.acc}')
                            #print(protID2geneid[jsseq.acc])
                            if("a" in protID2geneid[jsseq.acc]):
                                numa+=1
                            elif("b" in protID2geneid[jsseq.acc]):
                                numb+=1
                            else:
                                print(f'WARNING counts2: {protID2geneid[cursseq.acc]} of {cursseq.acc}')
                                #counts2: a272557 of BAA80813.1

                    i+=1
                else:
                    scluslist.append(cursclus)
                    #new cluster
                    if(dodetails):
                        outpat.write(f'{curfidstr}\t{cogid}\t{numa}\t{numb}\n')
                        numa = 0
                        numb = 0
                        if(jsseq.acc in protID2geneid):
                            if("a" in protID2geneid[jsseq.acc]):
                                numa+=1
                            elif("b" in protID2geneid[jsseq.acc]):
                                numb+=1
                            else:
                                print(f'WARNING counts3: {protID2geneid[cursseq.acc]} of {cursseq.acc}')
                        curfidstr = jsseq.createDIDStr()
                    cursclus = SClus(jsseq.sdomlist[0].low, jsseq.sdomlist[-1].up, [jsseq], uid, jsseq.ann)
                    
                    uid+=1
                    i+=1
            else:
                #append clus and stop inner loop to restart with a new cluster
                if(added == False):
                    scluslist.append(cursclus)
                    #new cluster
                    if(dodetails):
                        outpat.write(f'{curfidstr}\t{cogid}\t{numa}\t{numb}\n')
                        numa = 0
                        numb = 0
                        if(jsseq.acc in protID2geneid):
                            if("a" in protID2geneid[jsseq.acc]):
                                numa+=1
                            elif("b" in protID2geneid[jsseq.acc]):
                                numb+=1
                            else:
                                print(f'WARNING counts4: {protID2geneid[cursseq.acc]} of {cursseq.acc}')
                        curfidstr = jsseq.createDIDStr()
                    cursclus = SClus(jsseq.sdomlist[0].low, jsseq.sdomlist[-1].up, [jsseq], uid, jsseq.ann)
                    uid+=1
                    i+=1
                    break
    scluslist.append(cursclus)
    if(dodetails):
        #check for numa and numb!

        #write to outfile
        outpat.write(f'{curfidstr}\t{cogid}\t{numa}\t{numb}\n')
    
    #print(f'number of Sclusters: {len(scluslist)}')
    #print(f'number of SSequences: {len(sseqlist)}')


    #format: COGid numprot numdoms L L% NC NC% CP CP% IS IS% sdoms clus cats



    #prepare outline:
    ltot = numd - nctot - cptot - istot
    prelperc = 0
    if(numd>0):
        prelperc = round(ltot/numd * 100,2)
    lperc = str(prelperc)
    resratio = 0
    if(nump >0):
        resratio = round(len(sseqlist)/nump,2)

    print(f'{cogid}\t{nump}\t{numd}\t{ltot}({lperc}%)\t{nctot}{ncperc}\t{cptot}{cpperc}\t{istot}{isperc}\t{len(scluslist)}\t{countsdomains}\t{len(sseqlist)}\t{resratio}\t{numa}\t{numb}')



    #NEXT
    #add some parameters to the plot in order to dynamically set the plot size
    #output is the number of a or b proteins
    (sseqa,sseqb) = writeSClustofileverbose(scluslist,outputname,pid2length,protID2geneid)

    #check for numbers

    if(numa != sseqa or numb != sseqb):
        print(f'WARNING: counts for archaeal {numa},{sseqa} and bacterial {numb},{sseqb} do not match!')

    #sdomlist = [] #number of sdomains, however accessions can appear several times
#cluslist = [] #number of clusters, however accessions can appear several times
