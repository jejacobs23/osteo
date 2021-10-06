import Reactome_Tools2
import Process_Pathway_Analysis
import os
import sys

COM_DIR = "/home/exacloud/lustre1/jjacobs"
all_runs = ["01101_M1", "01101_M2", "01101_M3", "01101_M4", "01101_M5", "01101_M6", \
            "01101_M7", "01101_M8", "01103_D1", "01104_M1", "01105_D1", "01106_D1", \
            "01107_M1", "01107_M2", "01108_M1", "01109_D1", "01110_D1", "01111_M1", \
            "01112_D1", "01112_M1", "01112_M2", "01115_D1", "01117_D1", "01118_D1", \
            "01119_D1", "01120_D1", "01123_D1", "01124_D1", "01125_D1", "01126_D1", \
            "01127_D1", "01128_D1", "01131_D1", "01_M", "02_D", "03_D", "04_D", "05_D", \
            "06_D", "07_D", "08_D", "09_D", "10_D", "10_M", "11_D", "11_M", "12_D", \
            "12_M", "13_D", "13_M", "14_D", "15_D", "16_D", "17_D", "18_D", "19_D"]

def get_survival_data_Sample_IDs():
    L = []
    COM_DIR = "/home/exacloud/lustre1/jjacobs"
    infile = COM_DIR + "/data/osteo/Survival_Info.txt"
    fi = open(infile, 'r')
    next(fi)
    for i in fi:
        i = i.split()
        sID = i[0][:-1]#removes the last character from the subject ID (D or M)
        if sID not in L:
            L.append(sID)
    fi.close()
    print(L)
    return(L)

#Function to looks at a pathway and determins (using D) what level of separation exists between the
#different groups.  In the case of the osteo samples, this will be separation between patients with
#aberrations in that pathway vs patients without aberrations in that pathway.  There are a subset of
#samples that have survival data so only those will be considered.  However, any clinical group can
#be submitted for sepValue assessment by adjusting the GroupList variable.
def sepValue(D, Pathway, GroupList):
    print("Yokl!  Let's make some sepValues.")
    print(GroupList)
    print(len(GroupList))
    counter = 0
    dups = set()
    for j in D[Pathway]:
        Sample_Name = j[0]
        patID = Sample_Name[:-1] #removes the last character from the subject ID (D or M)
        if patID in GroupList:
            if patID not in dups:
                counter += 1
                dups.add(patID)
    sep = len(GroupList) - counter
    return(sep)

def daddy(G, Pathway):
    print("Looking at pathway " + Pathway)
    parents = []
    for i in G:
        if Pathway in G[i]:
            if i not in parents:
                parents.append(i)
    return(parents)

def checkSep(sepD, Pathway, upperCutOff, lowerCutOff):
    sep = sepD[Pathway]
    if sep >= lowerCutOff:
        if sep <= upperCutOff:
            return("winner")
        else:
            return("dud")
    else:
        return("dud")

def checkNextStep(pathwaysKeep, Roots, G, sepD, famEval):
    L = []
    for i in pathwaysKeep:
        if i not in Roots:
            rents = daddy(G=G, Pathway=i)
            for r in rents:
                if r in sepD:
                    fam = r + ":" + i
                    if fam not in famEval:
                        L.append(r)
                    else:
                        print("Busted!!  This parent:child combo has already been evaluted. " + fam)
    return(L)

def separation(D, G, upperCutOff, lowerCutOff):
    pathEval = set() #To keep track of individual pathways that have already been evaluated
    famEval = set() #To keep track of pathway family units that have already been evaluated
    sepD = {}#This will contain the separation info for each pathway in D
    pathwaysKeep = set()
    RL = Reactome_Tools2.RootsLeaves(G)
    Roots = RL[0]
    Leaves = RL[1]
    for i in D:
        sepD[i] = sepValue(D=D, Pathway=i, GroupList=GroupList)
    for f in Leaves:
        if f in sepD:
            if checkSep(sepD=sepD, Pathway=f, upperCutOff=upperCutOff, lowerCutOff=lowerCutOff) == "winner":
                pathwaysKeep.add(f)
    next = checkNextStep(pathwaysKeep=pathwaysKeep, Roots=Roots, G=G, sepD=sepD, famEval=famEval)
    while len(next) > 0:
        stuffToAdd = set()
        stuffToSubtract = set()
        for k in pathwaysKeep:
            if k not in Roots:
                pathEval.add(k)
                rents = daddy(G=G, Pathway=k)
                for r in rents:
                    if r in sepD:
                        fam = r + ":" + k
                        if fam not in famEval:
                            famEval.add(fam)
                            if checkSep(sepD=sepD, Pathway=r, upperCutOff=upperCutOff, lowerCutOff=lowerCutOff) == "winner":
                                stuffToSubtract.add(k)
                                stuffToAdd.add(r)
                        else:
                            print("DangerMouse!!  This parent:child combo has already been evaluted. " + fam)
        pathwaysKeep = pathwaysKeep.difference(stuffToSubtract)
        pathwaysKeep.update(stuffToAdd)
        next=checkNextStep(pathwaysKeep=pathwaysKeep, Roots=Roots, G=G, sepD=sepD, famEval=famEval)
    return(pathwaysKeep)

INFILE1=sys.argv[1]
OUTFILE1=sys.argv[2]
INFILE2=sys.argv[3]
OUTFILE2=sys.argv[4]
D = {}
for i in all_runs:
    ALIGNMENT_RUN = "SJOS0" + i
    print("Currently analyzing " + ALIGNMENT_RUN)
    INFILE = COM_DIR + "/data/osteo/" + ALIGNMENT_RUN + "/" + INFILE1
    if os.path.isfile(INFILE):
        D = Reactome_Tools2.read_results(Sample_Name=ALIGNMENT_RUN, infile=INFILE, D=D)
    else:
        print("What the fork! : " + INFILE)
S = Reactome_Tools2.subSamples(runs=all_runs)
Condensed = Reactome_Tools2.condense_D(D=D, S=S)
GroupList = get_survival_data_Sample_IDs()
upperCutOff = len(GroupList) - 10
lowerCutOff = 10
print(upperCutOff)
print(lowerCutOff)
Results = Reactome_Tools2.pathway_relations()
G = Results[0]
BestPathways = separation(D=Condensed, G=G, upperCutOff=upperCutOff, lowerCutOff=lowerCutOff)
outfile = COM_DIR + "/data/osteo/Reactome_Analyses/Pre_Filtered/" + OUTFILE1
with open(outfile, 'w') as fo:
    for p in BestPathways:
        fo.write(p + '\n')
infile = COM_DIR + "/data/osteo/Reactome_Analyses/Pre_Filtered/" + INFILE2
outfileP = COM_DIR + "/data/osteo/Reactome_Analyses/Pre_Filtered/" + OUTFILE2

Process_Pathway_Analysis.ProPath_LL(D=Condensed, infile=infile, outfile=outfileP)
