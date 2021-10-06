#Function to split the lines of a "results" file.  The reason this is neccessary is because there are
#commas inside quotes so splitting by "," doesn't work.
def special_line_split(s):
    L = []
    Q = []
    switch = "ON"
    for i in range(0,len(s)):
        if s[i] =='"':
            if switch == "ON":
                switch = "OFF"
            elif switch == "OFF":
                switch = "ON"
        elif s[i] == ',' and switch == "ON":
            L.append(i)
    L.append(len(s))
    u = 0
    for x in L:
        Q.append(s[u:x])
        u = x + 1
    return(Q)

#Function to read the "results" output from a Reactome analysis run.  This is a comma-separated-file.
#
#File Header Reference:
#      Pathway identifier,Pathway name,#Entities found,#Entities total,#Interactors found,#Interactors total,
#      Entities ratio,Entities pValue,Entities FDR,#Reactions found,#Reactions total,Reactions ratio,Species identifier,
#      Species name,Submitted entities found,Mapped entities,Submitted entities hit interactor,Interacts with,Found reaction identifiers
#
def read_results(Sample_Name, infile, D):
    fi = open(infile, 'r')
    print("Reading from " + infile)
    for i in fi:
        if i[0] == "R":
            Q = special_line_split(s=i)
            if len(Q) != 19:
                print("Something is fishy here " + Q[0] + " " + str(len(Q)))
                print(i)
            pID = Q[0]
            PName = Q[1]
            NumEF = Q[2]
            NumET = Q[3]
            InteractF = Q[4]
            InteractT = Q[5]
            Eratio = Q[6]
            EpVal = Q[7]
            FDR = Q[8]
            RF = Q[9]
            RT = Q[10]
            Rratio = Q[11]
            SpeciesID = Q[12]
            SpeciesN = Q[13]
            SubmitedE = Q[14].strip('"')
            GeneList = SubmitedE.split(";")
            if GeneList == ['']:
                GeneList = []
            GenesF = len(GeneList)
            MappedE = Q[15]
            SubmittedEhit = Q[16]
            InteractsWith = Q[17]
            FoundReactionIDs = Q[18]
            if pID in D:
                D[pID].append([Sample_Name, GeneList, GenesF, NumEF, NumET, InteractF, InteractT])
            else:
                D[pID] = [[Sample_Name, GeneList, GenesF, NumEF, NumET, InteractF, InteractT]]
    fi.close()
    return(D)

def subSamples(runs):
    S = {} #This will contain the baseline length for each uniqe sample ID (i.e. how many subsamples fall under that ID)
    for r in runs:
        r = "SJOS0" + r
        r = r.split("_")
        sp = r[0]
        sx = r[1][0]
        sID = sp + sx
        if sID in S:
            S[sID] = S[sID] + 1
        else:
            S[sID] = 1
    return(S)

#Function to take the final dictionary (D) with all the samples and condensing it to include only individual patients at specific timepoints.
#i.e. if there are multiple samples taken at one timepoint, it will condense those values into one entry.
def condense_D(D, S):
    E = {} #This dictionary will eventually contain the condensed version of D
    for i in D: #i interates through each pathway
        C = {} #C will contain the different types of samples for that individual: idealy, some number of diagnostic (D) and some number of metastatic (M) samples
        V = {} #V will be the condensed form of C
        for s in D[i]: #Iterate through the various sample ID's that have mutations in that particular pathway, i.
            SN = s[0]
            SN = SN.split("_")
            p = SN[0] #This is the number assigned to each subject
            x = SN[1][0] #This is the subclassification (M or D)
            ID = p + x
            if ID in C:
                C[ID].append([s[0], s[1], s[2], s[3], s[4], s[5], s[6]])
            else:
                C[ID] = [[s[0], s[1], s[2], s[3], s[4], s[5], s[6]]]
        for j in C: #Iterate through the unique sample ID's that have mutations in pathway i.  C resets for each new pathway in D
            G = {} #This will contain the gene name as the key and the number of samples that gene is mutated in as the value.  
            aveNumEF = 0.0
            aveInteractF = 0.0
            aveGenesF = 0.0
            if S[j] > 1: #This identifyies uniqe sample ID's (j) that have multiple subsample ID (such as 1101_M1 and 1101_M2).
                GeneList = [] #initialize an empty gene list
                for u in C[j]: #Iterate through all the subsamples in j
                    aveNumEF = aveNumEF + float(u[3])
                    aveInteractF = aveInteractF + float(u[5])
                    for g in u[1]:
                        if g not in G:
                            G[g] = 1
                        else:
                            G[g] = G[g] + 1
                        if g not in GeneList:
                            GeneList.append(g)
                for h in G:
                    aveGenesF = aveGenesF + (1.0 * float(G[h]))/float(S[j]) #This will account for genes that don't overlap with other subsamples
                aveNumEF = float(aveNumEF)/float(S[j])
                aveInteractF = float(aveInteractF)/float(S[j])
            else:
                aveNumEF = float(C[j][0][3])
                aveInteractF = float(C[j][0][5])
                aveGenesF = float(len(C[j][0][1]))
                GeneList = C[j][0][1]
            V[j] = [GeneList, str(aveGenesF), str(aveNumEF), C[j][0][4], str(aveInteractF), C[j][0][6]]
        E[i] = []
        for L in V:
            E[i].append([L, V[L][0], V[L][1], V[L][2], V[L][3], V[L][4], V[L][5]])
    return(E)


#Function to write a new .txt file from the Reactome output dictionary
def write_txt_from_dict(D, outfile):
    import os
    with open(outfile, 'w') as fo:
        fo.write('\t'.join(["Reaction_ID", "#Samples", "Sample_ID's", "Mutation_Counts", "Interactor_Counts", "Total_Entities", "Total_Interactors", "Gene_Counts", "Genes"]) + '\n')
        for u in D:
            S = ""
            C = ""
            R = ""
            G = ""
            GC = ""
            TC = ""
            TR = ""
            for i in D[u]:
                S = S + i[0] + ","
                C = C + str(i[3]) + ","
                R = R + str(i[5]) + ","
                genes = ""
                for g in i[1]:
                    genes = genes + g + ";"
                genes = genes.rstrip(";")
                G = G + genes + ","
                GC = GC + str(i[2]) + ","
                if len(TC) > 0:
                    if str(i[4]) != TC:
                        print("Oops! Different number of entities in different samples : " + u + ";" + i)
                else:
                    TC = TC + str(i[4])
                if len(TR) > 0:
                    if str(i[6]) != TR:
                        print("Oops! Different number of interactors in different samples : " + u + ";" + i)
                else:
                    TR = TR + str(i[6])
            S = S.rstrip(",")
            C = C.rstrip(",")
            R = R.rstrip(",")
            G = G.rstrip(",")
            GC = GC.rstrip(",")
            fo.write('\t'.join([u, str(len(D[u])), S, C, R, TC, TR, GC, G]) + '\n')
    sorted = outfile.rstrip(".txt") + "_sorted.txt"
    command = "(head -n 1 " + outfile + " && tail -n +2 " + outfile + " |sort -k 2nr) > " + sorted
    os.system(command)

#function to create a dictionary representing the hierarchy of the Reactome pathways.  Each node is 
#represented as a key and all other nodes its connected to are the associated values.  This will
#store all the information needed to represent a hierarchical tree of pathways.  Of note, the 
#lowest level pathways are not included in the "ReactomePathwayRelation.txt" file that was downloaded
#from the Reactome website. The second dictionary, V, contains all the needed information about a 
#specific pathway.  It holds a list of lists, the later of which contains a counter of how many 
#mutations exist in that pathway and the sample those mutations came from.  Each sample will have
#it's own count so they can be differentiated during analysis.
#
#Updated version of pathway relationship file downloaded from Reactome on 6/21/20.   
#    
def pathway_relations():
    infile = "/home/exacloud/lustre1/jjacobs/Reactome_Downloads/ReactomePathwaysRelation_06_21_20.txt"
    with open(infile, 'r') as fi:
        G = {}
        V = {}
        for i in fi:
            i = i.split()
            parent = i[0]
            child = i[1]
            if parent in G:
                G[parent].append(child)
            else:
                G[parent] = [child]
            if parent not in V:
                V[parent] = []
            if child not in V:
                V[child] = []
    return([G, V])

#Function to find the preceeding node or nodes.  This will return a list just in case
#there are child nodes that have more than one parent.
def preNode(G, curNode):
    preN = []
    for i in G:
        for u in G[i]:
            if u == curNode:
                preN.append(i)
    if len(preN) > 1:
        print("This is odd.  There is a child node with more than 1 parent :" + curNode)
    return(preN)

#Function to move one level down in the pathway hierarchy.  Accespt multiple parent nodes and returns
#all child nodes of each parent.  Since a child node can have multiple parents, if L is a list instead
#of a set(), it could potentially contain duplicate child nodes.  
def NextLevel(G, nodes):
    L = set() #This set will hold all of child nodes of all the parents (without duplicates)
    for i in nodes:
        if i in G: #This will identify nodes that have child nodes
            for b in G[i]:
                L.add(b)
    return(L)

#Function to update the mutation counts of the parent and child nodes in a given hierarchy tree.  Evaluates all
#the child nodes of each parent.  Three interesting quirks in the pathway hierarchy: 1) A child node can have multiple parent
#nodes, and 2) a node can be evaluated as a parent multiple times at different "levels".  This requires keeping track
#of the child:parent combos as you loop through the nodes so as to not subtract mutations in a child node more
#than once from its parent.  and finally 3) There can be overlaping child pathways.  That is to say child pathways that
#share parts of their parent pathway in common.  Thus, simply subtracting the raw entities counts won't work.  You have
#to subtract the mutated gene counts from the partent pathway for only those genes that were mutated in the child.
def UpdateMutationCounts(P, PV, G, V, D, nodes, S):
    L = NextLevel(G=G, nodes=nodes) #remember that L is a set(), not a list.
    print("Let's update some mutation counts" + '\n')
    for i in nodes:
        if i in G: #This will identify nodes that have child nodes
            for u in G[i]: #Iterate through the child nodes in order to add the mutation counts
                fam = u + ":" + i #This marks a specific child:parent combo since there are child nodes with multiple parent nodes
                print("New child pathway : " + fam)
                if u in D and u not in P: #This child pathway has mutations and hasn't been evaluated yet.
                    P.add(u)
                    for s in D[u]:
                        V[u].append([s[1], s[2], s[3], s[0], s[4], s[5], s[6]]) #This adds the mutated gene lists, ave_mutated_gene_counts, mutation counts, sample names, total entities, interactor mutation count and total interactos to the information for that pathway in V
                    else:
                        print("There are multiple occurences of this child pathway: " + u + " , but that's ok.")
                for q in V[i]: #Iterate through the samples that have mutations in the parent pathway
                    print("New Sample in the current parent pathway : " + q[3])
                    print("This is value of q[1]")
                    print(q[1])
                    STAMP = float(q[1]) #Need to lock in the value of mutated genes found as it will be updated in each iteration
                    child_overlapers = set() #For each sample in the parent pathway, we need to keep track of mutated genes as they may overlap between multiple child pathways.
                    if len(V[u]) > 0: #No sense in counting mutations in pathways that have none
                        for w in V[u]: #Iterate through the samples that have mutations in the child pathway
                            fam_Sample = fam + ";" + w[3] + ";" + q[3] #This marks the child:parent:ChildSample:ParentSample quartet
                            if fam_Sample not in PV: #This child:parent combo has not yet been evaluated in the context of these samples
                                PV.add(fam_Sample)
                                for g in w[0]: #These are the genes in the child pathway that have mutations in this particular sample.
                                    if w[3] == q[3] and g in q[0] and g not in child_overlapers: #Genes -- thus far not counted -- that are mutated in both parent and child from the same sample
                                        child_overlapers.add(g)
                                        aveOG = STAMP/float(len(q[0])) #Use the geneF from the parent pathway as this reflecs the value of this mutaiton in the parent
                                        q[1] = float(q[1]) - aveOG #Update the mutation counts of the parent pathway in that sample.
                                        if i == "R-HSA-109582":
                                            print("Donkeyballs!!  Here is R-HSA-109582 as we update the mutation counts")
                                            print(q[1])
                                            print(w[1])
                                        if q[1] < 0:
                                            print("Monkeyballs!! The parent node has a negative value for its mutated gene counts : " + i + ";" + str(q[1]) + " Child = " + u + "; sample = " + w[3])
                                            print(q[0])
                                            print(w[0])
                                            print(g)
                                            print(q[3])
                                            print(child_overlapers)
        if i == "R-HSA-109582":
            print("FunTimes!  Here's R-HSA-109582 at the end of UpdateMutationCounts")
            test = 0.0
            for y in V["R-HSA-109582"]:
                test = test + float(y[1])
            print("FunTimes!  Here's the total genesF for R-HSA-109582: " + str(test))
    return[P, PV, V, L]

#Function to identify roots and leaves within a hierarchy
def RootsLeaves(G):
    Roots = set() #Set of highest level pathways that have no incoming edges
    Leaves = set() #Set of lowest level pathways that have no outgoing edges
    for i in G:
        marker = "t"
        for u in G[i]:
            if u not in G: #This will identify terminal nodes as they will not be in the list of parent nodes
                Leaves.add(u)
        for j in G:
            for w in G[j]:
                if w == i:
                    marker = "f" #If any of the child pathways are the same as i, then i can't be a root
        if marker == "t":
            Roots.add(i)
    return[Roots, Leaves]

#Function to identify the parent node, given its child nodes
def FindParent(G, children):
    parents = set()
    for i in G:
        C = set()
        for j in G[i]:
            C.add(j)
        for x in children:
            if x in C:
                parents.add(i)
    return(parents)

#Function to identify the lowest level pathways that account for all (or most) of the submitted gene names
#in a Reactome analysis.
def LowestLevel(D, G, ALL_RUNS):
    LL = {}
    Pathways = set()
    RL = RootsLeaves(G)
    Leaves = RL[1]
    for r in ALL_RUNS:
        r = "SJOS0" + r
        print("Currently processing " + r)
        Keep_P = set()
        Keep_G = set()
        Exclude_P = set()
        Exclude_G = set()
        for p in D:
            for s in D[p]:
                if (r == s[0]) and (len(s[1]) > 0): #Identifies mutations in sample r from pathway p
                    print("Sucka!  This is goin' down! " + p)
                    Keep_P.add(p)
                    for g in s[1]:
                        Exclude_G.add(g)
        for i in Keep_P:
            if i not in Leaves: #if it's not already the lowest possible level
                if i in G:
                    children = G[i]
                    for c in children:
                        if c in D:
                            for xs in D[c]:
                                if (r == xs[0]) and (len(xs[1]) > 0): #The child pathway also has mutations in it from the sample, r.
                                    Exclude_P.add(i)
                else:
                    print("This is probably a disease pathway: " + i)
        Keep_P = Keep_P - Exclude_P
        for k in Keep_P:
            for ks in D[k]:
                if r == ks[0]:
                    for kg in ks[1]:
                        Keep_G.add(kg)
        Exclude_G = Exclude_G - Keep_G
        LL[r] = [Keep_P, Exclude_P, Keep_G, Exclude_G]
        total_P = len(Keep_P) + len(Exclude_P)
        percent_G = float(len(Keep_G))/(float(len(Keep_G)) + float(len(Exclude_G)))
        print("Number of kept pathways for " + r + ": " + str(len(Keep_P)) + " out of " + str(total_P) + " total")
        print("Accounting for " + str(percent_G) + " of the genes.")
    print("Well, we got this far.")
    print("This should be the number of samples: " + str(len(LL)))
    for x in LL:
        evaluated_p = set()
        evaluated_c = set()
        print("Currently attempting to improve " + x)
        print("Number of excluded genes: " + str(len(LL[x][3])))
        print("Number of kept genes: " + str(len(LL[x][2])))
        while len(LL[x][3]) > 0:
            G_current = len(LL[x][3])
            Add = set()
            Subtract = set()
            for child in LL[x][0]:
                if child not in evaluated_c:
                    parents = FindParent(G=G, children=set([child]))
                    for y in parents:
                        if y not in evaluated_p:
                            evaluated_p.add(y)
                            print("Child:Parent = " + child + ":" + y)
                            if y in D:
                                for z in D[y]:
                                    if x == z[0]:
                                        current = len(LL[x][3])
                                        for t in z[1]:
                                            LL[x][2].add(t)
                                            LL[x][3].discard(t)
                                        updated = len(LL[x][3])
                                        if updated < current: #including this parent adds genes to the Kept_G set
                                            print("Moving up one level improved things!")
                                            print("Number of excluded genes: " + str(len(LL[x][3])))
                                            print("Number of kept genes: " + str(len(LL[x][2])))
                                            Add.add(y)
                                            Subtract.add(child)
                                            LL[x][1].discard(y)
                                            LL[x][1].add(child)
                                            for branch in G[y]:
                                                Subtract.add(branch)
                                                LL[x][1].add(branch)
                                        else:
                                            evaluated_c.add(child) #don't need to evaluate it multiple times
            LL[x][0].update(Add)
            LL[x][0] = LL[x][0] - Subtract
            G_updated = len(LL[x][3])
            if G_updated == G_current:
                print("Went through a full loop without any benifit :(")
                break
    for f in LL:
        Pathways.update(LL[f][0])
    return(Pathways)
                    

#Function to iterate through the pathway hierarchy dictionary created from "pathway_relations"
#and populate the V dictinary with information about how many genes in that pathway contain mutations
#and which samples those mutations came from.
def pathway_aborrations(G, V, D, S):
    P = set() #Keep track of pathways that have already been traversed so that we don't count their mutations multiple times
    PV = set() #Keep track of pathways that have already been traversed so that we don't count their mutations multiple times
    Results = RootsLeaves(G=G)
    Roots = Results[0]
    Leaves = Results[1]
    for r in Roots:
        P.add(r) #Add the pathway to the set of transversed pathways
        PV.add(r)
        if r in D: #This will identify pathways that have mutations in them
            for s in D[r]:
                V[r].append([s[1], s[2], s[3], s[0], s[4], s[5], s[6]]) #This adds the mutated gene lists, ave_mutated_gene_counts, mutation counts, sample names, total entities, interactor mutations counts and total interactors to the information for that pathway in V
    nodes = Roots
    while len(nodes) > 0:
        O = UpdateMutationCounts(P=P, PV=PV, G=G, V=V, D=D, nodes=nodes, S=S)
        P = O[0]
        PV = O[1]
        V = O[2]
        nodes = O[3]
    print("Xena!  Here's the V entry for R-HSA-109582")
    print(V["R-HSA-109582"])
    test = 0.0
    for y in V["R-HSA-109582"]:
        test = test + float(y[1])
    print("Xena!  Here's the total genesF for R-HSA-109582: " + str(test))
    return(V)

#Function to evaluate the pathways in a given hierarchical level and exclude those that don't have sufficient
#mutational burden.
def levelBurden(level, D, CO):
    print("Friday! Here is the winning level")
    winner = []
    Dcount = []
    Pexcluded = set()
    for p in level:
        if p not in D:
            Pexcluded.add(p)
            print("This pathway was excluded from this level because it's not in D: " + p)
            continue #This will skip the rest of the lines in this iteration of the for loop and return to the top of the loop for the next interation.
        Scount = 0.0
        winner.append(p)
        for i in D[p]:
            Scount = Scount + float(i[2])
        Dcount.append(str(Scount))
    level = level - Pexcluded #gets rid of pathway ID's that are not in D.  This happens sometimes with the infectious disease pathways
    print('\t'.join(winner))
    print('\t'.join(Dcount))
    total = 0.0
    Dcount = [float(x) for x in Dcount]
    for j in Dcount:
        total = total + j
    lco = CO * total
    ordered = sorted(Dcount, reverse = True) #This will sort the genesF values from highest to lowest
    idx = []
    for n in ordered:
        dex = Dcount.index(n)
        idx.append(dex) #This will create a list with the index values of the pathway list if it were sorted
    sls = 0.0
    position = 0
    while sls < lco: #add the ordered values of genesF until you get to the cutoff value of 95% of the total count for that level
        sls = sls + float(ordered[position])
        position = position + 1
    if len(winner) > position: #identifies levels in which some pathways will be excluded
        excluded = idx[position:]
        Pexcluded = set()
        for x in excluded:
            Pexcluded.add(winner[x])
        for y in Pexcluded:
            print("This pathway was excluded from this level: " + y)
        pathway = level - Pexcluded
    else:
        pathway = level
    return(pathway)

#Function to analyze the mutation data in the pathway hierarchy dictionary (G) and corresponding values dictionary (V)
#and determin the most distal branch point that contains a certain percentage of the total mutations in that 
#particular tree.  The percentage is provided by the user.  For now, I'm not considering interactors in this function.  
def branchPoint(G, V, D, CutOff):
    Pathways = [] #This list will contain the pathways that are proximal to the mutation cutoffs
    P = set() #This set will contain the traversed pathways.  It can be used to see if a pathway is traversed more than once.
    Results = RootsLeaves(G=G)
    Roots = Results[0]
    Leaves = Results[1]
    for r in Roots:
        Levels = [] #This list will contain a list of nodes at each level
        Numbers = [] #This list will contain the mutated gene count at each level
        nodes = [r]
        while len(nodes) > 0:
            Levels.append(nodes)#The first time around, this will just add in the root as a list
            count = 0.0
            for i in nodes:
                for j in V[i]:
                    count = count + float(j[1])
            Numbers.append(count)
            nodes = list(NextLevel(G=G, nodes=nodes))
        print(Levels)
        print(Numbers)
        Total = 0.0
        for t in Numbers:
            Total = Total + t
        s = Total * CutOff
        if s != 0:
            A = 0#This will count the mutations starting from the bottom level and continuing up
            l = (len(Levels) - 1)#Start at the last level
            while A <= s:
                A = A + Numbers[l]
                l = l - 1
            if l != 0:
                pathway = FindParent(G=G, children=Levels[l])
            else:
                print("Looks like we made it all the way back to the root: " + Levels[l][0])
                print(len(Levels[l]))
                pathway = set()
                pathway.add(Levels[l][0])
            pathway = levelBurden(level=pathway, D=D, CO=0.95)
            for z in pathway: #Remember that FindParent returns a set
                if z not in P:
                    P.add(z)
                    if z in D:
                        Pathways.append(z)
                else:
                    print("We have duplicate convergence points! " + z)
                    print("Here are the nodes in the level just below the parent")
                    print(Levels[l])
    return(Pathways)
