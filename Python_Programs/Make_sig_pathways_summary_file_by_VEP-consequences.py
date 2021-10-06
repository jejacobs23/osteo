
#Python program to take the pathway FDR's and create a summary file of the
#significant pathways
#
#P holds the pathway information in a dictionary: P[pID] = [Pnamd, aberrant, non-aberrant, p_val, FDR, Subject_List, Gene_List]

CO = "0.3"

def compile_pathway_data(P, COM_DIR, CO):
    import requests
    P_dict = P
    #INFILE1 = COM_DIR + "/data/osteo/Reactome_Analyses/Pre_Filtered/Cut_Off_" + CO + "/Pathway_p-values_by_VEP-consequences.txt"
    INFILE1 = COM_DIR + "/data/osteo/Reactome_Analyses/Pre_Filtered/Pathway_p-values_by_VEP-consequences.txt"
    infile2 = COM_DIR + "/data/osteo/Reactome_Analyses/Pre_Filtered/Reactome_Pathways_by_Sample_VEP_prefiltered_condensed_by_VEP-consequences_sorted.txt"
    fi1 = open(INFILE1, 'r')
    next(fi1)
    for i in fi1:
        i = i.split()
        pID = i[0]
        r = requests.get("https://reactome.org/ContentService/data/discover/" + pID) #retrieve the pathway name from Reactome API
        Pname = r.json()["name"]
        Pname = Pname.replace(" ", "_")
        print(r.url)
        print(r.status_code)
        aberrant = i[1]
        non_aberrant = i[2]
        p_val = i[3]
        FDR = i[4]  
        if FDR != "NA":
            if float(FDR) < 0.06:
                P[pID] = [Pname, aberrant, non_aberrant, p_val, FDR]
                P_dict = add_subject_IDs_genes(P=P, L=L, Pathway=pID, infile=infile2)
    fi1.close()
    return(P_dict)

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

def screen_Sample_IDs(Sample_ID, L):
    Patient_ID = Sample_ID[:-1]
    print(Patient_ID)
    if Patient_ID in L:
        screen = "Yes"
    else:
        screen = "No"
    print(screen)
    return(screen)

def pull_data(Infile, Pathway):
    D = {}
    fi = open(Infile, 'r')
    next(fi)
    for i in fi:
        i = i.split()
        Reaction_ID = i[0]
        Sample_IDs = i[2]
        if len(i) == 9:
            Genes = i[8]
            if Reaction_ID == Pathway:
                Sample_IDs = Sample_IDs.split(",")
                Genes = Genes.split(",")
                for idx in range(0,len(Sample_IDs)):
                    if len(Genes) > idx:
                        gene_list = Genes[idx].split(";")
                        D[Sample_IDs[idx]] = gene_list
                    else:
                        D[Sample_IDs[idx]] = ["NA"]
        else:
            Genes = "Likely only interactors"
            if Reaction_ID == Pathway:
                Sample_IDs = Sample_IDs.split(",")
                for idx in range(0,len(Sample_IDs)):
                    D[Sample_IDs[idx]] = Genes
    fi.close()
    return(D)

def add_result(D, L, P, Pathway):
    print("Adding stuff")
    SL = ""
    GL = ""
    for i in D:
        screen = screen_Sample_IDs(Sample_ID=i, L=L)
        if screen == "Yes":
            SL += i + ","
            for j in D[i]:
                if j == "":
                    GL += "NA;"
                else:
                    GL += j + ";"
            GL = GL.rstrip(";")
            GL = GL + ","
    SL = SL.rstrip(",")
    GL = GL.rstrip(",")
    print(SL + '\t' + GL)
    P[Pathway].extend([SL, GL])
    return(P)

def add_subject_IDs_genes(P, L, Pathway, infile):
    D = pull_data(Infile=infile, Pathway=Pathway)
    P_dict = add_result(D=D, L=L, P=P, Pathway=Pathway)
    return(P_dict)

def write_output(P, COM_DIR, CO, Notes):
    import os 
    #P[pID] = [Pname, aberrant, non_aberrant, p_val, FDR, Subject_List, Gene_List]
    #OUTFILE = COM_DIR + "/data/osteo/Reactome_Analyses/Pre_Filtered/Cut_Off_" + CO + "/Significant_Pathways_by_VEP-consequences.txt"
    OUTFILE = COM_DIR + "/data/osteo/Reactome_Analyses/Pre_Filtered/Significant_Pathways_by_VEP-consequences_SeparationPathways.txt"
    fo = open(OUTFILE, 'w')
    fo.write('\t'.join(["pID", "Pathway_Name", "Aberrant", "Non-Aberrant", "p-value", "FDR", "Subject_IDs", "Genes", "Pathway_Notes"]) + '\n')
    for i in P:
        Pname = P[i][0]
        aberrant = P[i][1]
        non_aberrant = P[i][2]
        p_val = P[i][3]
        FDR = P[i][4]
        Sam_List = P[i][5]
        Gene_List = P[i][6]
        Pathway_Notes = i + Notes[i]
        fo.write('\t'.join([i, Pname, str(aberrant), str(non_aberrant), str(p_val), str(FDR), Sam_List, Gene_List, Pathway_Notes]) + '\n')
    fo.close()
    sorted = OUTFILE.rstrip(".txt") + "_sorted.txt"
    command = "(head -n 1 " + OUTFILE + " && tail -n +2 " + OUTFILE + " |sort -k 1) > " + sorted
    os.system(command)

def family_tree():
    infile = "/home/exacloud/lustre1/jjacobs/Reactome_Downloads/ReactomePathwaysRelation_06_21_20.txt"
    with open(infile, 'r') as fi:
        G = {}
        for i in fi:
            i = i.split()
            parent = i[0]
            child = i[1]
            if parent in G:
                G[parent].append(child)
            else:
                G[parent] = [child]
        return(G)

def whos_yo_daddy(G, P):
    L = []
    for j in P:
        for i in G:
            if j in G[i]:
                if i not in L:
                    L.append(i)
    return(L)

def climb_the_tree(G, Pathways):
    D = {}
    for i in Pathways:
        p = i
        s = []
        v = whos_yo_daddy(G=G, P=p)
        while len(v) > 0:
            for j in v:
                if j not in s:
                    s.append(j)
            v = whos_yo_daddy(G=G, P=v)
        n = ""
        for t in Pathways:
            if t in s:
                n += t + ","
        n = n.rstrip(",")
        if len(n) > 0:
            D[i[0]] = " is embeded in pathway(s) " + n
        else:
            D[i[0]] = " is all alone."
    return(D)

COM_DIR = "/home/exacloud/lustre1/jjacobs"

L = get_survival_data_Sample_IDs()
P = {}
P = compile_pathway_data(P=P, COM_DIR=COM_DIR, CO=CO)
G = family_tree()
Notes = climb_the_tree(G=G, Pathways=P)
write_output(P=P, COM_DIR=COM_DIR, CO=CO, Notes=Notes)
