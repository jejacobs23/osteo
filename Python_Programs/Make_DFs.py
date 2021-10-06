import sys

#Program to take pathway IDs as an input and create data frames that R can
#use to produce a Kaplan-Meier curves.
#
#This version takes into account the fact that individual patients may
#have multiple samples (M and D)

infile2_dir=sys.argv[1]
infile2_name=sys.argv[2]
outfile_dir=sys.argv[3]
outfile_name=sys.argv[4]

def Make_DF():
    S = {}
    COM_DIR = "/home/exacloud/lustre1/jjacobs"
    INFILE1 = COM_DIR + "/data/osteo/Survival_Info.txt"
    INFILE2 = COM_DIR + "/data/osteo/Reactome_Analyses/" + infile2_dir + infile2_name
    print("Reading from " + INFILE1 + " and " + INFILE2)
    fi1 = open(INFILE1, 'r')
    fi2 = open(INFILE2, 'r')
    next(fi1)
    for i in fi1:
        i = i.split()
        ID = i[0][:-1]#removes the last character from the ID
        TTE = i[1]
        SS = i[2]
        if SS == "A":
            status = 0
        elif SS == "D":
            status = 1
        cohort = i[4]
        sample_type = i[5]
        if ID not in S:#Even if there are multiple samples for a single patient, the survival info is all the same.
            S[ID] = [TTE, status, cohort, sample_type]
    next(fi2)
    for j in fi2:
        P = set()
        j = j.split()
        pID = j[0]
        OUTFILE = COM_DIR + "/data/osteo/Reactome_Analyses/" + outfile_dir + "/Survival_DataFrame_" + pID + outfile_name
        fo = open(OUTFILE, 'w')
        Sinfo = j[-2]
        Sinfo = Sinfo.split(";")
        for x in Sinfo:
            x = x.split(":")
            P.add(x[0][:-1])#removes the last character from the ID
        print("Here are the samples that have mutations in the pathway, " + pID)
        print(P)
        fo.write('\t'.join(["Sample_ID", "death", "status", "group", "cohort", "sample_type"]) + '\n')
        for u in S:
            if u in P:
                group = "aberrant"
            else:
                group = "non-aberrant"
            fo.write('\t'.join([u, str(S[u][0]), str(S[u][1]), group, S[u][2], S[u][3]]) + '\n')
        fo.close()
    fi1.close()
    fi2.close()

Make_DF()
