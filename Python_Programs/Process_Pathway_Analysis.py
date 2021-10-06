#This program takes the pathways outputed by the pathway hierarchy analysis and adds in the 
#mutation information associated with each pathway (including the samples the contributed 
#mutaitons)

def ProPath(D, V, infile, outfile):
    import os
    import requests
    if os.path.isfile(infile):
        print("Currently processing " + infile)
        fi = open(infile, 'r')
        fo = open(outfile, 'w')
        fo.write('\t'.join(["Reaction_ID", "Pathway_name", "Num_MutatedGenes", "Num_MutatedEntities", "Num_MutatedInteractors", "Num_Samples", "Total_Entities", "aveMutGenes/Entities", "aveMutEntities/Entities", "Total_Interactors", "aveiMutations/Interactors", "Sample_IDs:mGenes:mEntities:mInteractors", "Genes:SampleCounts"]) + '\n')
        for j in fi:
            j = j.rstrip() 
            if j in V and len(V[j]) > 0:
                print(j)
                r = requests.get("https://reactome.org/ContentService/data/discover/" + j) #retrieve the pathway name from Reactome API
                print(r.status_code)
                Pname = r.json()["name"]
                Pname = Pname.replace(" ", "_")
                Pname = Pname.replace("'", "_prime")
                N = [] #The list of mutated gene counts in the same order as the sample names
                L = [] #The list of mutated entities counts in the same order as the sample names
                S = [] #The list of sample names
                R = [] #The list of interactor mutation counts in the same order as the sample names
                gtotal = 0.0
                etotal = 0.0
                itotal = 0.0
                w = ""
                genes = ""
                Genes = {}
                for g in D[j]: #Use D instead of V because we want to know the overall pathway burden.
                    N.append(g[2])
                    L.append(g[3])
                    S.append(g[0])
                    R.append(g[5])
                    for h in g[1]:
                        if h in Genes:
                            Genes[h] = Genes[h] + 1
                        else:
                            Genes[h] = 1
                for k in Genes:
                    genes = genes + k + ":" + str(Genes[k]) + ";"
                genes = genes.rstrip(";")
                for z in N:
                    print(z)
                    gtotal = gtotal + float(z)
                for u in L:
                    print(u)
                    etotal = etotal + float(u)
                for q in R:
                    itotal = itotal + float(q)
                for v in range(0, len(S)):
                    w = w + S[v] + ":" + str(N[v]) + ":" + str(L[v]) + ":" + str(R[v]) + ";"
                w = w.rstrip(";")
                print("Why is this an issue?")
                print(V[j])
                print("This is the total entities : " + str(D[j][0][4]))
                print("This is the total interactors : " + str(D[j][0][6]))
                if float(gtotal) != 0.0:
                    agM = (float(gtotal))/(float(len(D[j]))) #Gives the average mutated genes per sample accross the samples that contained mutations
                    gMpE = (float(agM))/(float(D[j][0][4]))
                    gMpE = "{:.2e}".format(gMpE) #formats it in scientific notation with 2 decimals
                else:
                    gMpE = 0.0
                if float(etotal) != 0.0:
                    aM = (float(etotal))/(float(len(D[j]))) #Gives the average mutated entities per sample accross the samples that contained mutations
                    MpE = (float(aM))/(float(D[j][0][4]))
                    MpE = "{:.2e}".format(MpE) #formats it in scientific notation with 2 decimals
                else:
                    MpE = 0.0
                if float(itotal) != 0.0:
                    aiM = (float(itotal))/(float(len(D[j]))) #The average interactor mutations per sample accross the samples that contianed mutations (maybe not the best measure?)
                    iMpI = (float(aiM))/(float(D[j][0][6]))
                    iMpI = "{:.2e}".format(iMpI)
                else:
                    iMpI = 0.0
                fo.write('\t'.join([j, Pname, str(gtotal), str(etotal), str(itotal), str(len(D[j])), str(D[j][0][4]), str(gMpE), str(MpE), str(D[j][0][6]), str(iMpI), w, genes]) + '\n')
            else:
                print("Warning!!  You have pathways without a corresponding entry in V. " + j)
        fi.close()
        fo.close()
        sorted = outfile.rstrip(".txt") + "_sorted.txt"
        command = "(head -n 1 " + outfile + " && tail -n +2 " + outfile + " |sort -k 3nr) > " + sorted
        os.system(command)


#This program takes the pathways outputed by the pathway LowestLevel analysis and adds in the
#mutation information associated with each pathway (including the samples the contributed
#mutaitons).  This can also be used for Pathway_Separation_Analysis

def ProPath_LL(D, infile, outfile):
    import os
    import requests
    if os.path.isfile(infile):
        print("Currently processing " + infile)
        fi = open(infile, 'r')
        #fo = open(outfile, 'w', encoding="utf-8")#This little bit of code takes care of Reactome pathway names that have special characters in them. Only works with Python 3
        fo = open(outfile, 'wb')#This is needed for Python 2.7 since the "write" function doesn't have an "encoding" options 
        #fo.write('\t'.join(["Reaction_ID", "Pathway_name", "#MutatedGenes", "#MutatedEntities", "#MutatedInteractors", "#Samples", "Total_Entities", "aveMutGenes/Entities", "aveMutEntities/Entities", "Total_Interactors", "aveiMutations/Interactors", "Sample_IDs:mGenes:mEntities:mInteractors", "Genes:SampleCounts"]) + '\n')
        fo.write('\t'.join([b"Reaction_ID", b"Pathway_name", b"Num_MutatedGenes", b"Num_MutatedEntities", b"Num_MutatedInteractors", b"Num_Samples", b"Total_Entities", b"aveMutGenes/Entities", b"aveMutEntities/Entities", b"Total_Interactors", b"aveiMutations/Interactors", b"Sample_IDs:mGenes:mEntities:mInteractors", b"Genes:SampleCounts"]) + b'\n')
        for j in fi:
            j = j.rstrip()
            if j in D and len(D[j]) > 0:
                print(j)
                r = requests.get("https://reactome.org/ContentService/data/discover/" + j) #retrieve the pathway name from Reactome API
                Pname = r.json()["name"]
                Pname = Pname.replace(" ", "_")
                Pname = Pname.replace("'", "_prime")
                Pname = Pname.encode("utf8")#Need this line for Python 2.  In Python 3, this is unneccesasry
                N = [] #The list of mutated gene counts in the same order as the sample names
                L = [] #The list of mutated entities counts in the same order as the sample names
                S = [] #The list of sample names
                R = [] #The list of interactor mutation counts in the same order as the sample names
                gtotal = 0.0
                etotal = 0.0
                itotal = 0.0
                w = ""
                genes = ""
                Genes = {}
                for g in D[j]: #Use D instead of V because we want to know the overall pathway burden.
                    N.append(g[2])
                    L.append(g[3])
                    S.append(g[0])
                    R.append(g[5])
                    for h in g[1]:
                        if h in Genes:
                            Genes[h] = Genes[h] + 1
                        else:
                            Genes[h] = 1
                for k in Genes:
                    genes = genes + k + ":" + str(Genes[k]) + ";"
                genes = genes.rstrip(";")
                for z in N:
                    print(z)
                    gtotal = gtotal + float(z)
                for u in L:
                    print(u)
                    etotal = etotal + float(u)
                for q in R:
                    itotal = itotal + float(q)
                for v in range(0, len(S)):
                    w = w + S[v] + ":" + str(N[v]) + ":" + str(L[v]) + ":" + str(R[v]) + ";"
                w = w.rstrip(";")
                print("This is the total entities : " + str(D[j][0][4]))
                print("This is the total interactors : " + str(D[j][0][6]))
                if float(gtotal) != 0.0:
                    agM = (float(gtotal))/(float(len(D[j]))) #Gives the average mutated genes per sample accross the samples that contained mutations
                    gMpE = (float(agM))/(float(D[j][0][4]))
                    gMpE = "{:.2e}".format(gMpE) #formats it in scientific notation with 2 decimals
                else:
                    gMpE = 0.0
                    genes = "*"
                if float(etotal) != 0.0:
                    aM = (float(etotal))/(float(len(D[j]))) #Gives the average mutated entities per sample accross the samples that contained mutations
                    MpE = (float(aM))/(float(D[j][0][4]))
                    MpE = "{:.2e}".format(MpE) #formats it in scientific notation with 2 decimals
                else:
                    MpE = 0.0
                if float(itotal) != 0.0:
                    aiM = (float(itotal))/(float(len(D[j]))) #The average interactor mutations per sample accross the samples that contianed mutations (maybe not the best measure?)
                    iMpI = (float(aiM))/(float(D[j][0][6]))
                    iMpI = "{:.2e}".format(iMpI)
                else:
                    iMpI = 0.0
                print(j)
                print(Pname)
                print(str(gtotal))
                print(str(etotal))
                print(str(itotal))
                print(str(len(D[j])))
                print(str(D[j][0][4]))
                print(str(gMpE))
                print(str(MpE))
                print(str(D[j][0][6]))
                print(str(iMpI))
                print(w)
                print(genes)
                #fo.write('\t'.join([j, Pname, str(gtotal), str(etotal), str(itotal), str(len(D[j])), str(D[j][0][4]), str(gMpE), str(MpE), str(D[j][0][6]), str(iMpI), w, genes]) + '\n')
                fo.write(b'\t'.join([j.encode("utf8"), Pname, str(gtotal).encode("utf8"), str(etotal).encode("utf8"), str(itotal).encode("utf8"), str(len(D[j])).encode("utf8"), str(D[j][0][4]).encode("utf8"), str(gMpE).encode("utf8"), str(MpE).encode("utf8"), str(D[j][0][6]).encode("utf8"), str(iMpI).encode("utf8"), w.encode("utf8"), genes.encode("utf8")]) + b'\n')
            else:
                print("Warning!!  You have pathways without a corresponding entry in D. " + j)
        fi.close()
        fo.close()
        sorted = outfile.rstrip(".txt") + "_sorted.txt"
        command = "(head -n 1 " + outfile + " && tail -n +2 " + outfile + " |sort -k 3nr) > " + sorted
        os.system(command)
