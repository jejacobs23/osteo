import Reactome_Tools
import os

all_runs = ["01101_M1", "01101_M2", "01101_M3", "01101_M4", "01101_M5", "01101_M6", \
            "01101_M7", "01101_M8", "01103_D1", "01104_M1", "01105_D1", "01106_D1", \
            "01107_M1", "01107_M2", "01108_M1", "01109_D1", "01110_D1", "01111_M1", \
            "01112_D1", "01112_M1", "01112_M2", "01115_D1", "01117_D1", "01118_D1", \
            "01119_D1", "01120_D1", "01123_D1", "01124_D1", "01125_D1", "01126_D1", \
            "01127_D1", "01128_D1", "01131_D1", "01_M", "02_D", "03_D", "04_D", "05_D", \
            "06_D", "07_D", "08_D", "09_D", "10_D", "10_M", "11_D", "11_M", "12_D", \
            "12_M", "13_D", "13_M", "14_D", "15_D", "16_D", "17_D", "18_D", "19_D"]


#This program compiles all the "results.cvs" files from the Reactome analysis and creates
#an output file with the Reaction ID's and how many samples had aborations associated with
#that particular pathway.

COM_DIR="/home/exacloud/lustre1/jjacobs"

S = Reactome_Tools.subSamples(runs=all_runs)

D = {}
for i in all_runs:
    ALIGNMENT_RUN = "SJOS0" + i
    print("Currently analyzing " + ALIGNMENT_RUN)
    INFILE = COM_DIR + "/data/osteo/" + ALIGNMENT_RUN + "/VEP_prefiltered_by_VEP-consequences_result.csv"
    if os.path.isfile(INFILE):
        D = Reactome_Tools.read_results(Sample_Name=ALIGNMENT_RUN, infile=INFILE, D=D)
Condensed = Reactome_Tools.condense_D(D=D, S=S)
OUTFILE_1 = COM_DIR + "/data/osteo/Reactome_Analyses/Pre_Filtered/Reactome_Pathways_by_Sample_VEP_prefiltered_by_VEP-consequences.txt"
OUTFILE_2 = COM_DIR + "/data/osteo/Reactome_Analyses/Pre_Filtered/Reactome_Pathways_by_Sample_VEP_prefiltered_condensed_by_VEP-consequences.txt"
Reactome_Tools.write_txt_from_dict(D=D, outfile=OUTFILE_1)
Reactome_Tools.write_txt_from_dict(D=Condensed, outfile=OUTFILE_2)
