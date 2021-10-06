import sys

#Program to take the output of variant filtering and create a list of genes for input
#into Reactome.
#

CONS = ["3_prime_UTR_variant", "splice_acceptor_variant", "TF_binding_site_variant", \
        "protein_altering_variant", "splice_donor_variant", "stop_lost", "5_prime_UTR_variant", \
        "splice_region_variant", "stop_gained", "coding_sequence_variant", "start_lost", \
        "frameshift_variant", "regulatory_region_variant"]

all_runs = ["01101_M1", "01101_M2", "01101_M3", "01101_M4", "01101_M5", "01101_M6", \
            "01101_M7", "01101_M8", "01103_D1", "01104_M1", "01105_D1", "01106_D1", \
            "01107_M1", "01107_M2", "01108_M1", "01109_D1", "01110_D1", "01111_M1", \
            "01112_D1", "01112_M1", "01112_M2", "01115_D1", "01117_D1", "01118_D1", \
            "01119_D1", "01120_D1", "01123_D1", "01124_D1", "01125_D1", "01126_D1", \
            "01127_D1", "01128_D1", "01131_D1", "01_M", "02_D", "03_D", "04_D", "05_D", \
            "06_D", "07_D", "08_D", "09_D", "10_D", "10_M", "11_D", "11_M", "12_D", \
            "12_M", "13_D", "13_M", "14_D", "15_D", "16_D", "17_D", "18_D", "19_D"]


def write_output(com_dir, D, sample):
    output_file = com_dir + "/data/osteo/SJOS0" + sample + "/reactome_input_by_VEP-consequences.txt"
    fo = open(output_file, 'w')
    for h in D:
        fo.write(D[h] + '\n')
    fo.close()

def filter_by_CAT(TOOL, value1, value2, gate):
    print(TOOL)
    if gate == "stop":
        return(gate)
    else:
        TOOL = TOOL.split("(")
        if TOOL[0] == value1:
            gate = "go"
        elif TOOL[0] == value2:
            gate = "go"
        else:
            gate = "stop"
        return(gate)

for x in all_runs:
    D = {}
    sample = x
    common_dir="/home/exacloud/lustre1/jjacobs"
    input_file=common_dir + "/data/osteo/SJOS0" + sample + "/VEP_prefiltered.txt"
    fi = open(input_file, 'r')
    print("reading from " + input_file)

    for l in fi:
        h = []
        if l[0] != "#":
            i = l.split()
            if len(i) != 14:
                print("Abnormal line :" + l)
            else:
                ivar = i[0]
                Feature_type = i[5]
                Consequence = i[6]
                Consequence = Consequence.split(",")
                Extra = i[13]
                Extra = Extra.split(";")
                for u in Extra:
                    u = u.split("=")
                    h.append(u[0])
                    if u[0] == "SYMBOL":
                        SYMBOL = u[1]
                    elif u[0] == "NEAREST":
                        NSYMBOL = u[1]
                        NSYMBOL = NSYMBOL.split(",")
                    elif u[0] == "SIFT":
                        SIFT = u[1]
                    elif u[0] == "PolyPhen":
                        POLYPHEN = u[1]
                    elif u[0] == "CANONICAL":
                        CANONICAL = u[1]
                print(h)
                for c in Consequence:
                    if c in CONS:
                        if "NEAREST" in h:
                            if len(NSYMBOL) > 1:
                                print("This " + c + " variant is nearby more than one gene")
                                print(ivar)
                                print(NSYMBOL)
                            for n in NSYMBOL:
                                uvar = ivar + "_" + n
                                if uvar not in D:
                                    D[uvar] = n
                                elif uvar in D and "CANONICAL" in h and CANONICAL == "YES":
                                    D[uvar] = n
                                else:
                                    print("Non-unique " + c + " variant: " + uvar)
                        else:
                            print("Memphis, we have a problem: " + l)

    write_output(com_dir=common_dir, D=D, sample=sample)

    fi.close()                
