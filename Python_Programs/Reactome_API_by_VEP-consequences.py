import requests
import os

all_runs = ["01101_M1", "01101_M2", "01101_M3", "01101_M4", "01101_M5", "01101_M6", \
            "01101_M7", "01101_M8", "01103_D1", "01104_M1", "01105_D1", "01106_D1", \
            "01107_M1", "01107_M2", "01108_M1", "01109_D1", "01110_D1", "01111_M1", \
            "01112_D1", "01112_M1", "01112_M2", "01115_D1", "01117_D1", "01118_D1", \
            "01119_D1", "01120_D1", "01123_D1", "01124_D1", "01125_D1", "01126_D1", \
            "01127_D1", "01128_D1", "01131_D1", "01_M", "02_D", "03_D", "04_D", "05_D", \
            "06_D", "07_D", "08_D", "09_D", "10_D", "10_M", "11_D", "11_M", "12_D", \
            "12_M", "13_D", "13_M", "14_D", "15_D", "16_D", "17_D", "18_D", "19_D"]


#This program uses the "requests" module to send http calls to the Reactome API and retrieve the
#results files from the analyses.

COM_DIR="/home/exacloud/lustre1/jjacobs"

for i in all_runs:
    INFILE = COM_DIR + "/data/osteo/SJOS0" + i + "/reactome_input_by_VEP-consequences.txt"
    if os.path.isfile(INFILE):
        statinfo = os.stat(INFILE)
        print("Size of infile: " + str(statinfo.st_size))
        if statinfo.st_size != 0:
            print("Currently analyzing " + INFILE + '\n')
            OUTFILE_1 = COM_DIR + "/data/osteo/SJOS0" + i + "/VEP_prefiltered_by_VEP-consequences_result.csv"
            OUTFILE_2 = COM_DIR + "/data/osteo/SJOS0" + i + "/VEP_prefiltered_by_VEP-consequences_not_found.csv"
            OUTFILE_3 = COM_DIR + "/data/osteo/SJOS0" + i + "/VEP_prefiltered_by_VEP-consequences_result.json.gz"
            nitbits = {"interactors": "true", "sortBy": "ENTITIES_PVALUE", "order": "ASC", "resource": "TOTAL", "pValue": 1, "includeDisease": "true"}
            with open(INFILE) as fh:
                r = requests.post("https://reactome.org/AnalysisService/identifiers/form/projection", data = nitbits, files={'file': fh})
            print(r.url)
            toke = r.json()["summary"]["token"]
            with open(OUTFILE_1, 'wb') as fo1:
                u = requests.get("https://reactome.org/AnalysisService/download/" + toke + "/pathways/TOTAL/result.csv")
                fo1.write(u.content)
            with open(OUTFILE_2, 'wb') as fo2:
                v = requests.get("https://reactome.org/AnalysisService/download/" + toke + "/entities/notFound/not_found.csv")
                fo2.write(v.content)
            with open(OUTFILE_3, 'wb') as fo3:
                w = requests.get("https://reactome.org/AnalysisService/download/" + toke + "/result.json.gz")
                fo3.write(w.content)
