import pandas as pd
import os
import re
from numpy import nanmean
from numpy import nanstd


#Required Columns "PSMs Workflow ID"	"PSMs Peptide ID"	"Confidence"	"Identifying Node"	"PSM Ambiguity"	"Annotated Peptide"
# 	"Modifications"	"# Proteins"	"Master Protein Symbols"	"Protein Symbols"	"# Missed Cleavages"	"Charge"	"Rank"	"m/z [Da]"	
# "Contaminant"	"MH+ [Da]"	"Theo. MH+ [Da]"x	"Activation Type"	"NCE [%]"	"Ion Inject Time [ms]"	"RT [min]"	"First Scan"	"Spectrum File"
# 	"Percolator q-Value"	"Percolator PEP"

#Change these each time
#PSM_FILE_PATH_AND_NAME = "TMT_HCD_Study_PSMs.txt"    #TMT
#INC_OUTPUT_FILE_NAME = "TMT_Inclusion_List.csv" 
#EXC_OUTPUT_FILE_NAME = "TMT_Exclusion_List.csv" 
FragPipe_Results_PATH = "FragPipe_Output"  #LFQ
INC_OUTPUT_FILE_NAME = "test-TMT_Inclusion_List_90s.csv"
EXC_OUTPUT_FILE_NAME = "test-TMT_Exclusion_List_90s.csv"
PROT_OUTPUT_FILE_NAME = "proteins_symbols_v21.csv"
PEAK_OUTPUT_FILE_NAME = "test-TMT_Peaks.csv"
#Compound names
PEPTIDES_PER_PROTEIN = 2

#output assumptions
ADDUCT = "+H" #Almost always true
RT_WINDOW = 1.5 #Dr. Kelly says this is a good amount
GRADIENT_LENGTH = 20
MAX_MASS = 1600
MIN_MASS = 375
NCE = 40
#Import data

run_folders = [x[0] for x in os.walk(FragPipe_Results_PATH)]
# print(run_folders)
run_folders.remove(FragPipe_Results_PATH)
run_folders.remove(FragPipe_Results_PATH+"\\tmt-report")
allPSM = pd.DataFrame()
total_files = 0
for eachFolder in run_folders:
    currentPSM = pd.read_table(eachFolder+"\\psm.tsv",sep="\t")
    allPSM = pd.concat([allPSM,currentPSM])
    total_files = total_files + 1

print(allPSM.shape)
# filter data
#exclusion list
badPSM = allPSM[(allPSM["Protein"].str.contains("contam_sp")) |
                (allPSM["Is Unique"] == False) |
                (allPSM["Number of Missed Cleavages"] > 0)]

#inclusion list
allPSM = allPSM[~(allPSM["Protein"].str.contains("contam_sp")) &
                # (allPSM["PeptideProphet Probability"] >= 0.99) &
                (allPSM["Is Unique"] == True) &
                (allPSM["Number of Missed Cleavages"] == 0)]


print(allPSM.shape)
print(badPSM.shape)

# print(allPSM["Retention"])
#combining data for replicate psms for each peptide
allPSM = allPSM.groupby("Peptide").agg( protein = ("Gene", "first"),
                                        conf = ("PeptideProphet Probability",nanmean),                                       
                                        retention = ("Retention", nanmean),
                                        rt_stdev = ("Retention", nanstd),
                                        rt_max = ("Retention", "max"),
                                        rt_min = ("Retention", "min"),
                                        mz = ("Observed M/Z", "first"),
                                        mz_calc = ("Calculated M/Z", "first"),
                                        charge = ("Charge", "first"),
                                        count = ("Spectrum File", "count")).reset_index()
# print(allPSM["retention"])
allPSM["rt_range"] = allPSM["rt_max"] - allPSM["rt_min"]
# allPSM["missingValueRate"] = 100* (total_files - allPSM["count"]) / total_files
allPSM["mz_error"] = (allPSM["mz"] - allPSM["mz_calc"]) / allPSM["mz_calc"]
allPSM = allPSM[(allPSM["mz"]>MIN_MASS)&(allPSM["mz"]<MAX_MASS)]

# allPSM = allPSM[(allPSM["missingValueRate"] >= 0)]
allPSM = allPSM.sort_values(by="conf",ascending=False)
allPSM.to_csv(PEAK_OUTPUT_FILE_NAME,index=False)

badPSM = badPSM.groupby("Peptide").agg( protein = ("Gene", "first"),
                                        conf = ("PeptideProphet Probability",nanmean),                                       
                                        retention = ("Retention", nanmean),
                                        rt_stdev = ("Retention", nanstd),
                                        rt_max = ("Retention", "max"),
                                        rt_min = ("Retention", "min"),
                                        mz = ("Observed M/Z", "first"), 
                                        mz_calc = ("Calculated M/Z", "first"),
                                        charge = ("Charge", "first"),
                                        count = ("Spectrum File", "count")).reset_index()

badPSM["rt_range"] = badPSM["rt_max"] - badPSM["rt_min"]
badPSM["mz_error"] = (badPSM["mz"] - badPSM["mz_calc"]) / badPSM["mz_calc"]
badPSM = badPSM[(badPSM["mz"]>MIN_MASS)&(badPSM["mz"]<MAX_MASS)]


#Adds the proteins and peptides in order of increasing confidence
InclusionList = pd.DataFrame({"Compound": [],
                              "Formula": [],
                              "Adduct": [],
                              "m/z": [],
                              "z": [],
                              "RT Time (min)": [],
                              "Window (min)": [],
                              "HCD Collision Energies (%)": []
                              })
Proteins_Included = []
Proteins_Full = []
Peptides_Included = []
Peptides_Excluded = []
TimesSeen = {}


#populate inclusion list
for index, eachRow in allPSM.iterrows():
    currentProtein = eachRow["protein"]
    currentProtein = re.sub("sp\\|","", currentProtein)
    currentProtein = re.sub("\\|.*","",currentProtein)
    currentSequence = eachRow["Peptide"]
    if currentProtein not in Proteins_Included and currentSequence not in Peptides_Included:
        Proteins_Included.append(currentProtein)
        Peptides_Included.append(currentSequence)
        TimesSeen[currentProtein] = 1
        if TimesSeen[currentProtein] >= PEPTIDES_PER_PROTEIN:
            Proteins_Full.append(currentProtein)
        newRow = {"Compound": [str(currentProtein) + "_1"],
                  "Formula": [""],
                  "Adduct": [ADDUCT],
                  "m/z": [eachRow["mz"]],
                  "z": [int(eachRow["charge"])],
                  "RT Time (min)": [eachRow["retention"]/60],
                  "Window (min)": [RT_WINDOW],
                  "HCD Collision Energies (%)": [NCE]
                  }
        if newRow["RT Time (min)"][0] <= RT_WINDOW/2: #ELUTES TOO SOON
            newRow["RT Time (min)"] = [RT_WINDOW/2 + 0.01]
        if newRow["RT Time (min)"][0] >= GRADIENT_LENGTH - RT_WINDOW/2: #ELUTES TOO LATE
            newRow["RT Time (min)"] = [GRADIENT_LENGTH - (RT_WINDOW/2 + 0.01)]
        InclusionList = pd.concat([InclusionList, pd.DataFrame(newRow)], ignore_index = True)
    elif currentProtein not in Proteins_Full and currentSequence not in Peptides_Included:
        TimesSeen[currentProtein] = TimesSeen[currentProtein] + 1
        Peptides_Included.append(currentSequence)
        if TimesSeen[currentProtein] >= PEPTIDES_PER_PROTEIN:
            Proteins_Full.append(currentProtein)
        newRow = {"Compound": [str(currentProtein)+ "_" +str(TimesSeen[currentProtein])],
                  "Formula": [""],
                  "Adduct": [ADDUCT],
                  "m/z": [eachRow["mz"]],
                  "z": [eachRow["charge"]],
                  "RT Time (min)": [eachRow["retention"]/60],
                  "Window (min)": [RT_WINDOW],
                  "HCD Collision Energies (%)": [NCE]
                  }
        if newRow["RT Time (min)"][0] <= RT_WINDOW/2: #ELUTES TOO SOON
            newRow["RT Time (min)"] = [RT_WINDOW/2 + 0.01]
        elif newRow["RT Time (min)"][0] >= GRADIENT_LENGTH - RT_WINDOW/2: #ELUTES TOO LATE
            newRow["RT Time (min)"] = [GRADIENT_LENGTH - (RT_WINDOW/2 + 0.01)]
        else:
            # print(newRow["RT Time (min)"][0])
            pass
        InclusionList = pd.concat([InclusionList, pd.DataFrame(newRow)], ignore_index = True)
    else:
        pass
        
ExclusionList = pd.DataFrame({"Compound": [],
                              "Formula": [],
                              "Adduct": [],
                              "m/z": [],
                              "z": [],
                              "RT Time (min)": [],
                              "Window (min)": [],
                              "HCD Collision Energies (%)": []
                              })
#populate exclusion list
for index, eachRow in badPSM.iterrows():
    currentProtein = eachRow["protein"]
    currentProtein = re.sub("contam_","", str(currentProtein))
    currentProtein = re.sub("sp\\|","", currentProtein)
    currentProtein = re.sub("\\|.*","",currentProtein)
    currentSequence = eachRow["Peptide"]
    if currentSequence not in Peptides_Included:
        if currentProtein not in TimesSeen.keys():
            TimesSeen[currentProtein] = 0
        TimesSeen[currentProtein] = TimesSeen[currentProtein] + 1
        Peptides_Included.append(currentSequence)
        newRow = {"Compound": [str(currentProtein)+ "_"+ currentSequence + "_" + str(TimesSeen[currentProtein])],
                  "Formula": [""],
                  "Adduct": [ADDUCT],
                  "m/z": [eachRow["mz"]],
                  "z": [int(eachRow["charge"])],
                  "RT Time (min)": [eachRow["retention"]/60],
                  "Window (min)": [RT_WINDOW],
                  "HCD Collision Energies (%)": [NCE]
                  }
        if newRow["RT Time (min)"][0] <= RT_WINDOW/2: #ELUTES TOO SOON
            newRow["RT Time (min)"] = [RT_WINDOW/2 + 0.01]
        if newRow["RT Time (min)"][0] >= GRADIENT_LENGTH - RT_WINDOW/2: #ELUTES TOO LATE
            newRow["RT Time (min)"] = [GRADIENT_LENGTH - (RT_WINDOW/2 + 0.01)]
        ExclusionList = pd.concat([ExclusionList, pd.DataFrame(newRow)], ignore_index = True)

    

#Print number of proteins with each number of identifications
print(newRow) #to see format
print(len(Proteins_Included)) # 1st peptide
print(len(Proteins_Full)) # 2nd peptide
print(len(InclusionList.index)) # Total peptides
print(len(ExclusionList)) #peptides excluded

InclusionList["z"] = InclusionList["z"].astype(str).str.replace('\.[0-9]+', '',regex=True)
ExclusionList["z"] = ExclusionList["z"].astype(str).str.replace('\.[0-9]+', '',regex=True)

#Export to csv
pd.DataFrame({"Symbol": Proteins_Included}).to_csv(PROT_OUTPUT_FILE_NAME, index=False,encoding="ASCII")
InclusionList.to_csv(INC_OUTPUT_FILE_NAME, index=False,encoding="ASCII")
ExclusionList.to_csv(EXC_OUTPUT_FILE_NAME, index=False,encoding="ASCII")
