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
PEPTIDE_INPUT_FILE = "msms 1.txt"  #LFQ
INC_OUTPUT_FILE_NAME = "test-MQ-Inclusion_90s.csv"
EXC_OUTPUT_FILE_NAME = "test-MQ-Exclusion_90s.csv"
PROT_OUTPUT_FILE_NAME = "proteins_symbols_MQ.csv"
PEAK_OUTPUT_FILE_NAME = "test-MQ-90s_Peaks.csv"
#Compound names
PEPTIDES_PER_PROTEIN = 2

#output assumptions
ADDUCT = "+H" #Almost always true
RT_WINDOW = 1.5 #Dr. Kelly says this is a good amount
GRADIENT_LENGTH = 60
MAX_MASS = 1600
MIN_MASS = 375
NCE = 45
#Import data

allPSM = pd.read_table(PEPTIDE_INPUT_FILE,sep="\t")
    
print(allPSM.shape)
# filter data
#exclusion list
badPSM = allPSM[((allPSM["Contaminant"]=="+") |
                # (allPSM["Unique (Proteins)"] == "no") |
                (allPSM["Missed cleavages"] > 0))] # multiple Charge are too hard to deal with, even if they are bad
print(badPSM.shape)


# print(allPSM["Contaminant"]=="+")

#inclusion list
allPSM = allPSM[(allPSM["Contaminant"]!="+") &
                # (allPSM["PeptideProphet Probability"] >= 0.99) &
                # (allPSM["Unique (Proteins)"] == "yes") &
                (allPSM["Missed cleavages"] == 0) &                
                ~(allPSM["Gene Names"].str.contains(";", na=False))]


print(allPSM.shape)

# print(allPSM["Retention"])
#combining data for replicate psms for each peptide
allPSM = allPSM.groupby("Sequence").agg( protein = ("Gene Names", "first"),
                                        conf = ("Score",nanmean),                                       
                                        retention = ("Retention time", nanmean),
                                        rt_stdev = ("Retention time", nanstd),
                                        rt_max = ("Retention time", "max"),
                                        rt_min = ("Retention time", "min"),
                                        charge = ("Charge", "first"),
                                        mass = ("Mass", "first")).reset_index()
# print(allPSM["retention"])
allPSM["rt_range"] = allPSM["rt_max"] - allPSM["rt_min"]
# allPSM["missingValueRate"] = 100* (total_files - allPSM["count"]) / total_files
allPSM["mz"] = allPSM["mass"] / allPSM["charge"].astype(int)
allPSM = allPSM[(allPSM["mz"]>MIN_MASS)&(allPSM["mz"]<MAX_MASS)]

# allPSM = allPSM[(allPSM["missingValueRate"] >= 0)]
allPSM = allPSM.sort_values(by="conf",ascending=False)
allPSM.to_csv(PEAK_OUTPUT_FILE_NAME,index=False)
allPSM = allPSM.dropna(subset=["protein"])

badPSM = badPSM.groupby("Sequence").agg( protein = ("Gene Names", "first"),
                                        conf = ("Score",nanmean),                                       
                                        retention = ("Retention time", nanmean),
                                        rt_stdev = ("Retention time", nanstd),
                                        rt_max = ("Retention time", "max"),
                                        rt_min = ("Retention time", "min"),
                                        charge = ("Charge", "first"),
                                        mass = ("Mass", "first")).reset_index()

badPSM["rt_range"] = badPSM["rt_max"] - badPSM["rt_min"]
badPSM["mz"] = (badPSM["mass"]) / badPSM["charge"].dropna().astype(int)
badPSM = badPSM[(badPSM["mz"]>MIN_MASS)&(badPSM["mz"]<MAX_MASS)]
badPSM = badPSM.dropna(subset=["protein"])


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
    currentSequence = eachRow["Sequence"]
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
                  "RT Time (min)": [eachRow["retention"]], #already in minutes
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
        print(eachRow["retention"])
        newRow = {"Compound": [str(currentProtein)+ "_" +str(TimesSeen[currentProtein])],
                  "Formula": [""],
                  "Adduct": [ADDUCT],
                  "m/z": [eachRow["mz"]],
                  "z": [eachRow["charge"]],
                  "RT Time (min)": [eachRow["retention"]],
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
    currentSequence = eachRow["Sequence"]
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
                  "RT Time (min)": [eachRow["retention"]],
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
