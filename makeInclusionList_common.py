import pandas as pd
import re

#Required Columns "PSMs Workflow ID"	"PSMs Peptide ID"	"Confidence"	"Identifying Node"	"PSM Ambiguity"	"Annotated Sequence"
# 	"Modifications"	"# Proteins"	"Master Protein Accessions"	"Protein Accessions"	"# Missed Cleavages"	"Charge"	"Rank"	"m/z [Da]"	
# "Contaminant"	"MH+ [Da]"	"Theo. MH+ [Da]"x	"Activation Type"	"NCE [%]"	"Ion Inject Time [ms]"	"RT [min]"	"First Scan"	"Spectrum File"
# 	"Percolator q-Value"	"Percolator PEP"

#Change these each time
#PSM_FILE_PATH_AND_NAME = "TMT_HCD_Study_PSMs.txt"    #TMT
#INC_OUTPUT_FILE_NAME = "TMT_Inclusion_List.csv" 
#EXC_OUTPUT_FILE_NAME = "TMT_Exclusion_List.csv" 
PSM_FILE_PATH_AND_NAME = "LFQ_HCD_Scan-(1)_PSMs.txt"  #LFQ
PEPTIDE_FILE_PATH_AND_NAME = "LFQ_HCD_Scan-(1)_PeptideGroups.txt"
INC_OUTPUT_FILE_NAME = "LFQ_Common_Inclusion_List.csv"
EXC_OUTPUT_FILE_NAME = "LFQ_Common_Exclusion_List.csv" 

#Compound names
MIN_PEPTIDES = 2
MAX_PEPTIDES = 20
PERCENT_MISSING_VALUES = 30

#output assumptions
ADDUCT = "+H" #Almost always true
RT_WINDOW = 5 #Dr. Kelly says this is a good amount
GRADIENT_LENGTH = 65



#Import data

# make sure we aren't choosing uncommon peptides
Good_Peptides = pd.read_table(PEPTIDE_FILE_PATH_AND_NAME, delimiter="\t").filter(regex="Found in File:|Annotated Sequence")
Good_Peptides["Missing Value Rate"] = (len(Good_Peptides.filter(regex="Found in File:").columns) - Good_Peptides[(Good_Peptides.filter(regex="Found in File:") == "Peak Found") | (Good_Peptides.filter(regex="Found in File:") == "High")].count(axis=1)) / len(Good_Peptides.filter(regex="Found in File:").columns) *100
Good_Peptides = Good_Peptides[Good_Peptides["Missing Value Rate"] < PERCENT_MISSING_VALUES]["Annotated Sequence"]
print(len(Good_Peptides))

#PSMs and remove uncommon psms
AllData = pd.read_table(PSM_FILE_PATH_AND_NAME, delimiter="\t")
print(AllData.shape[0])
AllData = AllData[AllData["Annotated Sequence"].isin(Good_Peptides)]
print(AllData.shape[0])


#Filter
AllData = AllData[(AllData["Contaminant"] == False) &
                  (AllData["# Proteins"] == 1) &
                  (AllData["# Missed Cleavages"] == 0) &
                  (AllData["Confidence"] == "High")]


#Finds an average QValue and sorts from lowest to highest
AllData["Average q-value"] = AllData.groupby("NCE [%]")["Percolator q-Value"].transform("mean")

AllData.sort_values(by="Average q-value")


#Adds the proteins and peptides in order of increasing confidence
InclusionList = pd.DataFrame({"Compound": [],
                              "Formula": [],
                              "Adduct": [],
                              "m/z": [],
                              "z": [],
                              "RT Time (min)": [],
                              "Window (min)": []
                              })
Proteins_Included = []
Proteins_Full = []
Peptides_Included = []
Peptides_Excluded = []
TimesSeen = {}

ExclusionList = pd.DataFrame({"Compound": [],
                              "Formula": [],
                              "Adduct": [],
                              "m/z": [],
                              "z": [],
                              "RT Time (min)": [],
                              "Window (min)": []
                              })

for index, eachRow in AllData.iterrows():
    currentProtein = eachRow["Master Protein Accessions"]
    currentSequence = re.sub("\\]\\.","",eachRow["Annotated Sequence"])
    currentSequence = re.sub("\\.\\[[A-z\\-\\+]\\]","",currentSequence)
    currentSequence = re.sub("[a-z]","",currentSequence)
    currentSequence = re.sub("\\[","",currentSequence)
    currentSequence = re.sub("\\.","",currentSequence)
    currentSequence = re.sub("\\]","",currentSequence)
    currentSequence = re.sub("\\+","",currentSequence)
    currentSequence = re.sub("\\-","",currentSequence)
    if currentProtein not in Proteins_Included and currentSequence not in Peptides_Included:
        Proteins_Included.append(currentProtein)
        Peptides_Included.append(currentSequence)
        TimesSeen[currentProtein] = 1
        if TimesSeen[currentProtein] >= MIN_PEPTIDES:
            Proteins_Full.append(currentProtein)
        newRow = {"Compound": [re.sub(",","", str(currentProtein)) + "_1"],
                  "Formula": [currentSequence],
                  "Adduct": [ADDUCT],
                  "m/z": [eachRow["m/z [Da]"]],
                  "z": [int(eachRow["Charge"])],
                  "RT Time (min)": [eachRow["RT [min]"]],
                  "Window (min)": [RT_WINDOW],
                  "HCD Collision Energies (%)": [eachRow["NCE [%]"]]
                  }
        if newRow["RT Time (min)"][0] <= RT_WINDOW/2: #ELUTES TOO SOON
            newRow["RT Time (min)"] = [RT_WINDOW/2 + 0.01]
        if newRow["RT Time (min)"][0] >= GRADIENT_LENGTH - RT_WINDOW/2: #ELUTES TOO LATE
            newRow["RT Time (min)"] = [GRADIENT_LENGTH - (RT_WINDOW/2 + 0.01)]
        InclusionList = pd.concat([InclusionList, pd.DataFrame(newRow)], ignore_index = True)
    elif currentProtein not in Proteins_Full and currentSequence not in Peptides_Included:
        TimesSeen[currentProtein] = TimesSeen[currentProtein] + 1
        Peptides_Included.append(currentSequence)
        if TimesSeen[currentProtein] >= MIN_PEPTIDES:
            Proteins_Full.append(currentProtein)
        newRow = {"Compound": [str(currentProtein) + str(TimesSeen[currentProtein])],
                  "Formula": [currentSequence],
                  "Adduct": [ADDUCT],
                  "m/z": [eachRow["m/z [Da]"]],
                  "z": [int(eachRow["Charge"])],
                  "RT Time (min)": [eachRow["RT [min]"]],
                  "Window (min)": [RT_WINDOW],
                  "HCD Collision Energies (%)": [eachRow["NCE [%]"]]
                  }
        if newRow["RT Time (min)"][0] <= RT_WINDOW/2: #ELUTES TOO SOON
            newRow["RT Time (min)"] = [RT_WINDOW/2 + 0.01]
        if newRow["RT Time (min)"][0] >= GRADIENT_LENGTH - RT_WINDOW/2: #ELUTES TOO LATE
            newRow["RT Time (min)"] = [GRADIENT_LENGTH - (RT_WINDOW/2 + 0.01)]
        InclusionList = pd.concat([InclusionList, pd.DataFrame(newRow)], ignore_index = True)
    else:
        if currentProtein in TimesSeen.keys() and currentSequence not in Peptides_Excluded and currentSequence not in Peptides_Included:
            TimesSeen[currentProtein] = TimesSeen[currentProtein] + 1
        elif currentProtein in TimesSeen.keys(): #peptide already counted
            pass
        else:
            TimesSeen[currentProtein] = 1
        if TimesSeen[currentProtein] > MAX_PEPTIDES and currentSequence not in Peptides_Excluded and currentSequence not in Peptides_Included:
            Peptides_Excluded.append(currentSequence)
            newRow = {"Compound": [str(currentProtein) + "_" + str(TimesSeen[currentProtein])],
                  "Formula": [currentSequence],
                  "Adduct": [ADDUCT],
                  "m/z": [eachRow["m/z [Da]"]],
                  "z": [int(eachRow["Charge"])],
                  "RT Time (min)": [eachRow["RT [min]"]],
                  "Window (min)": [RT_WINDOW],
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

#Export to csv
InclusionList.to_csv(INC_OUTPUT_FILE_NAME, index=False)
ExclusionList.to_csv(EXC_OUTPUT_FILE_NAME, index=False)