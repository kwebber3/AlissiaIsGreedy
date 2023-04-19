import pandas as pd
import re


#Change these each time
FILE_PATH_AND_NAME = "TMT_HCD_Study_PSMs.txt"
OUTPUT_FILE_NAME = "TMT_Inclusion_List.csv"

#Compound names
FIRST_EXTENSION = "_01"
FULL_EXTENSION = "_02"

#output assumptions
ADDUCT = "+H" #Almost always true
RT_WINDOW = 5 #Dr. Kelly says this is a good amount

#Import data

AllData = pd.read_table(FILE_PATH_AND_NAME, delimiter="\t")

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

for index, eachRow in AllData.iterrows():
    currentProtein = eachRow["Master Protein Accessions"]
    currentSequence = re.sub("\\]\\.","",eachRow["Annotated Sequence"])
    currentSequence = re.sub("\\.\\[[a-z\\-\\+]\\]","",currentSequence)
    currentSequence = re.sub("[a-z]","",currentSequence)
    currentSequence = re.sub("\\[","",currentSequence)
    if currentProtein not in Proteins_Included and currentSequence not in Peptides_Included:
        Proteins_Included.append(currentProtein)
        Peptides_Included.append(currentSequence)
        newRow = {"Compound": [str(currentProtein) + FIRST_EXTENSION],
                  "Formula": [currentSequence],
                  "Adduct": [ADDUCT],
                  "m/z": [eachRow["Theo. MH+ [Da]"]],
                  "z": [eachRow["Charge"]],
                  "RT Time (min)": [eachRow["RT [min]"]],
                  "Window (min)": [RT_WINDOW],
                  "HCD Collision Energies (%)": [eachRow["NCE [%]"]]
                  }
        InclusionList = pd.concat([InclusionList, pd.DataFrame(newRow)], ignore_index = True)
    elif currentProtein not in Proteins_Full and currentSequence not in Peptides_Included:
        Peptides_Included.append(currentSequence)
        Proteins_Full.append(currentProtein)
        newRow = {"Compound": [str(currentProtein) + FULL_EXTENSION],
                  "Formula": [currentSequence],
                  "Adduct": [ADDUCT],
                  "m/z": [eachRow["Theo. MH+ [Da]"]],
                  "z": [eachRow["Charge"]],
                  "RT Time (min)": [eachRow["RT [min]"]],
                  "Window (min)": [RT_WINDOW],
                  "HCD Collision Energies (%)": [eachRow["NCE [%]"]]
                  }
        InclusionList = pd.concat([InclusionList, pd.DataFrame(newRow)], ignore_index = True)

#Print number of proteins with each number of identifications
print(newRow) #to see format
print(len(Proteins_Included)) # 1st peptide
print(len(Proteins_Full)) # 2nd peptide
print(len(InclusionList.index)) # Total peptides

#Export to csv
InclusionList.to_csv(OUTPUT_FILE_NAME, index=False)