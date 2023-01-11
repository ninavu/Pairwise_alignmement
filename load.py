import readline

# a function that takes the name of the FASTA file and returns all the DNA/protein sequences in a string
def loadSeq(fileName):
    fin = open(fileName, "r")
    first_line = fin.readline() # skip the first line
    lines = fin.readlines()
    
    str = ""
    for line in range(len(lines)):
        str += lines[line].strip()
    fin.close()
    return str

# a function that takes the name of a dnaMatrix file and returns scores and gap penalty in a dictionary
def loadMatrix(fileName):
    fin = open(fileName, "r")
    first_line = fin.readline() # skip the first line
    sec_line = fin.readline() # skip the 2nd line
    lines = fin.readlines()

    adict = {}    
    adict["match score"]= int(lines[1].split()[1])
    adict["mismatch score"] = int(lines[1].split()[2])
    adict["gap penalty"] = int(lines[-1].strip())

    fin.close()
    return adict

# a function that takes the name of a BLOSUM file and returns scores and gap penalty in a dictionary of dictionaries
def loadBLOSUM(fileName):
    fin = open(fileName, "r")
    first_line = fin.readline() # skip the first line
    sec_line = fin.readline() # skip the 2nd line
    lines = fin.readlines()

    pairList = {}
    proteinList = lines[0].split()
    for i in range(len(proteinList)):
        scoreDict = {}
        for j in range(len(proteinList)):
            scoreDict[proteinList[j]] = int(lines[i+1].split()[j+1])
        pairList[proteinList[i]] = scoreDict
    
    pairList["gap penalty"] = int(lines[-1].strip())

    fin.close()
    return pairList

# def main():
#     print(loadSeq("dnaseq1.txt"))
#     # print(loadSeq("dnaseq2.txt"))
#     # print(loadSeq("proteinseq1.txt"))
#     # print(loadSeq("proteinseq2.txt"))
#     # print(loadSeq("proteinseq3.txt"))
#     # print(loadSeq("proteinseq4.txt"))

#     # print(loadMatrix("dnaMatrix.txt"))
#     # print(loadMatrix("dnaMatrix2.txt"))
#     # print(loadBLOSUM("BLOSUM45.txt"))

# main()