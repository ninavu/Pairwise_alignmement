import sys
from turtle import back
from load import *

# this function creates a 2D matrix with pre-determined scores for certain rows
def scoreMatrix(x, y, score):
    matrix = [[0] * (len(y)+1) for i in range(len(x)+1)]
    for i in range(len(x)+1):
        matrix[i][0] = i * score
        for j in range(len(y)+1):
            matrix[0][j] = j * score
    return matrix

# this function creates a 2D matrix to backtrack previous locations
def backMatrix(x, y):
    matrix = [[0] * (len(y)+1) for i in range(len(x)+1)]
    for i in range(len(x)+1):
        matrix[i][0] = "up"
        for j in range(len(y)+1):
            matrix[0][j] = "left"
    matrix[0][0] = 0
    return matrix

# this function decides which score to use for DNA scoring
def DNAscore(match, mismatch, seq1, seq2):
    if seq1 == seq2:
        return match
    else: 
        return mismatch

# this function calculates identities between 2 sequences
def identities(seq1, seq2):
    total = len(seq1)
    count = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            count += 1
    
    num = str(count) + "/" + str(total)
    percent = "(" + str(int(count/total * 100)) + "%)"
    return num + " " + percent

# this function takes 2 DNA sequences and use global alignment to decide the optimal score
def DNAglobal(scores, seq1, seq2):

    # get scores from the dnaMatrix file
    matchScore = scores.get("match score")
    mismatchScore = scores.get("mismatch score")
    gapPenalty = scores.get("gap penalty")

    # create a 2D matrix using list of list and fill in the basic gap penalties 
    matrix = scoreMatrix(seq1, seq2, gapPenalty)

    # create another 2D matrix to save the information for backtracking
    backtrack = backMatrix(seq1, seq2)

    # fill in all the scores and backtracking info
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            diagonal = matrix[i-1][j-1] + DNAscore(matchScore, mismatchScore, seq1[i-1], seq2[j-1])
            left = matrix[i][j-1] + gapPenalty
            up = matrix[i-1][j] + gapPenalty
            matrix[i][j] = max(diagonal, left, up)

            if matrix[i][j] == diagonal:
                backtrack[i][j] = "diagonal"
            elif matrix[i][j] == left:
                backtrack[i][j] = "left"
            else:
                backtrack[i][j] = "up"

    final_score = matrix[len(seq1)][len(seq2)]
    return backtrack, final_score

# this function takes 2 DNA sequences and use semi-global alignment to decide the optimal score
def DNAsemi_global(scores, seq1, seq2):

    # get scores from the dnaMatrix file
    matchScore = scores.get("match score")
    mismatchScore = scores.get("mismatch score")
    gapPenalty = scores.get("gap penalty")

    matrix = scoreMatrix(seq1, seq2, 0)
    backtrack = backMatrix(seq1, seq2)

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            diagonal = matrix[i-1][j-1] + DNAscore(matchScore, mismatchScore, seq1[i-1], seq2[j-1])

            if i == len(seq1):
                left = matrix[i][j-1]
            else:
                left = matrix[i][j-1] + gapPenalty

            if j == len(seq2):
                up = matrix[i-1][j]
            else:
                up = matrix[i-1][j] + gapPenalty
            
            matrix[i][j] = max(diagonal, left, up)

            if matrix[i][j] == diagonal:
                backtrack[i][j] = "diagonal"
            elif matrix[i][j] == left:
                backtrack[i][j] = "left"
            else:
                backtrack[i][j] = "up"
  
    final_score = matrix[len(seq1)][len(seq2)]
    return backtrack, final_score

# this function takes 2 DNA sequences and uses local alignment to decide the optimal score
def DNAlocal(scores, seq1, seq2):
    matchScore = scores.get("match score")
    mismatchScore = scores.get("mismatch score")
    gapPenalty = scores.get("gap penalty")

    matrix = scoreMatrix(seq1, seq2, 0)
    backtrack = backMatrix(seq1, seq2)
    max_score = 0

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            diagonal = matrix[i-1][j-1] + DNAscore(matchScore, mismatchScore, seq1[i-1], seq2[j-1])
            left = matrix[i][j-1] + gapPenalty
            up = matrix[i-1][j] + gapPenalty
            matrix[i][j] = max(diagonal, left, up)

            if matrix[i][j] == diagonal:
                backtrack[i][j] = "diagonal"
            elif matrix[i][j] == left:
                backtrack[i][j] = "left"
            else:
                backtrack[i][j] = "up"

            if matrix[i][j] < 0:
                matrix[i][j] = 0

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_index = [i,j]

    final_score = matrix[len(seq1)][len(seq2)]
    return backtrack, final_score, max_score, max_index, matrix

# this function uses backtracking info to create new aligned sequences
def aligned(seq1, seq2, src, type):
    l1 = len(seq1)
    l2 = len(seq2)
    back_matrix = src[0]
    s1 = ""
    s2 = ""
    
    if type == "G" or type == "S":
        score = src[1]
        while (l1 > 0 or l2 > 0):
            if back_matrix[l1][l2] == "diagonal":
                s1 = seq1[l1-1] + s1
                s2 = seq2[l2-1] + s2
                l1 -= 1
                l2 -= 1
            elif back_matrix[l1][l2] == "left":
                s1 = "-" + s1
                s2 = seq2[l2-1] + s2
                l2 -= 1
            else:
                s1 = seq1[l1-1] + s1
                s2 = "-" + s2
                l1 -= 1
    else:
        score = src[2]
        max_start = src[3]
        l1 = max_start[0]
        l2 = max_start[1]
        score_matrix = src[4]
        while (score_matrix[l1][l2] != 0):
            if back_matrix[l1][l2] == "diagonal":
                s1 = seq1[l1-1] + s1
                s2 = seq2[l2-1] + s2
                l1 -= 1
                l2 -= 1
            elif back_matrix[l1][l2] == "left":
                s1 = "-" + s1
                s2 = seq2[l2-1] + s2
                l2 -= 1
            else:
                s1 = seq1[l1-1] + s1
                s2 = "-" + s2
                l1 -= 1
            
    return s1, s2, score, l1, l2

# this function determines which score to use from the BLOSUM file
def proteinScore(blosum, seq1, seq2):
    score = blosum.get(seq1).get(seq2)
    return score

# this function creates a scoring matrix using global alignment for 2 protein sequences
def proteinGlobal(scores, seq1, seq2):
    gapPenalty = scores.get("gap penalty")
    matrix = scoreMatrix(seq1, seq2, gapPenalty)
    backtrack = backMatrix(seq1, seq2)
    
    # fill in all the scores and backtracking info
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            diagonal = matrix[i-1][j-1] + proteinScore(scores, seq1[i-1], seq2[j-1])
            left = matrix[i][j-1] + gapPenalty
            up = matrix[i-1][j] + gapPenalty
            matrix[i][j] = max(diagonal, left, up)

            if matrix[i][j] == diagonal:
                backtrack[i][j] = "diagonal"
            elif matrix[i][j] == left:
                backtrack[i][j] = "left"
            else:
                backtrack[i][j] = "up"

    final_score = matrix[len(seq1)][len(seq2)]
    return backtrack, final_score

# this function creates a scoring matrix using semi-global alignment for 2 protein sequences
def proteinSemi_global(scores, seq1, seq2):
    matrix = scoreMatrix(seq1, seq2, 0)
    backtrack = backMatrix(seq1, seq2)
    gapPenalty = scores.get("gap penalty")

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            diagonal = matrix[i-1][j-1] + proteinScore(scores, seq1[i-1], seq2[j-1])

            if i == len(seq1):
                left = matrix[i][j-1]
            else:
                left = matrix[i][j-1] + gapPenalty

            if j == len(seq2):
                up = matrix[i-1][j]
            else:
                up = matrix[i-1][j] + gapPenalty
            
            matrix[i][j] = max(diagonal, left, up)

            if matrix[i][j] == diagonal:
                backtrack[i][j] = "diagonal"
            elif matrix[i][j] == left:
                backtrack[i][j] = "left"
            else:
                backtrack[i][j] = "up"
  
    final_score = matrix[len(seq1)][len(seq2)]
    return backtrack, final_score

# this function takes 2 DNA sequences and uses local alignment to decide the optimal score
def proteinLocal(scores, seq1, seq2):
    matrix = scoreMatrix(seq1, seq2, 0)
    backtrack = backMatrix(seq1, seq2)
    gapPenalty = scores.get("gap penalty")
    max_score = 0

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            diagonal = matrix[i-1][j-1] + proteinScore(scores, seq1[i-1], seq2[j-1])
            left = matrix[i][j-1] + gapPenalty
            up = matrix[i-1][j] + gapPenalty
            matrix[i][j] = max(diagonal, left, up)

            if matrix[i][j] == diagonal:
                backtrack[i][j] = "diagonal"
            elif matrix[i][j] == left:
                backtrack[i][j] = "left"
            else:
                backtrack[i][j] = "up"

            if matrix[i][j] < 0:
                matrix[i][j] = 0

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_index = [i,j]

    final_score = matrix[len(seq1)][len(seq2)]
    return backtrack, final_score, max_score, max_index, matrix

# this function deals with all command-line arguments and put them in an ordered list
def param(seq1, seq2, output, proseq, atype):
    argv = sys.argv[1:]
    pars = ["-i", "-j", "-o", "-p", "-atype"]
    alist = [seq1, seq2, output, proseq, atype]

    for i in range(len(argv)):
        for j in range(len(pars)):
            if argv[i] == pars[j]:
                try:
                    argv[i+1]
                except IndexError:
                    print("Missing argument after the last index " +  pars[j] + ". Can't run the program!")
                else:
                    if argv[i+1] not in pars:
                        alist[j] = argv[i+1]

    for i in range(len(alist)):
        if alist[i] == "":
            print("Missing " + pars[i] + " or missing argument after " +  pars[i] + ". Can't run the program!")
            return None
        if alist[3] != 'T' and alist[3] != 'F':
            print("Wrong argument, can't run the program! It should be either T or F after '-p'.")
            return None
            break
        if alist[4] != 'G' and alist[4] != 'S' and alist[4] != 'L':
            print("Wrong argument, can't run the program! It should be either G or S after '-atype'.")
            return None
            break

    return alist

def main():
    alist = param(None, None, None, None, None)
    # print(alist)
    i = loadSeq(alist[0])
    j = loadSeq(alist[1])
    output = alist[2]
    proseq = alist[3]
    atype = alist[4]

    if alist != None:
        fout = open(output, "w")
        fout.write("\n \n")
        if proseq == "F":
            scoring = loadMatrix("dnaMatrix.txt")
            if atype == "G":
                source = DNAglobal(scoring, i, j)
            elif atype == "S":
                source = DNAsemi_global(scoring, i, j)
            else:
                source = DNAlocal(scoring, i, j)
        else:
            scoring = loadBLOSUM("BLOSUM45.txt")
            if atype == "G":
                source = proteinGlobal(scoring, i, j)
            elif atype == "S":
                source = proteinSemi_global(scoring, i, j)
            else: 
                source = proteinLocal(scoring, i, j)

        result = aligned(i, j, source, atype)
        align1 = result[0]
        align2 = result[1]
        score = str(result[2])
        if atype == "G" or atype == "S":
            fout.write("seq1: " + str(1) + " " + align1 + " " + str(len(i)) + "\n")
            fout.write("\n")
            fout.write("seq2: " + str(1) + " " + align2 + " " + str(len(j)) + "\n")   
        else:
            last_idx = source[3]
            fout.write("seq1: " + str(result[3]+1) + " " + align1 + " " + str(last_idx[0]+1) + "\n")
            fout.write("\n")
            fout.write("seq2: " + str(result[4]+1) + " " + align2 + " " + str(last_idx[1]+1) + "\n") 
        fout.write("\n")
        fout.write("Score: " + score + "\n")
        fout.write("Identities: " + identities(align1, align2) + "\n")
        fout.close()
main()