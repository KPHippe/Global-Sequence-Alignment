import random
import sys
import math
import time
import csv
#scoring global variables
match = 1
mismatch = -1
space = -2
'''testing/debug and extra info specifiers'''
#debug prints majority of steps in solving matrix and backtracing
debug = False
#test initiates the test mode
#test mode generates many sequences in a range and writes
#times into csv
test = True
#printing prints the final alignment of two sequencese
#turn off for large sequences
printing = True

class Cell: 
    def __init__(self, value, parentCell, coordinate):
        self.value = value
        self.parentCell = list(parentCell)
        self.coordinate = coordinate


def main(argv):
    sequences = []
    #default sequence if the user does not do the test function 
    s1 = "GATCGGCAT"
    s2 = "CAATGTGAATC"

    #worst case would be an n x n board, so we will measure with n x n boards for worst case scenarios
    if test:
        #10 trials per (for small n) for an accurate average
        trials = 10
        for n in range(10, 2001, 10):
            averageTime = 0
            totalTime = 0
            for trial in range(trials): 
                sequences = generateWorstRandomSequences(n)
                s1 = sequences[0]
                s2 = sequences[1]
                if debug: print(*s1)
                if debug: print(*s2)
                print("n: " , n, "trial: ", trial)
                startTime = time.time()
                scoreMatrix = score(s1, s2)
                results = backtrace(scoreMatrix, s1, s2)
                #still need to find and print final tree
                endTime = time.time()
                totalTime += endTime - startTime
                if printing: printResuls(results)
                
                
            averageTime = totalTime / trials
            
            print ("n: ", n)
            print("average time: ", averageTime)
            #write it to csv here
            with open('GeneAlignmentDataNew.csv', 'a') as file : 
                writer = csv.writer(file)
                writer.writerow([n, averageTime])

    else:
        scoreMatrix = score(list(s1), list(s2))
        results = backtrace(scoreMatrix, s1, s2)
        if printing: printResuls(results)

        
''' Start of analysis '''
def score(seqX, seqY):
    n = len(seqX) + 1
    m = len(seqY) + 1
    if debug: print (n)
    if debug: print (m)
    #create board
    result = [[0 for x in range (m)] for y in range(n)]
    for x in range(n):
        for y in range(m):
            result[x][y] = Cell(0, list(), (x,y))
    #initialize board, row 0 and column 0 are initialized to mismatch accumulation
    for x in range(n):
        for y in range(m):
            if y == 0 and x >= 1:
                result[x][y].value = result[x-1][y].value + space
                result[x][y].parentCell.append(result[x-1][y])
            
            if x == 0 and y >= 1:
                result[x][y].value = result[x][y-1].value + space
                result[x][y].parentCell.append(result[x][y-1])

    #solve table for each value according to sequences x and y
    for x in range(n):
        for y in range(m):
            #verify not in first column or row
            if x > 0 and y > 0:
                previousDiagonal = result[x-1][y-1]
                previousDiagonalVal = result[x-1][y-1].value
                previousAbove = result[x][y-1]
                previousAboveVal = result[x][y-1].value
                previousLeft = result[x-1][y]
                previousLeftVal = result[x-1][y].value
                if debug: print(previousLeftVal)
                if debug: print(previousDiagonalVal)
                if debug: print(previousAboveVal)

                if seqX[x-1] == seqY[y-1]:
                    maxScore = max(previousLeftVal + space, previousAboveVal + space, previousDiagonalVal + match)
                    maxCells = list()
                    if debug: print("max score = ", maxScore)
                    if maxScore == previousDiagonalVal + match:
                        if debug: print("inside equals diagonal parent")
                        maxCells.append(previousDiagonal)
                        if debug: print("maxCells len equals" , len(maxCells))
                    if maxScore == previousLeftVal + space:
                        if debug: print("inside equals left parent")
                        maxCells.append(previousLeft)
                        if debug: print("maxCells len equals" , len(maxCells))
                    if maxScore == previousAboveVal + space:
                        if debug: print("inside equals above parent")
                        maxCells.append(previousAbove)
                        if debug: print("maxCells len equals" , len(maxCells))
                    
                    result[x][y].value = maxScore
                    result[x][y].parentCell = maxCells

                elif seqX[x-1] != seqY[y-1]:
                    maxScore = max(previousLeftVal + space, previousAboveVal + space, previousDiagonalVal + mismatch)
                    maxCells = list()
                    if debug: print("max score = ", maxScore)
                    if maxScore == previousDiagonalVal + mismatch:
                        if debug: print("inside not equals diagonal parent")
                        maxCells.append(previousDiagonal)
                        if debug: print("maxCells len not equals" , len(maxCells))
                    if maxScore == previousLeftVal + space:
                        if debug: print("inside not equals left parent")
                        maxCells.append(previousLeft)
                        if debug: print("maxCells len not equals" , len(maxCells))
                    if maxScore == previousAboveVal + space:
                        if debug: print("inside not equals above parent")
                        maxCells.append(previousAbove)
                        if debug: print("maxCells len not equals" , len(maxCells))

                    result[x][y].value = maxScore
                    result[x][y].parentCell = maxCells
                
                if debug: print("parents: " , len(maxCells))
                if debug: print("comparing " + seqX[x-1] + " and " + seqY[y-1] + " result is " + str(result[x][y].value),
                     " at index " + str(result[x][y].coordinate) , "parents values are " , [parent.value for parent in result[x][y].parentCell])
                if debug: print("\n\n")

                
    
    if debug: print("\n**************************\n")
    if debug: printMatrix(result)
    return result

'''Final backtrace methods'''
def backtrace(matrix, s1, s2):
    finalS1 = []
    finalS2 = []
    symbols = []
    score = 0
    curCell = matrix[len(matrix)-1][len(matrix[0])-1]
    while len(curCell.parentCell) > 0:
        parentX = curCell.parentCell[0].coordinate[0]
        parentY = curCell.parentCell[0].coordinate[1]
        #diagonal
        if parentX + 1 == curCell.coordinate[0] and parentY + 1 == curCell.coordinate[1]:
            if s1[parentX] == s2[parentY]: 
                score += match
                symbols.insert(0, '+')
            elif s1[parentX] != s2[parentY]: 
                score += mismatch
                symbols.insert(0, '-')
            finalS1.insert(0, s1[parentX])
            finalS2.insert(0, s2[parentY])
        #above skip
        elif parentX + 1 == curCell.coordinate[0] and parentY == curCell.coordinate[1]:
            finalS1.insert(0, s1[parentX])
            finalS2.insert(0, '-')
            symbols.insert(0, '*')
            score += space
        elif parentX == curCell.coordinate[0] and parentY + 1 == curCell.coordinate[1]:
            finalS1.insert(0, '-')
            finalS2.insert(0, s2[parentY])
            symbols.insert(0, '*')
            score += space

        
        curCell = curCell.parentCell[0]
    return [finalS1, finalS2, symbols, score]

''' Helper methods for generating inputs'''
#if no second length provided it will be a square
def generateRandomSequences(n1, n2 = -1):
    basePairs = ["A", "T", "G", "C"]
    sequence1 = []
    sequence2 = []
    if n2 == -1:
        n2 = random.randint(1,n1-1)

    for x in range(n1):
        sequence1.append(basePairs[random.randint(0,3)])
    for x in range(n2):
        sequence2.append(basePairs[random.randint(0,3)])

    return [sequence1, sequence2]

def generateWorstRandomSequences(n):
    basePairs = ["A", "T", "G", "C"]
    sequence1 = []
    sequence2 = []

    for x in range(n):
        sequence1.append(basePairs[random.randint(0,3)])
        sequence2.append(basePairs[random.randint(0,3)])
    
    return [sequence1, sequence2]


def generateReverseSequences(n1):
    basePairs = ["A", "T", "G", "C"]
    sequence1 = []
    sequence2 = []

    for x in range(n1):
        sequence1.append(basePairs[random.randint(0,3)])
        sequence2.insert(0,basePairs[random.randint(0,3)])

    return [sequence1, sequence2]


''' helper methods for debugging/printing '''
def printMatrix(arr):
    print('\n'.join([''.join(['{:3}'.format(item.value) for item in row])
    for row in arr]))

def printResuls(results):
    #results[0] is sequence 1 with correct alignment
    #results[1] is sequence 2 with correct alignment 
    #results[2] is the sequence of symbols
    #results[3] is the total score
    print("\n*************\n")
    print(*results[0])
    print(*results[1])
    print(*results[2], "({})".format(results[3]))


if __name__ == '__main__':
    main(sys.argv[1:])