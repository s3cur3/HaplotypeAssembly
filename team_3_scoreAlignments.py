'''
A program to find the best overlap between all sequence fragments and write them
to disk.

Tyler Young and Bob Barnhart
Written for Python 3
'''


# If a fragment in our file is less than or equal to MIN_LENGTH, we ignore it
# entirely (with 5x coverage, it may wind up being a red herring)
MIN_LENGTH = 3

# If there are greater than this percentage of errors among the "matched" base pairs,
# we reject this as a possible alignment
ALLOWED_ERROR_RATE = 0.8

import math

def overlap( s1, s2 ):
    '''
    Calculates the overlap between 2 strings. The overlap between strings 1 and 2
    is the length of the longest prefix of string 2 that matches a suffix of
    string 1.

    E.g., overlap( ABCDEFG, EFGHIJK ) = 3

    Thus, *you* are responsible for deciding whether overlap(s1,s2) is greater
    than overlap(s2,s1). Furthermore, you are responsible for handling any
    reverse-complimentarity issues.
    '''

    maxSoFar = 0
    s2pos = 0
    alignmentStart = len(s1)
    for s1pos in range(len(s1)-1, -1, -1):
        errorsSoFar = 0

        maxAtThisSize = 0
        potentialAlignmentStart = s1pos
        if s1pos < 0: # no suffix of S1 can begin left of its first character
            break
        for i in range(0, s2pos + 1):
            #allowedErrorsHere = math.ceil( (i+1)*ALLOWED_ERROR_RATE )
            allowedErrorsHere = math.ceil( math.sqrt(i)*ALLOWED_ERROR_RATE )

            if i >= len(s2): # if s2 matches completely with an internal section of s1
                break

            if s2[i] == s1[s1pos + i]:
                maxAtThisSize += 1
            else:
                errorsSoFar += 1
                maxAtThisSize += 1
                #print("Found an error! This is number ", errorsSoFar, "for this alignment, of a maximum of ", maxAllowedErrors)
                if errorsSoFar >= allowedErrorsHere:
                    # Stop considering this as a possible alignment
                    break

        s2pos += 1
        if(maxAtThisSize > maxSoFar) and (errorsSoFar < allowedErrorsHere):
            maxSoFar = maxAtThisSize
            alignmentStart = potentialAlignmentStart

    #print("Found max overlap of", string1, "and", string2, "to be", maxSoFar)
    return int(maxSoFar), alignmentStart


def getReverseCompliment(string):
    '''
    @param string Iterable The DNA sequence for which we should calculate the
                           reverse compliment
    @return Array The reverse compliment of the input string
    '''
    revComp = []
    for letter in string:
        if letter == 'A':
            revComp.append('T')
        elif letter == 'T':
            revComp.append('A')
        elif letter == 'C':
            revComp.append('G')
        elif letter == 'G':
            revComp.append('C')
        else:
            print("ERROR! Non-DNA letter",letter,"found in your FASTA data.")
            exit()
    revComp.reverse()
    revComp = "".join(revComp)
    return revComp

def getFragments(file):
    f = open(file, 'rU')
    resultList = []

    # initialize the sequence accumulator
    sequence = ''

    # process all the rest of the lines in the file
    for line in f:
        line = line.rstrip()

        # ignore blank lines
        if line == '':
            continue

        # if it's a header line, finish the previous sequence
        # and start a new one
        if line[0] == '>':
            if sequence != "" and len(sequence) > MIN_LENGTH:
                resultList.append(sequence)
            sequence = ''

        # if we're here, we must be in letters of the sequence
        else:
            sequence += line

    # Add the final sequence to the list (whatever's remaining)
    if sequence != '':
        resultList.append(sequence)
    return resultList

def getDuplicateMatrix( matrixToCopy ):
    m = []
    for row in matrixToCopy:
        m.append( list(row) )
    return m


def getOverlapMatrix(fragments):
    numSeqs = len(fragments)
    overlapMatrix = [ [0] * numSeqs for i in range(numSeqs) ]

    # Creating a matrix of the overlap distances
    # Comparing all fragments to one another
    for i in range(numSeqs):
        for j in range(numSeqs):
            overlapMatrix[i][j], theOffset = overlap(fragments[i],fragments[j])

    return overlapMatrix

def getRevCompMatrix(fragments):
    numSeqs = len(fragments)
    revCompMatrix = [ [0] * numSeqs for i in range(numSeqs) ]

    # Create a matrix of the overlap distances if each fragment, in turn,
    # is treated as part of the anti-sense strand
    for i in range(len(fragments)):
        revCompOfI = getReverseCompliment(fragments[i])
        for j in range(len(fragments)):
            revCompMatrix[i][j], theOffset = overlap(revCompOfI,fragments[j])

    return revCompMatrix

def negateMainDiagonal(squareMatrix):
    for i in range(len(squareMatrix)):
        squareMatrix[i][i] = -squareMatrix[i][i]
    return squareMatrix

def removeDuplicates(fragments):
    workingFragmentList = list(fragments)
    for i in range(len(fragments)):
        for j in range(len(fragments)):
            if fragments[i] in fragments[j] and i != j:
                workingFragmentList.remove(fragments[i])
                break
    return workingFragmentList

def main():
    workingFragmentList = getFragments('fragments.fasta')
    workingFragmentList = removeDuplicates(workingFragmentList)
    fragments = list(workingFragmentList)

    #### TODO: Should we actually be iterating this a few times? ####
    # Iterate to make sure we remove anything that aligns better on the antisense
    # strand than the sense strand
    for iteration in range(4):
        overlapMatrix = negateMainDiagonal( getOverlapMatrix( workingFragmentList ) )
        revCompMatrix = negateMainDiagonal( getRevCompMatrix( workingFragmentList ) )

        # Remove from the list of fragments any fragment which aligns better as part
        # of the anti-sense strand
        for i in range(len(fragments)):
            maxWhenFirst = max(overlapMatrix[i])
            maxWhenSecond = -1
            for row in overlapMatrix:
                maxWhenSecond = max(maxWhenSecond, row[i])

            revCompMaxWhenFirst = max(revCompMatrix[i])
            revCompMaxWhenSecond = -1
            for row in revCompMatrix:
                revCompMaxWhenSecond = max(revCompMaxWhenSecond, row[i])

            #if (maxWhenFirst + maxWhenSecond) < (revCompMaxWhenFirst + revCompMaxWhenSecond):
            if (maxWhenFirst < revCompMaxWhenFirst) \
                    and (maxWhenSecond < revCompMaxWhenSecond):
                workingFragmentList.remove(fragments[i])

        # "Delete" our knowledge of any fragments which we've decided to treat as part
        # of the antisense strand
        numDeletedFragments = len(fragments) - len(workingFragmentList)
        if iteration == 0:
            print("Decided that",numDeletedFragments,
                  "fragments belong on the anti-sense strand.")
        else:
            print("Decided another",numDeletedFragments,
                  "fragments belong on the anti-sense strand.")
        fragments = list(workingFragmentList)

        #### End possible for loop ####

    print("After removing duplicated content and fragments that fit better on the "
          + "anti-sense strand, we have", len(fragments), "fragments.")

    # Using that new, trimmed-down list of fragments, recreate the overlap matrix
    overlapMatrix = getOverlapMatrix( fragments )


    # Write the matrix to disk
    overlapFile = open("overlap.txt", "w")
    for row in overlapMatrix:
        for col in row:
            overlapFile.write(str(col))
            overlapFile.write(" ")
        overlapFile.write("\n")

    # Write the fragment file to disk
    fragmentFile = open("fragments.txt", "w")
    for f in fragments:
        fragmentFile.write(f)
        fragmentFile.write("\n")

    print("\nWrote the output file. Now run the team_3_tsp.py program in python2.")


    #print("Testing: best score of alignment for GGATGTCCTGATCCAACATCGAGGTCGTAAACCCTATTGTTGA and TCCAACATC is:")
    #print(overlap("GGATGTCCTGATCCAACATCGAGGTCGTAAACCCTATTGTTGA", "TCCAACATC"))

    print("Testing: best score of alignment for ACAAGTCCAAATTTTTGG and GGGGG is:")
    print(overlap("ACAAGTCCAAATTTTTGG", "GGGGG"))

    print("Testing: best score of alignment for TTGACCAACGCAACAAGTTAGCCTAGCGATGGCAGCGC and TGG is:")
    print(overlap("TTGACCAACGCAACAAGTTAGCCTAGCGATGGCAGCGC", "TGG"))



    # Count the base pairs:
    count = 0
    for fragment in fragments:
        count += len(fragment)
    print("\n\nNumber of base pairs is",count)
    print("Expected length of the final sequence is",count,"/ 5 = ",count/5)

if __name__ == "__main__":
    main()