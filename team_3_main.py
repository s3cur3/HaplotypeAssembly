
def overlap( string1, string2 ):
    '''
    Calculates the overlap between 2 strings. The overlap between strings 1 and 2
    is the length of the longest prefix of string 2 that matches a suffix of
    string 1.

    Thus, *you* are responsible for deciding whether overlap(s1,s2) is greater
    than overlap(s2,s1). Furthermore, you are responsible for handling any
    reverse-complimentarity issues.
    '''
    maxAllowedErrors = (max(len(string1), len(string2)) * 0.1)

    maxSoFar = 0
    errorsSoFar = 0
    posInS1 = 0
    maxAtThisSize = 0
    for posInS2 in range(len(string2)-1, 0, -1):
        maxAtThisSize = 0
        if posInS1 > len(string1) or posInS2 < 0:
            break
        for i in range(0, posInS1 + 1):
            if i >= len(string1):
                break

            #print(string1," compared to ", string2)
            #print( "i is", i, " and posInS2+i is", posInS2+i)
            if string1[i] == string2[posInS2 + i]:
                maxAtThisSize += 1
            else:
                errorsSoFar += 1
                if errorsSoFar >= maxAllowedErrors:
                    break

        posInS1 += 1
        if(maxAtThisSize > maxSoFar):
            maxSoFar = maxAtThisSize

    #print("Found max overlap of", string1, "and", string2, "to be", maxSoFar)
    return maxSoFar


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
        else:
            revComp.append('C')
    revComp.reverse()
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
            if sequence != "":
                resultList.append(sequence)
            sequence = ''

        # if we're here, we must be in letters of the sequence
        else:
            sequence += line
    return resultList



def main():
    f = open('fragments.txt', 'rU')
    fragments = []
    # Read in each line in the file
    for line in f:
        # Only get the sequences, not fragment # or newline char
        if(line[0] != '>' and line != '\n'):
            fragments.append(line.rstrip())

    fragments = getFragments('fragments.txt')

    numSeqs = len(fragments)
    overlapMatrix = [ [0] * numSeqs for i in range(numSeqs) ]

    # Creating a matrix of the overlap distances
    # Comparing all fragments to one another
    for i in range(len(fragments)-1):
        for j in range(len(fragments)-1):
            overlapMatrix[i][j] = overlap(fragments[i],fragments[j])

    revCompMatrix = [ [0] * numSeqs for i in range(numSeqs) ]

    # Creating a matrix of the overlap distances
    # Comparing all fragments to one another
    for i in range(len(fragments)-1):
        revCompOfI = getReverseCompliment(fragments[i])
        for j in range(len(fragments)-1):
            revCompMatrix[i][j] = overlap(revCompOfI,fragments[j])


    alignments = []
    currentRow = 0
    while len(alignments) < len(fragments):
        for row in overlapMatrix:
            row[currentRow] = -1
        for row in revCompMatrix:
            row[currentRow] = -1

        currentRowMax = max(overlapMatrix[currentRow])
        currentRowMaxIfRevComp = max(revCompMatrix[currentRow])
        nextRow = overlapMatrix[currentRow].index(currentRowMax)

        if currentRowMaxIfRevComp > currentRowMax:
            # If this fragment should be treated as a reverse complement of
            # another, we ignore it (since we have 5x coverage anyway when going
            # the "right" way)

            # Since "nextRow" is better treated as a rev comp, eliminate it from
            # future consideration and continue on -- next time we'll find the "real"
            # best match for currentRow
            for row in overlapMatrix:
                row[nextRow] = -1
            for row in revCompMatrix:
                row[nextRow] = -1
            continue


        alignments.append([currentRow,nextRow])
        currentRow = nextRow



    print(alignments)

    for alignment in alignments:
        print(fragments[alignment[0]] )

main()