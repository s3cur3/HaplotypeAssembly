
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
    maxAllowedErrors = (max(len(s2), len(s1)) * 0.05)

    maxSoFar = 0
    errorsSoFar = 0
    s2pos = 0
    alignmentStart = 0
    for s1pos in range(len(s1)-1, 0, -1):
        maxAtThisSize = 0
        potentialAlignmentStart = s1pos
        if s1pos < 0: # no suffix of S1 can begin left of its first character
            break
        for i in range(0, s2pos + 1):
            if i >= len(s2): # if s2 matches completely with an internal section of s1
                break

            #print(string1," compared to ", string2)
            #print( "i is", i, " and posInS2+i is", posInS2+i)
            if s2[i] == s1[s1pos + i]:
                maxAtThisSize += 1
            else:
                errorsSoFar += 1
                if errorsSoFar >= maxAllowedErrors:
                    break

        s2pos += 1
        if(maxAtThisSize > maxSoFar):
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
            overlapMatrix[i][j], theOffset = overlap(fragments[i],fragments[j])

    revCompMatrix = [ [0] * numSeqs for i in range(numSeqs) ]

    # Creating a matrix of the overlap distances
    # Comparing all fragments to one another
    for i in range(len(fragments)-1):
        revCompOfI = getReverseCompliment(fragments[i])
        for j in range(len(fragments)-1):
            revCompMatrix[i][j], theOffset = overlap(revCompOfI,fragments[j])


    alignments = []
    currentRow = 0
    normalFragments = [] # *not* reverse-complimentary fragments!
    while len(alignments) < len(fragments):
        for row in overlapMatrix:
            row[currentRow] = -1
        for row in revCompMatrix:
            row[currentRow] = -1

        currentRowMax = max(overlapMatrix[currentRow])
        currentRowMaxIfRevComp = max(revCompMatrix[currentRow])
        nextRow = overlapMatrix[currentRow].index(currentRowMax)
        nextRowIfRevComp = revCompMatrix[currentRow].index(currentRowMaxIfRevComp)

        if currentRowMaxIfRevComp > currentRowMax and currentRow not in normalFragments:
            # If this fragment should be treated as a reverse complement of
            # another, we ignore it (since we have 5x coverage anyway when going
            # the "right" way)

            # Since "nextRow" is better treated as a rev comp, eliminate it from
            # future consideration and continue on -- next time we'll find the "real"
            # best match for currentRow
            normalFragments.append( nextRowIfRevComp )
            continue


        alignments.append([currentRow,nextRow])
        currentRow = nextRow

    offset = 0
    for pair in alignments:
        if pair[0] == pair[1] == 0:
            break
        print( (" "*offset), fragments[pair[0]], sep="" )
        theOverlap, addToOffset = overlap(fragments[pair[0]], fragments[pair[1]])
        #print("Overlap between",fragments[pair[0]], "and", fragments[pair[1]], " was", theOverlap)
        offset += addToOffset

    offset = 0
    csvFile = open("alignments.csv", "w")
    for pair in alignments:
        if pair[0] == pair[1] == 0:
            break
        csvFile.write( (" ,"*offset) )
        for letter in fragments[pair[0]]:
            csvFile.write( letter + "," )
        csvFile.write("\n")
        theOverlap, addToOffset = overlap(fragments[pair[0]], fragments[pair[1]])
        #print("Overlap between",fragments[pair[0]], "and", fragments[pair[1]], " was", theOverlap)
        offset += addToOffset


    #for alignment in alignments:
    #    print(fragments[alignment[0]] )

main()