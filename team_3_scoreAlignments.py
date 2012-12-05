'''
A program to find the best overlap between all sequence fragments and write them
to disk.

Tyler Young and Bob Barnhart
Written for Python 3
'''


# If a fragment in our file is less than or equal to MIN_LENGTH, we ignore it
# entirely (with 5x coverage, it may wind up being a red herring)
MIN_LENGTH = 5


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
    maxAllowedErrors = (min(len(s2), len(s1)) * 0.2) # 10% of the shortest seq.

    maxSoFar = 0
    s2pos = 0
    alignmentStart = len(s1)
    for s1pos in range(len(s1)-1, 0, -1):
        errorsSoFar = 0

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
                #print("Found an error! This is number ", errorsSoFar, "for this alignment, of a maximum of ", maxAllowedErrors)
                if errorsSoFar >= maxAllowedErrors:
                    # Stop considering this as a possible alignment
                    break

        s2pos += 1
        if(maxAtThisSize > maxSoFar) and (errorsSoFar < maxAllowedErrors):
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
            if sequence != "" and len(sequence) > MIN_LENGTH:
                resultList.append(sequence)
            sequence = ''

        # if we're here, we must be in letters of the sequence
        else:
            sequence += line
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
    for i in range(numSeqs-1):
        for j in range(numSeqs-1):
            overlapMatrix[i][j], theOffset = overlap(fragments[i],fragments[j])

    return overlapMatrix

def getRevCompMatrix(fragments):
    numSeqs = len(fragments)
    revCompMatrix = [ [0] * numSeqs for i in range(numSeqs) ]

    # Create a matrix of the overlap distances if each fragment, in turn,
    # is treated as part of the anti-sense strand
    for i in range(len(fragments)-1):
        revCompOfI = getReverseCompliment(fragments[i])
        for j in range(len(fragments)-1):
            revCompMatrix[i][j], theOffset = overlap(revCompOfI,fragments[j])

    return revCompMatrix


def main():
    workingFragmentList = getFragments('fragments.fasta')
    fragments = list(workingFragmentList)

    #### TODO: Should we actually be iterating this a few times? ####
    # Iterate to make sure we remove anything that aligns better on the antisense
    # strand than the sense strand
    for iteration in range(5):
        overlapMatrix = getOverlapMatrix( workingFragmentList )
        revCompMatrix = getRevCompMatrix( workingFragmentList )

        # Remove from the list of fragments any fragment which aligns better as part
        # of the anti-sense strand
        for i in range(len(fragments)):
            if max(overlapMatrix[i]) > max(revCompMatrix[i]):
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


    ##### NOTE: All code that follows is a remnant of when we were trying to #####
    ##### solve TSP by hand. Its output should be ignored, but I've left the #####
    ##### code for the sake of posterity (i.e., when we need to do similar   #####
    ##### things in the future.                                              #####


    alignments = []
    currentFragment = 0
    senseFragments = [] # *not* reverse-complimentary fragments!
    antisenseFragments = [] # reverse-complimentary fragments
    prevOverlapMatrix = getDuplicateMatrix(overlapMatrix)
    while len(alignments) < len(workingFragmentList):
        for row in overlapMatrix:
            row[currentFragment] = -1
        for row in revCompMatrix:
            row[currentFragment] = -1

        currentFragmentOverlap = max(overlapMatrix[currentFragment])
        revCompOverlap = max(revCompMatrix[currentFragment])
        nextFragment = overlapMatrix[currentFragment].index(currentFragmentOverlap)
        nextFragmentIfRevComp = revCompMatrix[currentFragment].index(revCompOverlap)

        if revCompOverlap > currentFragmentOverlap \
                and currentFragment not in senseFragments\
                and currentFragment != 0:
            # If this fragment should be treated as a reverse complement of
            # another, we ignore it (since we have 5x coverage anyway when going
            # the "right" way)

            # Since the current fragment is better treated as a rev comp, eliminate
            # it from future consideration and continue on -- next time we'll find
            # the "real" best match for currentRow

            # Force the sense version of the current fragment to go in the "sense"
            # direction
            senseFragments.append( nextFragmentIfRevComp )

            # For posterity, remember that we forced this to be an antisense fragment
            antisenseFragments.append(currentFragment)
            # Delete the current fragment from the list of fragments (due to the
            # loop's exit condition
            workingFragmentList.remove(fragments[currentFragment])

            # Undo matching the previous strand with this one
            prevMatch = alignments.pop()
            currentFragment = prevMatch[0] # reset current fragment

            # The overlap matrix after undoing the previous alignment
            #overlapMatrix = prevOverlapMatrix # reset the overlap matrix
        else:
            alignments.append([currentFragment,nextFragment])
            currentFragment = nextFragment
            prevOverlapMatrix = getDuplicateMatrix(overlapMatrix)

    print("\nWrote the output file. Now run the team_3_tsp.py program in python2.")


    #print("Testing: best score of alignment for GGATGTCCTGATCCAACATCGAGGTCGTAAACCCTATTGTTGA and TCCAACATC is:")
    #print(overlap("GGATGTCCTGATCCAACATCGAGGTCGTAAACCCTATTGTTGA", "TCCAACATC"))

    print("Testing: best score of alignment for ACAAGTCCAAATTTTTGG and GGGGG is:")
    print(overlap("ACAAGTCCAAATTTTTGG", "GGGGG"))

if __name__ == "__main__":
    main()