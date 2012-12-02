'''
Jon Beck -- Two routines to use to read a fasta file
'''
'''
parseHeader - split out the label from the header line
Parameter: a string starting with ">" and ending without a newline
Return: the first item in the string, after the ">", up to the first space
'''
def parseHeaderLine(line):
    header = line[1:]

    label = line[1:].split(' ')[0]
    return label

'''
readfasta - the subroutine that reads the fasta file
Parameter: a filename that must be in fasta format.  The file is assumed to have:
1. arbitrary blank lines, but every line (especially including the last)
   is terminated by a line terminator (carriage return)
2. no line has only spaces on it
3. a header line as the first line
Return: a list of lists. Each inner list will have three elements:
1. the sequence identifier, the characters between the leading ">"
   and the first space
2. the entire header, the entire first line not including the ">"
3. the sequence, a single string of all the letters with no line terminators
'''
def readfasta(filename):
    resultList = []
    infile = open(filename, 'rU')

    # process the first line, which must be a header line
    line = infile.readline()
    headerLine = line
    label = parseHeaderLine(headerLine)

    # initialize the sequence accumulator
    sequence = ''

    # process all the rest of the lines in the file
    for line in infile:
        line = line.rstrip()

        # ignore blank lines
        if line == '':
            continue

        # if it's a header line, finish the previous sequence
        # and start a new one
        if line[0] == '>':
            resultList.append([label, headerLine, sequence])

            label = parseHeaderLine(line)
            headerLine = line[1:]
            sequence = ''
            
        # if we're here, we must be in letters of the sequence
        else:
            sequence += line
            
    # we're done, so clean up, terminate the last sequence, and return
    infile.close()
    resultList.append([label, headerLine, sequence])
    return resultList
