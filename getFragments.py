# BioInformatics Project 3
# 11/26/12

from readfasta import readfasta

# A program to get the fragments out of the file
def main():
    sumlen = 0
    fragments = []
    
    f = open('fragments.txt', 'rU')
    
    # Read in each line in the file
    for line in f:
        # Only get the sequences, not fragment # or newline char
        if(line[0] != '>' and line != '\n'):
            fragments.append(line.rstrip())

    # Print out the fragments
    for frag in fragments:
        print(frag)
        print(len(frag))
        print(" ")
        # Count the lengths
        sumlen = sumlen + len(frag)


    # Prints the length of the two sequences
    # The length is 996
    print("Each sequence has a length of",(sumlen/4))





    # Find and print the possible matching fragments for each fragment
    while(len(fragments) != 0):
        # Find the lowest length fragment
        smallFrag = min(fragments)
        fragments.remove(smallFrag)

        # Find the complement of the smallest fragment
        fragToFind = revcomp(smallFrag)

        # Find what fragments match it
        matchFrags = []
        for frag in fragments:
            if(fragToFind in frag):
                matchFrags.append(frag)

        # If there are possible matches, print them out
        if(len(matchFrags) != 0):
            print(" ")
            print("Fragment: ",smallFrag)
            print("Matching part: ",fragToFind)
            print("Possible Matches in: ")
            for match in matchFrags:
                print(match)





# A function to get the reverse complement of a fragment
# @param fragment The fragment to find the complement of
# @return comp_strand The complement of fragment
def revcomp(fragment):
    sequence = fragment

    #create a dictionary for complementation
    comp_values = {'A':'T','C':'G','G':'C','T':'A'}

    #make a variable to store the complement
    comp_strand = ''

    #loop through the original strand and complement each letter
    for letter in sequence:
        #get the complement of the letter
        newletter = comp_values[letter]
    
        #append the new letter onto the growing complement strand
        comp_strand = comp_strand + newletter

    return comp_strand




main()
