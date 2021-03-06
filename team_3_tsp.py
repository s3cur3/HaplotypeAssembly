"""
A program to take the scored overlap between all sequence fragments and perform
the maximum Traveling Salesperson Problem on them. Outputs a file called
"alignmentOrder.txt" which is the best order to align the sequences in.

Relies on files output by team_3_scoreAlignments.py. First run that program, then
run this one in the same directory.

Tyler Young
Written for Python 2 (as pyevolve has not been updated
    for compatibility with Python3)
"""


# The following comment is from the original code and I'm too lazy
# to turn it into coherent English, but it seems wrong to delete it.
# There is nothing especially interesting at the site mentioned.

# The follow TSP routines was get from the above site, I'm too lazy to reinvent a new pretty wheel:
# http://www.psychicorigami.com/2007/04/17/tackling-the-travelling-salesman-problem-part-one/
# Routines:
# - cartesian_matrix
# - read_coords
# - tour_length
# - write_tour_to_img

from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import GAllele
from pyevolve import Mutators
from pyevolve import Initializators
from pyevolve import DBAdapters
from pyevolve import Crossovers
from pyevolve import Consts
import datetime, random
from math import sqrt
from team_3_scoreAlignments import overlap

PIL_SUPPORT = False

def read_coords(coord_file):
    """ Read the coords from file """
    coords=[]
    data = (coord_file.read()).split("\n")
    for row in data:
        theRow = []
        vals = (row.rstrip()).split()

        if len(vals) == 0:
            continue

        for val in vals:
            theRow.append(int(val))
        coords.append(theRow)
    return coords

def readFragmentFile(fragmentFile):
    allLines = (fragmentFile.read()).split("\n")

    fragments = []
    for line in allLines:
        if len(line.rstrip()) != 0:
            fragments.append(line)
    return fragments

memo = {}
def getAlignmentScore(matrix, tour):
    """ Returns the total score for this solution """
    global fragments, memo

    alignmentList = []
    for i in range( len(tour)-1 ):
        alignmentList.append( ( int(tour[i]), int(tour[i+1])) )

    offset = 0
    for pair in alignmentList:
        if pair[0] == pair[1] == 0:
            print("Error, error, error!")
            break
        if pair in memo:
            offset += memo[pair]
        else:
            theOverlap, addToOffset = overlap(fragments[pair[0]], fragments[pair[1]])
            offset += addToOffset
            memo[pair] = addToOffset

    # Return the length of the overall alignment
    # Offset tells us where the last thing lines up against the whole sequence
    # The full length of the aligned fragments is the offset plus the length of
    # the last fragment
    return offset + len(fragments[ tour[-1] ])

def G1DListTSPInitializator(genome, **args):
   """ The initializator for the TSP """
   genome.clearList()
   # Note: getListSize() used to be just a field called listSize
   lst = [i for i in xrange(genome.getListSize())]

   for i in xrange(genome.getListSize()):
      choice = random.choice(lst)
      lst.remove(choice)
      genome.append(choice)

cm = []
coords = []
fragments = []

def eval_func(chromosome):
   """ The evaluation function """
   global cm
   # return the inverse of the alignment score so that we do max. TSP
   return getAlignmentScore(cm, chromosome)


def main_run(distancesFileName, fragmentFileName, crossover_rate=1.0, mutation_rate=0.03, population_size=80):
    """
    @param distancesFileName String The file containing the pairwise distances of all
                                fragments
    """
    global cm, coords, fragments

    # Load the fragments
    fragmentFile = open(fragmentFileName, "r")
    fragments = readFragmentFile(fragmentFile)

    # load the tsp data file
    filehandle = open(distancesFileName, "r")
    coords = read_coords(filehandle)
    cm = coords

    # set the alleles to the cities numbers
    setOfAlleles = GAllele.GAlleles(homogeneous=True)
    lst = [ i for i in xrange(len(coords)) ]
    a = GAllele.GAlleleList(lst)
    setOfAlleles.add(a)

    genome = G1DList.G1DList(len(coords))
    genome.setParams(allele=setOfAlleles)

    genome.evaluator.set(eval_func)
    genome.mutator.set(Mutators.G1DListMutatorSwap)
    genome.crossover.set(Crossovers.G1DListCrossoverOX)
    genome.initializator.set(G1DListTSPInitializator)

    ga = GSimpleGA.GSimpleGA(genome)
    ga.setGenerations(5000) # 10000 is a good "real" value
    ga.setMinimax(Consts.minimaxType["minimize"])
    ga.setCrossoverRate(crossover_rate)
    ga.setMutationRate(mutation_rate)
    ga.setPopulationSize(population_size)

    ga.evolve(freq_stats=100)
    best = ga.bestIndividual()

    return eval_func(best), best.getInternalList()

def writeResults(fileNameNoExt, sequence, score, crossover, mutation, popSize):
    # Write the best sequence out to a file
    theFile = open(fileNameNoExt + ".txt", "w")
    for num in sequence:
        theFile.write(str(num))
        theFile.write(" ")

    # Write a second file in case we're running multiple tests at once
    theFile = open(fileNameNoExt + "At" + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ".txt", "w")
    for num in sequence:
        theFile.write(str(num) + " ")

    theFile.write("\n# Parameters were: Pop size: ")
    theFile.write(str(popSize))
    theFile.write(", Crossover rate: ")
    theFile.write(str(crossover))

    theFile.write(", Mutation rate: ")
    theFile.write(str(mutation))

    theFile.write("\n# Score was: ")
    theFile.write(str(score))



if __name__ == "__main__":
    bestScore = -1
    bestSequence = []
    bestCrossover = -1
    bestMutationRate = -1
    bestPopSize = -1

    # Repeat many times so that we vary the parameters of the model (results in much
    # better... results.)
    for i in range(1):
        crossover_rate = float(random.randrange(20, 100, 5))/100
        mutation_rate = float(random.randrange(0, 15, 1))/100
        population_size = random.randrange(10, 150, 10)
        score, sequence = main_run("overlap.txt", "fragments.txt", crossover_rate, mutation_rate, population_size)

        if score > bestScore:
            bestScore = score
            bestSequence = sequence

            # Remember the parameters that gave us this score
            bestCrossover = crossover_rate
            bestMutationRate =mutation_rate
            bestPopSize = population_size

        writeResults("bestAlignmentOrderSoFar", bestSequence, bestScore, bestCrossover, bestMutationRate, bestPopSize)

    print "Best score: ", bestScore
    print "Sequence: ", bestSequence

    writeResults("alignmentOrder", bestSequence, bestScore, bestCrossover, bestMutationRate, bestPopSize)
