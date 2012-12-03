Readme
---------------

This is a program to perform haplotype assembly, with some help from humans. 

To use the program, do the following:

1. Run `$ python3 team_3_scoreAlignments.py` from the command line. This will output two files in your working directory called `fragments.txt` and `overlap.txt`, which will be used implicitly in the following step. 
2. Run `$ python2 team_3_tsp.py` from the command line. Note that this requires **Python 2** due to the fact that the library which runs the max. TSP algorithm has not been updated for Python 3. This program will write a file called `alignmentOrder.txt` to your working directory, to be used in the following step.
3. Run `$ python3 team_3_prettyPrintTSPAlignments.py` from the command line. This will output the best alignment (as determined in the previous step) to the screen, as well as to a file called `alignments.csv` which you can open in Excel.
