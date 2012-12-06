rFile = open("alignment_fixed_by_hand_betterer.txt")
wFile = open("alignment_fixed_by_hand.csv", "w")

for line in rFile:
    for char in line:
        wFile.write(char + ",")
    wFile.write("\n")
