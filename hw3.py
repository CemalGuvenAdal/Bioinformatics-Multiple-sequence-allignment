import argparse
from typing import Match

parser = argparse.ArgumentParser()
parser.add_argument("--fasta")
parser.add_argument("--aln")
parser.add_argument("--out")
parser.add_argument("--match")
parser.add_argument("--mismatch")
parser.add_argument("--gap")

args = parser.parse_args()
file1Name = args.fasta
file2Name = args.aln
txtfilename = args.out
global match
match = int(args.match)
global mismatch
mismatch = int(args.mismatch)
global gap
gap = int(args.gap)

file1 = open(file1Name, 'r')
Lines = file1.readlines()
T = ""
count = 0
for line in Lines[1:]:
    count += 1
    T = T+line.strip()
file2Name = 'seqs.aln'


with open(file2Name) as f:
    lines = f.readlines()

seqlist = []
sequences = []
for line in lines:
    seqlist = line.split()
    sequences.append(seqlist[1])


def profile(sequence):
    numrows = len(sequence)
    numcol = len(sequence[0])
    rows, cols = (5, numcol)
    profil = [[0 for i in range(cols)] for j in range(rows)]

    for i in range(0, numcol):
        Ascore = 0
        Cscore = 0
        Gscore = 0
        Tscore = 0
        emptyscore = 0
        for j in range(0, numrows):
            if(sequence[j][i] == 'A'):
                Ascore += 1
            elif(sequence[j][i] == 'C'):
                Cscore += 1
            elif(sequence[j][i] == 'G'):
                Gscore += 1
            elif(sequence[j][i] == 'T'):
                Tscore += 1
            elif(sequence[j][i] == '-'):
                emptyscore += 1
        profil[0][i] = Ascore/numrows
        profil[1][i] = Cscore/numrows
        profil[2][i] = Gscore/numrows
        profil[3][i] = Tscore/numrows
        profil[4][i] = emptyscore/numrows
    return profil


def scoring(x, y, profil):
    score = 0
    if(x == 'A'):
        for i in range(0, len(profil)):
            if i == 0:
                score = score + match*profil[0][y]
            elif i == 4:
                score = score+gap*profil[4][y]
            else:
                score = score+mismatch*profil[i][y]
    if(x == 'C'):
        for i in range(0, len(profil)):

            if i == 1:
                score = score+match*profil[1][y]
            elif i == 4:
                score = score+gap*profil[4][y]
            else:
                score = score+mismatch*profil[i][y]

    if(x == 'G'):
        for i in range(0, len(profil)):
            if i == 2:
                score = score+match*profil[2][y]
            elif i == 4:
                score = score+gap*profil[4][y]
            else:
                score = score+mismatch*profil[i][y]
    if(x == 'T'):
        for i in range(0, len(profil)):
            if i == 3:
                score = score+match*profil[3][y]
            elif i == 4:
                score = score+gap*profil[4][y]
            else:
                score = score+mismatch*profil[i][y]

    return score


def allignment(profil, newseq):
    numrows = len(profil)+1
    numcol = len(profil[0])+1
    length = len(newseq)+1
    table = [[0 for i in range(numcol)] for j in range(length)]
    tablerow = len(table)
    tablecol = len(table[0])
    backtracetable = [[0 for i in range(numcol)] for j in range(length)]
    for j in range(1, tablecol):
        table[0][j] = table[0][j-1]+gap

    for j in range(1, tablerow):
        table[j][0] = table[j-1][0]+gap
    for i in range(1, tablerow):
        for j in range(1, tablecol):
            a = table[i-1][j-1]+scoring(newseq[i-1], j-1, profil)
            b = table[i-1][j]+gap
            c = table[i][j-1]+gap  # gap
            res = max(a, b, c)
            table[i][j] = res
            if(res == a):
                backtracetable[i][j] = 'd'
            elif(res == b):
                backtracetable[i][j] = 'u'
            else:
                backtracetable[i][j] = 'l'
    rows, cols = len(backtracetable), len(backtracetable[0])

    i = rows-1
    j = cols-1
    arr = []
    count = len(newseq)-1
    while(i > 0):
        while(j > 0):
            if(backtracetable[i][j] == 'd'):

                arr.insert(0, newseq[count])
                count = count-1
                i = i-1
                j = j-1
            elif(backtracetable[i][j] == 'u'):
                arr.insert(0, '-')
                i = i-1
            elif(backtracetable[i][j] == 'l'):
                arr.insert(0, '-')
                j = j-1
    return arr


def listtostr(l):
    str = ""
    for k in l:
        str = str+k
    return str


prof = profile(sequences)
result = allignment(prof, T)
strres = listtostr(result)


f = open(txtfilename, "w+")
for i in range(len(sequences)):
    f.write("sequence" + str(i+1)+" "+(sequences[i]))
    f.write('\n')
f.write("sequence  "+(strres))
f.close()
