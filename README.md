# Bioinformatics-Multiple-sequence-allignment
I have written Multiple sequence allignment in python. Program alligns the new sequence
with the other sequences according to the match, mismatch and gap scores provided by the user input. profile()
function creates the profile of multiple alligned sequences and returns the matrix. scoring() function calculates 
the score according to the coloumn, the letter and the profile. allignment() function creates the backtrack table using
profile and scoring, then this backtrack table is used for alligning the new sequence. alligned sequence and the old sequence
is written to a new txt file.
