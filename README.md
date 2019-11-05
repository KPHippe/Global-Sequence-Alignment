Kyle Hippe
August 19 2019
300518740

This program takes two randomly generated DNA sequences and uses a dynamic global alignment algorithm
to determine the optimal alignment of them. 

To run the test mode, set the test flag at the top of the program to True. 
This will run through worst case randomly generated sequences (both of n length since that is the worst case)
and determine the optimal alignment. Each length n is performed 10 times and the time it takes to align is recorded
and averaged and then written to a CSV file named GeneAlignmentData where data can be viewed and analyzed. 

If you would like to print both sequences, their optimal alignment, and their score, set the printing flad to true. 

If you would like the debug information, which includes the scoring matrix, parents of each cell, and backtracing steps, 
set the debug flag to True. 

If you would not like to run the test mode, leave the test flag false and input your own sequences into seq1 and seq2 variables. 
(I would like to add a command line arguement option in the future but ran out of time) If not, there are two sequences 
pre-loaded which are the two sequences given in the handout. 

In order to run, python 3 is required on the machine. 

Once in console, type 'python alignment.py' and the program will execute. 