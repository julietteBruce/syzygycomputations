Here are what the files do.


------ relevantRange.m2 ---------
computes the relevant range given D={d1,d2} and B={b1,b2}. Currently hardcoded for P1xP1.
main function: relevantRange(0,B,D).



------ step1.py -----------------
computes the betti data for (p,q) in the relevant range w/o condor. run from top directory as:

python3 src/step1.py output_dir d1 d2 b1 b2 char

where output_dir is a directory to store the matricies and betti data, D={d1, d2}, B={b1,b2}, and char is the characteristic of the field to compute matrices. 


------ step1a_condor.py ---------
makes matrices and submits jobs to condor (with 2 GB of RAM). used to compute betti data for (p,q) in the relevant range. run from top directory as:

python3 src/step1a_condor.py output_dir d1 d2 b1 b2 char

where output_dir is a directory to store the matricies and betti data, D={d1, d2}, B={b1,b2}, and char is the characteristic of the field to compute matrices.


------ step1b_condor.py ---------
determines which ranks that need to be computed and asks user to input memory to request, then resubmits jobs to condor. if all ranks are computed, then the betti data is computed. 

python3 src/step1a_condor.py output_dir d1 d2 b1 b2 char

where output_dir is a directory to store the matricies and betti data, D={d1, d2}, B={b1,b2}, and char is the characteristic of the field to compute matrices.



------ ranks.magma -------------
this is the magma script that computes the rank of the input matrix. run like this:

magma -b p:= matrixFile:= ranks.magma

input the caracteristic for p and the matrix file for matrixFile



