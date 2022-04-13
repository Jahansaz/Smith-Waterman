# Smith-Waterman
This function will align two sequences based on Smith-Waterman algorithm with an afine gap penalty.


Directions for running Smith-Waterman algorithm “SW_fun.m” from Terminal:

1.	First you need to export your MATLAB directory path to the terminal path. 
a.	For example, my MATLAB version is R2021a, and the path is: /Applications/MATLAB_R2021b.app/bin/matlab

b.	WRITE THIS COMMAND:
  export PATH=/Applications/MATLAB_R2021b.app/bin/:$PATH

2.	Now you can cd to the directory that your files are:
a.	In my case my files are in Desktop
i.	Input txt file       ==   input-sample1.txt
ii.	Similarity matrix    ==   blosum62.txt
iii.	Matlab function    ==   SW_fun.m

3.	WRITE THE FOLLOWING COMMAND:

matlab -nodisplay -r "SW_fun('sample-input1.txt','blosum62.txt')"

OR COMMOND BELOW IF YOU WANT TO SPECIFY THE OPEN GAP AND EXTENSION GAP:

matlab -nodisplay -r "SW_fun('sample-input1.txt','blosum62.txt',2,1)"


The result would be exported in your directory (in this case: Desktop) with the name of output.txt


Directions for running Smith-Waterman algorithm “SW_fun.m” from MATLAB user interface:

1.	First open your MATLAB (it should be version 2019a or later version like 2020 or 2021). 
2.	Put all the files (input file, similarity matrix, and the SW_fun.m function) in the same directory. Let’s say Desktop
3.	Open a new script in MATLAB and type one of the following code lines:

SW_fun(input_file, similarity_matrix, openGap, extGap)
SW_fun('input.txt', 'blosum62.txt', 2, 1)
SW_fun('input.txt', 'blosum62.txt')

The result would be exported in your directory (in this case: Desktop) with the name of output.txt

