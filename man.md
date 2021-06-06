COMPRESS(1)																User Commands												COMPRESS(1)

## NAME

​					DNA compress - compress DNA sequence files in three ways

## SYNOPSIS

​					**compress./** [FILE] [OPTION]

## DESCRIPTION

​					Get the compressed DNA files.

​					**algorithm 1**        	   combine binary compression with  segments compression

​					**algorithm 2**   			create an offline dictionary which contains all such repeats along with the details of mismatches

​					**huffman coding**       create a binary tree and its nodes based on weight 

## OPTIONS

​					**./compress.cc**

​						If the FASTA file is input, multiple files containing headers, lowercases, non-ATGC, 						subsegments, sequences will be output.  The files of subsegments and binary coded 						sequences represent main compressed files.

​				    **./alg2compress.c**

​						Output an offline dictionary which consists of "Extended seed","Type of repeat","Position of 						repeat","Length" and "Mismatch details".

​					**./Huffman_encode.c**

​						Each base is coded by 1-4 figures(0 or 1). The sequence is converted into these codes.

​					**./Huffman_decode.c**	

​						According to the coding of each base, the sequence is restored. 

## COPYRIGHT

​					Copyright  © 2021 Free Software Foundation, Inc.  License GPLv3+: GNU GPL version 3 or later 					<http://gnu.org/licenses/gpl.html>.
​					This is free software: you are free to change and redistribute it. 

​					There  is NO WARRANTY, to the extent permitted by law.

## SEE ALSO

​					You can get more details about algorithms from the websites below：

​					[SeqCompress: An algorithm for biological sequence compression - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0888754314001499?via%3Dihub)

​					[An Optimal Seed Based Compression Algorithm for DNA Sequences (hindawi.com)](https://www.hindawi.com/journals/abi/2016/3528406/)

