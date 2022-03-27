====================================================================
		    
 kngMap                             2022/03/28
 ver1.0.0                           Baoji, China
		    	 	       
====================================================================

System Requirements
===================

Software Requirements
---------------------
The kngMap is supported on the following operating systems 
   
   Linux
   Windows 7
   Windows 10 

Hardware Requirements/Recommendations
-------------------------------------
   Intel Core i3 2100 or later for Windows 
   Color display capable of 1024 X 768 pixel resolution
   

Memory Requirements/Recommendations
-------------------------------------
kngMap  requires approximately
800 M of available disk space on the drive or partition.

kngMap requires a minimum of
16G of RAM for mapping sequence datasets.


User's Guide
=================

Input Sequences Requirements
----------------------------
Input file requires fasta or fastq format

Execute Step
------------
Step 1: compile codes using make command
Step 2: ./locate or ./kngmap then come out help information and options reuired


A running example
-----------
./locate genome.fa reads.fa >pos.txt
./kngmap -g genome.fa -r reads.fa -p pos.txt -n reads_num -o kngmap_aligned.txt

The final mapping result is kngmap_aligned.txt


Copyright Notice
===================
Software, documentation and related materials:
Copyright (c) 2022-2025 
Baoji University of Arts and Sciences,China
All rights reserved.
