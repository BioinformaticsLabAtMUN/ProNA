# ProNA
A machine-learning based model for predicting interactions between bacterial sRNAs and proteins. ProNA takes two input fasta file: one with protein sequences and another with RNA sequences. RNA sequences must have only A, C, G, and T. The output is a csv file which the probabilities of interaction and not-interaction for all possible pairings from the two input fasta files (i.e., if there are two protein sequences and 10 sRNAs sequences then the output file will have 20 rows).
## Setup
For running ProNA, you can download Anaconda from  [here](https://www.anaconda.com/download/)
## Requirements
* Python 3
* os
* Collections
* Numpy
* Math
* Random
* Pandas
* Joblib
## Usage
To run ProNA, you need to first enter sRNA fasta file and then protein fasta file. For example sRNA fasta file with two sequences might look like this:
```
>2IY5-T
GCCGAGGTAGCTCAGTTGGTAGAGCATGCGACTGAAAATCGCAGTGTCCGCGGTTCGATTCCGCGCCTCGGCACCA
>1TRN-R
AATCCATTGCACTCCGGATTT
```
and a protein fasta file with two sequences migth look like this:
```
>3CW1-O
MKLVRFLMKLSHETVTIELKNGTQVHGTITGVDVSMNTHLKAVKMTLKNREPVQLETLSIRGNNIRYFILPDSLPLDTLLVDVEPKVKSKKREAVAGRGRGRGRGRGRGRGRGRGGPRR
>3CW1-N
MKLVRFLMKLSHETVTIELKNGTQVHGTITGVDVSMNTHLKAVKMTLKNREPVQLETLSIRGNNIRYFILPDSLPLDTLLVDVEPKVKSKKREAVAGRGRGRGRGRGRGRGRGRGGPRR
```
## Result
ProNA result will consists of:
* Rows (4 rows): sRNA ID, protein ID, Probability of interaction, and Propability of non-interaction. 
* Columns: RNA sequences times protein sequences.

**Note: The output is sorted based on probability of interaction.**
Here there are few sample rows from ProNA output:
```
2IY5-T  3CW1-O  0.9571418 0.0428581
2IY5-T  3CW1-N  0.6550656 0.3449343
1TRN-R  3CW1-O  0.4539958 0.5460042
1TRN-R  3CW1-N  0.1448176 0.8551823
```






