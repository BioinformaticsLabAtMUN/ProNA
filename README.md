# ProNA
Machine-learning based model for predicting interactions between bacterial sRNAs and proteins. ProNA accepts two different fasta file, one for protein sequences and the other for RNA sequences. The output is a csv file which has all the possible pairings taken from two input fasta files.
## Setup
To run this script, you can download Anaconda from  [here](https://www.anaconda.com/download/)
## Necessaries
* Python 3
* os
* Collections
* Numpy
* Math
* Random
* Pandas
* Joblib
## Usage
To run the script, you need to first enter sRNA fasta file and then protein fasta file. For example sRNA fasta file with two sequences would look like:
```
>2IY5-T
GCCGAGGTAGCTCAGTTGGTAGAGCATGCGACTGAAAATCGCAGTGTCCGCGGTTCGATTCCGCGCCTCGGCACCA
>1TRN-R
AATCCATTGCACTCCGGATTT
```
Also, protein fasta file with header looks like this:
```
>3CW1-O
MKLVRFLMKLSHETVTIELKNGTQVHGTITGVDVSMNTHLKAVKMTLKNREPVQLETLSIRGNNIRYFILPDSLPLDTLLVDVEPKVKSKKREAVAGRGRGRGRGRGRGRGRGRGGPRR
>3CW1-N
MKLVRFLMKLSHETVTIELKNGTQVHGTITGVDVSMNTHLKAVKMTLKNREPVQLETLSIRGNNIRYFILPDSLPLDTLLVDVEPKVKSKKREAVAGRGRGRGRGRGRGRGRGRGGPRR
```
## Result
ProNA result will be saved in a dataframe below characteristics:
* Rows (4 rows): sRNA ID, protein ID, Propability of interaction, and Propability of non-interaction. 
* Columns: RNA sequences times protein sequences.
**Note: The output is sorted based on probability of interaction.
