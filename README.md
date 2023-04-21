# ProNA
Machine-learning based model for predicting interactions between bacterial sRNAs and proteins. ProNA accepts two different fasta file, one for protein sequences and the other for RNA sequences. The output is a csv file which has all the possible pairings taken from two input fasta files.
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
'''
>2IY5-T
GCCGAGGTAGCTCAGTTGGTAGAGCATGCGACTGAAAATCGCAGTGTCCGCGGTTCGATTCCGCGCCTCGGCACCA
>1TRN-R
AATCCATTGCACTCCGGATTT

'''
Also, protein fasta file with header looks like this:
'''
>3CW1-O
MKLVRFLMKLSHETVTIELKNGTQVHGTITGVDVSMNTHLKAVKMTLKNREPVQLETLSIRGNNIRYFILPDSLPLDTLLVDVEPKVKSKKREAVAGRGRGRGRGRGRGRGRGRGGPRR
>3CW1-N
MKLVRFLMKLSHETVTIELKNGTQVHGTITGVDVSMNTHLKAVKMTLKNREPVQLETLSIRGNNIRYFILPDSLPLDTLLVDVEPKVKSKKREAVAGRGRGRGRGRGRGRGRGRGGPRR
'''
