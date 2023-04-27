#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
import re
from collections import Counter
import numpy as np
import math
import randomh
import pandas as pd
import joblib


class Sequence(object):
    def __init__(self, file):
        self.file = file  
        self.fasta_list = []  
        self.sequence_number = 0  
        self.sequence_type = ''  
        self.is_equal = False  
        self.minimum_length = 1  
        self.maximum_length = 0 
        self.minimum_length_without_minus = 1  
        self.maximum_length_without_minus = 0  
        self.error_msg = '' 

        self.fasta_list, self.error_msg = self.read_fasta(self.file)
        self.sequence_number = len(self.fasta_list)

        if self.sequence_number > 0:
            self.is_equal, self.minimum_length, self.maximum_length, self.minimum_length_without_minus, self.maximum_length_without_minus = self.sequence_with_equal_length()
            self.sequence_type = self.check_sequence_type()
        else:
            self.error_msg = 'File format error.'

    def read_fasta(self, file):
        global string_list
        global identifier
        global sequence
        global name
        string_list=[]
        identifier=pd.DataFrame()
        
        msg=""
        if not os.path.exists(self.file):
            msg = 'Error: file %s does not exist.' % self.file
            return [], None, msg
        with open(file) as f:
            records = f.read()
        records = records.split('>')[1:]
        fasta_sequences = []
        for fasta in records:
            array = fasta.split('\n')
            header= array[0]           
            sequence=array[1]
            sequence = re.sub('[^ACDEFGHIKLMNPQRSTUVWY-]', '-',sequence.upper())          
            name=header
            label = '0'
            label_train = "testing"
            fasta_sequences.append([name, sequence, label, label_train])    
            string_list.append(name)             
            identifier = pd.DataFrame(string_list,index=None)
        return fasta_sequences,msg
        

    def sequence_with_equal_length(self):
        length_set = set()
        length_set_1 = set()
        for item in self.fasta_list:
            length_set.add(len(item[1]))
            length_set_1.add(len(re.sub('-', '', item[1])))

        length_set = sorted(length_set)
        length_set_1 = sorted(length_set_1)
        if len(length_set) == 1:
            return True, length_set[0], length_set[-1], length_set_1[0], length_set_1[-1]
        else:
            return False, length_set[0], length_set[-1], length_set_1[0], length_set_1[-1]

    def check_sequence_type(self):
        tmp_fasta_list = []
        if len(self.fasta_list) < 100:
            tmp_fasta_list = self.fasta_list
        else:
            random_index = random.sample(range(0, len(self.fasta_list)), 100)
            for i in random_index:
                tmp_fasta_list.append(self.fasta_list[i])

        sequence = ''
        for item in tmp_fasta_list:
            sequence += item[1]

        char_set = set(sequence)
        if 5 < len(char_set) <= 21:
            for line in self.fasta_list:
                line[1] = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '-', line[1])
            return 'Protein'
        elif 0 < len(char_set) <= 5 and 'T' in char_set:
            return 'DNA'
        elif 0 < len(char_set) <= 5 and 'U' in char_set:
            for line in self.fasta_list:
                line[1] = re.sub('U', 'T', line[1])
            return 'RNA'
        else:
            return 'Unknown'


class Descriptor(Sequence):
    def __init__(self, file, kw):
        super(Descriptor, self).__init__(file=file)
        self.kw = kw  
        self.encoding_array = np.array([])  
        self.column = 0  
        self.row = 0  
        """ variable for ACC descriptors """
        self.myDiIndex = {
            'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3,
            'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
            'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11,
            'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15
        }
        self.myTriIndex = {
            'AAA': 0, 'AAC': 1, 'AAG': 2, 'AAT': 3,
            'ACA': 4, 'ACC': 5, 'ACG': 6, 'ACT': 7,
            'AGA': 8, 'AGC': 9, 'AGG': 10, 'AGT': 11,
            'ATA': 12, 'ATC': 13, 'ATG': 14, 'ATT': 15,
            'CAA': 16, 'CAC': 17, 'CAG': 18, 'CAT': 19,
            'CCA': 20, 'CCC': 21, 'CCG': 22, 'CCT': 23,
            'CGA': 24, 'CGC': 25, 'CGG': 26, 'CGT': 27,
            'CTA': 28, 'CTC': 29, 'CTG': 30, 'CTT': 31,
            'GAA': 32, 'GAC': 33, 'GAG': 34, 'GAT': 35,
            'GCA': 36, 'GCC': 37, 'GCG': 38, 'GCT': 39,
            'GGA': 40, 'GGC': 41, 'GGG': 42, 'GGT': 43,
            'GTA': 44, 'GTC': 45, 'GTG': 46, 'GTT': 47,
            'TAA': 48, 'TAC': 49, 'TAG': 50, 'TAT': 51,
            'TCA': 52, 'TCC': 53, 'TCG': 54, 'TCT': 55,
            'TGA': 56, 'TGC': 57, 'TGG': 58, 'TGT': 59,
            'TTA': 60, 'TTC': 61, 'TTG': 62, 'TTT': 63
        }

    
    #protein feature sets
    def Protein_AAC(self):
        try:
            self.encoding_array = np.array([])
            AA = 'ACDEFGHIKLMNPQRSTVWY'
            header = ['SampleName', 'label']
            encodings = []
            for i in AA:
                header.append(i)
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                count = Counter(sequence)
                for key in count:
                    count[key] = count[key] / len(sequence)
                code = [name, label]
                for aa in AA:
                    code.append(count[aa])
                encodings.append(code)

            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

   
    def Protein_DPC(self):
        try:
            self.encoding_array = np.array([])

            AA = 'ACDEFGHIKLMNPQRSTVWY'
            encodings = []
            diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
            header = ['SampleName', 'label'] + diPeptides
            encodings.append(header)

            AADict = {}
            for i in range(len(AA)):
                AADict[AA[i]] = i

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                tmpCode = [0] * 400
                for j in range(len(sequence) - 2 + 1):
                    tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j + 1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[
                        sequence[j + 1]]] + 1
                if sum(tmpCode) != 0:
                    tmpCode = [i / sum(tmpCode) for i in tmpCode]
                code = code + tmpCode
                encodings.append(code)

            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

    def Protein_DDE(self):
        try:
            self.encoding_array = np.array([])
            AA = 'ACDEFGHIKLMNPQRSTVWY'

            myCodons = {'A': 4, 'C': 2, 'D': 2, 'E': 2, 'F': 2, 'G': 4, 'H': 2, 'I': 3, 'K': 2, 'L': 6,
                        'M': 1, 'N': 2, 'P': 4, 'Q': 2, 'R': 6, 'S': 6, 'T': 4, 'V': 4, 'W': 1, 'Y': 2
                        }

            encodings = []
            diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
            header = ['SampleName', 'label'] + diPeptides
            encodings.append(header)

            myTM = []
            for pair in diPeptides:
                myTM.append((myCodons[pair[0]] / 61) * (myCodons[pair[1]] / 61))

            AADict = {}
            for i in range(len(AA)):
                AADict[AA[i]] = i

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                tmpCode = [0] * 400
                for j in range(len(sequence) - 2 + 1):
                    tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j + 1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[
                        sequence[j + 1]]] + 1
                if sum(tmpCode) != 0:
                    tmpCode = [i / sum(tmpCode) for i in tmpCode]

                myTV = []
                for j in range(len(myTM)):
                    myTV.append(myTM[j] * (1 - myTM[j]) / (len(sequence) - 1))

                for j in range(len(tmpCode)):
                    tmpCode[j] = (tmpCode[j] - myTM[j]) / math.sqrt(myTV[j])

                code = code + tmpCode
                encodings.append(code)

            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False


        
    def Protein_GAAC(self):
        try:
            self.encoding_array = np.array([])

            group = {
                'alphatic': 'GAVLMI',
                'aromatic': 'FYW',
                'postivecharge': 'KRH',
                'negativecharge': 'DE',
                'uncharge': 'STCPNQ'
            }

            groupKey = group.keys()

            encodings = []
            header = ['SampleName', 'label']
            for key in groupKey:
                header.append(key)
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                count = Counter(sequence)
                myDict = {}
                for key in groupKey:
                    for aa in group[key]:
                        myDict[key] = myDict.get(key, 0) + count[aa]

                for key in groupKey:
                    code.append(myDict[key] / len(sequence))
                encodings.append(code)

            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

 
    def Protein_GDPC(self):
        try:
            self.encoding_array = np.array([])

            group = {
                'alphaticr': 'GAVLMI',
                'aromatic': 'FYW',
                'postivecharger': 'KRH',
                'negativecharger': 'DE',
                'uncharger': 'STCPNQ'
            }

            groupKey = group.keys()
            baseNum = len(groupKey)
            dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]

            index = {}
            for key in groupKey:
                for aa in group[key]:
                    index[aa] = key

            encodings = []
            header = ['SampleName', 'label'] + dipeptide
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]

                code = [name, label]
                myDict = {}
                for t in dipeptide:
                    myDict[t] = 0

                sum = 0
                for j in range(len(sequence) - 2 + 1):
                    myDict[index[sequence[j]] + '.' + index[sequence[j + 1]]] = myDict[index[sequence[j]] + '.' + index[
                        sequence[j + 1]]] + 1
                    sum = sum + 1

                if sum == 0:
                    for t in dipeptide:
                        code.append(0)
                else:
                    for t in dipeptide:
                        code.append(myDict[t] / sum)
                encodings.append(code)

            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

    def Protein_GTPC(self):
        try:
            self.encoding_array = np.array([])

            group = {
                'alphaticr': 'GAVLMI',
                'aromatic': 'FYW',
                'postivecharger': 'KRH',
                'negativecharger': 'DE',
                'uncharger': 'STCPNQ'
            }

            groupKey = group.keys()
            baseNum = len(groupKey)
            triple = [g1 + '.' + g2 + '.' + g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

            index = {}
            for key in groupKey:
                for aa in group[key]:
                    index[aa] = key

            encodings = []
            header = ['SampleName', 'label'] + triple
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]

                code = [name, label]
                myDict = {}
                for t in triple:
                    myDict[t] = 0

                sum = 0
                for j in range(len(sequence) - 3 + 1):
                    myDict[index[sequence[j]] + '.' + index[sequence[j + 1]] + '.' + index[sequence[j + 2]]] = myDict[index[
                                                                                                                        sequence[
                                                                                                                            j]] + '.' +
                                                                                                                    index[
                                                                                                                        sequence[
                                                                                                                            j + 1]] + '.' +
                                                                                                                    index[
                                                                                                                        sequence[
                                                                                                                            j + 2]]] + 1
                    sum = sum + 1

                if sum == 0:
                    for t in triple:
                        code.append(0)
                else:
                    for t in triple:
                        code.append(myDict[t] / sum)
                encodings.append(code)

            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

    def Protein_CTDC(self):
        try:
            group1 = {
                'hydrophobicity_PRAM900101': 'RKEDQN',
                'hydrophobicity_ARGP820101': 'QSTNGDE',
                'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
                'hydrophobicity_PONP930101': 'KPDESNQT',
                'hydrophobicity_CASG920101': 'KDEQPSRNTG',
                'hydrophobicity_ENGD860101': 'RDKENQHYP',
                'hydrophobicity_FASG890101': 'KERSQD',
                'normwaalsvolume': 'GASTPDC',
                'polarity': 'LIFWCMVY',
                'polarizability': 'GASDT',
                'charge': 'KR',
                'secondarystruct': 'EALMQKRH',
                'solventaccess': 'ALFCGIVW'
            }
            group2 = {
                'hydrophobicity_PRAM900101': 'GASTPHY',
                'hydrophobicity_ARGP820101': 'RAHCKMV',
                'hydrophobicity_ZIMJ680101': 'HMCKV',
                'hydrophobicity_PONP930101': 'GRHA',
                'hydrophobicity_CASG920101': 'AHYMLV',
                'hydrophobicity_ENGD860101': 'SGTAW',
                'hydrophobicity_FASG890101': 'NTPG',
                'normwaalsvolume': 'NVEQIL',
                'polarity': 'PATGS',
                'polarizability': 'CPNVEQIL',
                'charge': 'ANCQGHILMFPSTWYV',
                'secondarystruct': 'VIYCWFT',
                'solventaccess': 'RKQEND'
            }
            group3 = {
                'hydrophobicity_PRAM900101': 'CLVIMFW',
                'hydrophobicity_ARGP820101': 'LYPFIW',
                'hydrophobicity_ZIMJ680101': 'LPFYI',
                'hydrophobicity_PONP930101': 'YMFWLCVI',
                'hydrophobicity_CASG920101': 'FIWC',
                'hydrophobicity_ENGD860101': 'CVLIMF',
                'hydrophobicity_FASG890101': 'AYHWVMFLIC',
                'normwaalsvolume': 'MHKFRYW',
                'polarity': 'HQRKNED',
                'polarizability': 'KMHFRYW',
                'charge': 'DE',
                'secondarystruct': 'GNPSD',
                'solventaccess': 'MSPTHY'
            }

            groups = [group1, group2, group3]
            property = (
                'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
                'hydrophobicity_PONP930101',
                'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
                'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

            encodings = []
            header = ['SampleName', 'label']
            for p in property:
                for g in range(1, len(groups) + 1):
                    header.append(p + '.G' + str(g))
            encodings.append(header)
            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                for p in property:
                    c1 = self.Count(group1[p], sequence) / len(sequence)
                    c2 = self.Count(group2[p], sequence) / len(sequence)
                    c3 = 1 - c1 - c2
                    code = code + [c1, c2, c3]
                encodings.append(code)

            self.encoding_array = np.array([])
            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

    def Protein_CTDT(self):
        try:
            group1 = {
                'hydrophobicity_PRAM900101': 'RKEDQN',
                'hydrophobicity_ARGP820101': 'QSTNGDE',
                'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
                'hydrophobicity_PONP930101': 'KPDESNQT',
                'hydrophobicity_CASG920101': 'KDEQPSRNTG',
                'hydrophobicity_ENGD860101': 'RDKENQHYP',
                'hydrophobicity_FASG890101': 'KERSQD',
                'normwaalsvolume': 'GASTPDC',
                'polarity': 'LIFWCMVY',
                'polarizability': 'GASDT',
                'charge': 'KR',
                'secondarystruct': 'EALMQKRH',
                'solventaccess': 'ALFCGIVW'
            }
            group2 = {
                'hydrophobicity_PRAM900101': 'GASTPHY',
                'hydrophobicity_ARGP820101': 'RAHCKMV',
                'hydrophobicity_ZIMJ680101': 'HMCKV',
                'hydrophobicity_PONP930101': 'GRHA',
                'hydrophobicity_CASG920101': 'AHYMLV',
                'hydrophobicity_ENGD860101': 'SGTAW',
                'hydrophobicity_FASG890101': 'NTPG',
                'normwaalsvolume': 'NVEQIL',
                'polarity': 'PATGS',
                'polarizability': 'CPNVEQIL',
                'charge': 'ANCQGHILMFPSTWYV',
                'secondarystruct': 'VIYCWFT',
                'solventaccess': 'RKQEND'
            }
            group3 = {
                'hydrophobicity_PRAM900101': 'CLVIMFW',
                'hydrophobicity_ARGP820101': 'LYPFIW',
                'hydrophobicity_ZIMJ680101': 'LPFYI',
                'hydrophobicity_PONP930101': 'YMFWLCVI',
                'hydrophobicity_CASG920101': 'FIWC',
                'hydrophobicity_ENGD860101': 'CVLIMF',
                'hydrophobicity_FASG890101': 'AYHWVMFLIC',
                'normwaalsvolume': 'MHKFRYW',
                'polarity': 'HQRKNED',
                'polarizability': 'KMHFRYW',
                'charge': 'DE',
                'secondarystruct': 'GNPSD',
                'solventaccess': 'MSPTHY'
            }

            groups = [group1, group2, group3]
            property = (
                'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
                'hydrophobicity_PONP930101',
                'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
                'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

            encodings = []
            header = ['SampleName', 'label']
            for p in property:
                for tr in ('Tr1221', 'Tr1331', 'Tr2332'):
                    header.append(p + '.' + tr)
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
                for p in property:
                    c1221, c1331, c2332 = 0, 0, 0
                    for pair in aaPair:
                        if (pair[0] in group1[p] and pair[1] in group2[p]) or (
                                pair[0] in group2[p] and pair[1] in group1[p]):
                            c1221 = c1221 + 1
                            continue
                        if (pair[0] in group1[p] and pair[1] in group3[p]) or (
                                pair[0] in group3[p] and pair[1] in group1[p]):
                            c1331 = c1331 + 1
                            continue
                        if (pair[0] in group2[p] and pair[1] in group3[p]) or (
                                pair[0] in group3[p] and pair[1] in group2[p]):
                            c2332 = c2332 + 1
                    code = code + [c1221 / len(aaPair), c1331 / len(aaPair), c2332 / len(aaPair)]
                encodings.append(code)

            self.encoding_array = np.array([])
            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False


    def Protein_CTDD(self):
        try:
            group1 = {
                'hydrophobicity_PRAM900101': 'RKEDQN',
                'hydrophobicity_ARGP820101': 'QSTNGDE',
                'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
                'hydrophobicity_PONP930101': 'KPDESNQT',
                'hydrophobicity_CASG920101': 'KDEQPSRNTG',
                'hydrophobicity_ENGD860101': 'RDKENQHYP',
                'hydrophobicity_FASG890101': 'KERSQD',
                'normwaalsvolume': 'GASTPDC',
                'polarity': 'LIFWCMVY',
                'polarizability': 'GASDT',
                'charge': 'KR',
                'secondarystruct': 'EALMQKRH',
                'solventaccess': 'ALFCGIVW'
            }
            group2 = {
                'hydrophobicity_PRAM900101': 'GASTPHY',
                'hydrophobicity_ARGP820101': 'RAHCKMV',
                'hydrophobicity_ZIMJ680101': 'HMCKV',
                'hydrophobicity_PONP930101': 'GRHA',
                'hydrophobicity_CASG920101': 'AHYMLV',
                'hydrophobicity_ENGD860101': 'SGTAW',
                'hydrophobicity_FASG890101': 'NTPG',
                'normwaalsvolume': 'NVEQIL',
                'polarity': 'PATGS',
                'polarizability': 'CPNVEQIL',
                'charge': 'ANCQGHILMFPSTWYV',
                'secondarystruct': 'VIYCWFT',
                'solventaccess': 'RKQEND'
            }
            group3 = {
                'hydrophobicity_PRAM900101': 'CLVIMFW',
                'hydrophobicity_ARGP820101': 'LYPFIW',
                'hydrophobicity_ZIMJ680101': 'LPFYI',
                'hydrophobicity_PONP930101': 'YMFWLCVI',
                'hydrophobicity_CASG920101': 'FIWC',
                'hydrophobicity_ENGD860101': 'CVLIMF',
                'hydrophobicity_FASG890101': 'AYHWVMFLIC',
                'normwaalsvolume': 'MHKFRYW',
                'polarity': 'HQRKNED',
                'polarizability': 'KMHFRYW',
                'charge': 'DE',
                'secondarystruct': 'GNPSD',
                'solventaccess': 'MSPTHY'
            }

            groups = [group1, group2, group3]
            property = (
                'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
                'hydrophobicity_PONP930101',
                'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
                'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

            encodings = []
            header = ['SampleName', 'label']
            for p in property:
                for g in ('1', '2', '3'):
                    for d in ['0', '25', '50', '75', '100']:
                        header.append(p + '.' + g + '.residue' + d)
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                for p in property:
                    code = code + self.Count1(group1[p], sequence) + self.Count1(group2[p], sequence) + self.Count1(
                        group3[p], sequence)
                encodings.append(code)

            self.encoding_array = np.array([])
            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

  
    def Protein_ASDC(self):
        try:
            # clear
            self.encoding_array = np.array([])

            AA = 'ACDEFGHIKLMNPQRSTVWY'
            encodings = []
            aaPairs = []
            for aa1 in AA:
                for aa2 in AA:
                    aaPairs.append(aa1 + aa2)

            header = ['SampleName', 'label']
            header += [aa1 + aa2 for aa1 in AA for aa2 in AA]
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                sum = 0
                pair_dict = {}
                for pair in aaPairs:
                    pair_dict[pair] = 0
                for j in range(len(sequence)):
                    for k in range(j + 1, len(sequence)):
                        if sequence[j] in AA and sequence[k] in AA:
                            pair_dict[sequence[j] + sequence[k]] += 1
                            sum += 1
                for pair in aaPairs:
                    code.append(pair_dict[pair] / sum)
                encodings.append(code)

            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

    #RNA feature sets
    def NAC(self):
        try:
            NA = 'ACGT'
            encodings = []
            header = ['SampleName', 'label']
            for i in NA:
                header.append(i)
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                count = Counter(sequence)
                for key in count:
                    count[key] = count[key] / len(sequence)
                code = [name, label]
                for na in NA:
                    code.append(count[na])
                encodings.append(code)
            self.encoding_array = np.array([])
            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

    def DNC(self):
        try:
            base = 'ACGT'
            encodings = []
            dinucleotides = [n1 + n2 for n1 in base for n2 in base]
            header = ['SampleName', 'label'] + dinucleotides
            encodings.append(header)

            AADict = {}
            for i in range(len(base)):
                AADict[base[i]] = i

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                tmpCode = [0] * 16
                for j in range(len(sequence) - 2 + 1):
                    tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j + 1]]] = tmpCode[AADict[sequence[j]] * 4 + AADict[
                        sequence[j + 1]]] + 1
                if sum(tmpCode) != 0:
                    tmpCode = [i / sum(tmpCode) for i in tmpCode]
                
                code = code + tmpCode
                encodings.append(code)
            self.encoding_array = np.array([])
            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

    def TNC(self):
        try:
            AA = 'ACGT'
            encodings = []
            triPeptides = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
            header = ['SampleName', 'label'] + triPeptides
            encodings.append(header)

            AADict = {}
            for i in range(len(AA)):
                AADict[AA[i]] = i

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                tmpCode = [0] * 64
                for j in range(len(sequence) - 3 + 1):
                    tmpCode[AADict[sequence[j]] * 16 + AADict[sequence[j + 1]] * 4 + AADict[sequence[j + 2]]] = tmpCode[
                                                                                                                    AADict[
                                                                                                                        sequence[
                                                                                                                            j]] * 16 +
                                                                                                                    AADict[
                                                                                                                        sequence[
                                                                                                                            j + 1]] * 4 +
                                                                                                                    AADict[
                                                                                                                        sequence[
                                                                                                                            j + 2]]] + 1
                if sum(tmpCode) != 0:
                    tmpCode = [i / sum(tmpCode) for i in tmpCode]
                code = code + tmpCode
                encodings.append(code)
            self.encoding_array = np.array([])
            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

 
    def CKSNAP(self):
        try:
            gap=3
            if self.minimum_length_without_minus < gap +2:
                self.error_msg = 'CKSNAP - all the sequence length should be larger than the (gap value) + 2 = %s.' % (
                        gap + 2)
                return False

            AA = 'ACGT'
            encodings = []
            aaPairs = []
            for aa1 in AA:
                for aa2 in AA:
                    aaPairs.append(aa1 + aa2)

            header = ['SampleName', 'label']
            for g in range(gap + 1):
                for aa in aaPairs:
                    header.append(aa + '.gap' + str(g))
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], i[1], i[2]
                code = [name, label]
                for g in range(gap + 1):
                    myDict = {}
                    for pair in aaPairs:
                        myDict[pair] = 0
                    sum = 0
                    for index1 in range(len(sequence)):
                        index2 = index1 + g + 1
                        if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[
                            index2] in AA:
                            myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                            sum = sum + 1
                    for pair in aaPairs:
                        code.append(myDict[pair] / sum)
                encodings.append(code)
            self.encoding_array = np.array([])
            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

  

    def ASDC(self):
        try:
            self.encoding_array = np.array([])

            AA = 'ACGT'
            encodings = []
            aaPairs = []
            for aa1 in AA:
                for aa2 in AA:
                    aaPairs.append(aa1 + aa2)

            header = ['SampleName', 'label']
            header += [aa1 + aa2 for aa1 in AA for aa2 in AA]
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                sum = 0
                pair_dict = {}
                for pair in aaPairs:
                    pair_dict[pair] = 0
                for j in range(len(sequence)):
                    for k in range(j + 1, len(sequence)):
                        if sequence[j] in AA and sequence[k] in AA:
                            pair_dict[sequence[j] + sequence[k]] += 1
                            sum += 1
                for pair in aaPairs:
                    code.append(pair_dict[pair] / sum)
                encodings.append(code)

            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False

 
    def MMI(self):
        try:
            NA = 'ACGT'
            dinucleotide_list = [a1 + a2 for a1 in NA for a2 in NA]
            trinucleotide_list = [a1 + a2 + a3 for a1 in NA for a2 in NA for a3 in NA]
            dinucleotide_dict = {}
            trinucleotide_dict = {}
            for elem in dinucleotide_list:
                dinucleotide_dict[''.join(sorted(elem))] = 0
            for elem in trinucleotide_list:
                trinucleotide_dict[''.join(sorted(elem))] = 0

            encodings = []
            header = ['SampleName', 'label']
            header += ['MMI_%s' % elem for elem in sorted(dinucleotide_dict.keys())]
            header += ['MMI_%s' % elem for elem in sorted(trinucleotide_dict.keys())]
            encodings.append(header)

            for i in self.fasta_list:
                name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
                code = [name, label]
                f1_dict = {
                    'A': 0,
                    'C': 0,
                    'G': 0,
                    'T': 0,
                }
                f2_dict = dinucleotide_dict.copy()
                f3_dict = trinucleotide_dict.copy()

                for elem in sequence:
                    if elem in f1_dict:
                        f1_dict[elem] += 1
                for key in f1_dict:
                    f1_dict[key] /= len(sequence)

                for i in range(len(sequence) - 1):
                    if ''.join(sorted(sequence[i: i + 2])) in f2_dict:
                        f2_dict[''.join(sorted(sequence[i: i + 2]))] += 1
                for key in f2_dict:
                    f2_dict[key] /= (len(sequence) - 1)

                for i in range(len(sequence) - 2):
                    if ''.join(sorted(sequence[i: i + 3])) in f3_dict:
                        f3_dict[''.join(sorted(sequence[i: i + 3]))] += 1
                for key in f3_dict:
                    f3_dict[key] /= (len(sequence) - 2)

                for key in sorted(f2_dict.keys()):
                    if f2_dict[key] != 0 and f1_dict[key[0]] * f1_dict[key[1]] != 0:
                        code.append(f2_dict[key] * math.log(f2_dict[key] / (f1_dict[key[0]] * f1_dict[key[1]])))
                    else:
                        code.append(0)
                for key in sorted(f3_dict.keys()):
                    element_1 = 0
                    element_2 = 0
                    element_3 = 0
                    if f2_dict[key[0:2]] != 0 and f1_dict[key[0]] * f1_dict[key[1]] != 0:
                        element_1 = f2_dict[key[0:2]] * math.log(f2_dict[key[0:2]] / (f1_dict[key[0]] * f1_dict[key[1]]))
                    if f2_dict[key[0] + key[2]] != 0 and f1_dict[key[2]] != 0:
                        element_2 = (f2_dict[key[0] + key[2]] / f1_dict[key[2]]) * math.log(
                            f2_dict[key[0] + key[2]] / f1_dict[key[2]])
                    if f2_dict[key[1:3]] != 0 and f3_dict[key] / f2_dict[key[1:3]] != 0:
                        element_3 = (f3_dict[key] / f2_dict[key[1:3]]) * math.log(f3_dict[key] / f2_dict[key[1:3]])
                    code.append(element_1 + element_2 - element_3)
                encodings.append(code)
            self.encoding_array = np.array([])
            self.encoding_array = np.array(encodings, dtype=str)
            self.column = self.encoding_array.shape[1]
            self.row = self.encoding_array.shape[0] - 1
            del encodings
            if self.encoding_array.shape[0] > 1:
                return True
            else:
                return False
        except Exception as e:
            self.error_msg = str(e)
            return False


    def get_data(self):
        if self.encoding_array.shape[0] >= 2:
            return self.encoding_array[1:, 1:]
        else:
            return None

        
if __name__ == '__main__':
    M = {
        'sliding_window': 2
    }

    seq=Descriptor(input("enter your srna file "),M)
    srna_identifier=identifier
    ssrna_identifier=srna_identifier
    pro=Descriptor(input("enter your protein file "),M)
    pro_identifier=identifier
    
    
    #get all sRNA features
    seq.NAC()
    NAC=seq.get_data()
    seq.DNC()
    DNC=seq.get_data() 
    seq.TNC()
    TNC=seq.get_data()
    seq.ASDC()
    ASDC=seq.get_data()
    seq.MMI()
    MMI=seq.get_data()
    seq.CKSNAP()
    CKSNAP=seq.get_data()

    
    #get all protein features
    pro.Protein_AAC()
    Protein_AAC=pro.get_data()
    pro.Protein_DPC()
    Protein_DPC=pro.get_data()
    pro.Protein_DDE()
    Protein_DDE=pro.get_data()
    pro.Protein_GAAC()
    Protein_GAAC=pro.get_data()
    pro.Protein_GDPC()
    Protein_GDPC=pro.get_data()
    pro.Protein_GTPC()
    Protein_GTPC=pro.get_data()
    pro.Protein_CTDC()
    Protein_CTDC=pro.get_data()
    pro.Protein_CTDT()
    Protein_CTDT=pro.get_data()
    pro.Protein_CTDD()
    Protein_CTDD=pro.get_data()
    pro.Protein_ASDC()
    Protein_ASDC=pro.get_data()

    
    #concatenate all sRNA feature sets and then select 15 needed features
    sRNA = np.concatenate((NAC,DNC,TNC,ASDC,MMI,CKSNAP),axis=1)
    DF_srna = pd.DataFrame(sRNA)
    col_list = [2,12,15,21,48,91,92,109,110,111,112,114,158,167,173]
    DF_srna = DF_srna[col_list]
    
    #concatenate all protein feature sets and then select 454 needed features
    protein=np.concatenate ((Protein_AAC,Protein_DPC,Protein_DDE,Protein_GAAC,Protein_GDPC,Protein_GTPC,Protein_CTDC,Protein_CTDT,Protein_CTDD,Protein_ASDC),axis=1)
    DF_protein=pd.DataFrame(protein)
    col_list1=[2,4,5,6,7,9,10,12,13,14,16,18,19,24,26,27,30,31,32,33,34,35,36,37,38,57,59,63,64,65,68,75,76,79,80,81,82,83,85,86,87,88,89,97,102,104,109,111,113,114,121,127,128,131,133,134,135,136,139,141,142,146,147,148,149,151,154,155,156,157,159,161,166,171,174,177,181,182,183,184,191,196,198,200,201,202,206,207,209,212,214,216,220,222,224,225,226,227,231,232,233,235,236,237,239,243,244,245,246,249,270,426,427,428,429,430,431,432,433,434,435,436,437,440,441,444,445,447,449,450,451,452,455,456,457,458,459,460,461,462,463,465,466,467,469,471,472,473,474,475,477,479,480,481,482,483,484,485,487,824,828,830,831,832,833,834,835,836,837,838,839,842,843,844,845,856,861,863,864,868,869,878,879,881,886,887,888,890,891,893,895,896,897,898,900,904,906,909,910,911,913,914,915,916,917,918,920,923,924,926,927,985,986,987,988,989,990,991,992,994,995,996,997,1024,1025,1027,1029,1030,1031,1032,1033,1034,1035,1036,1037,1038,1039,1040,1041,1042,1043,1044,1062,1063,1066,1067,1068,1069,1070,1073,1077,1078,1079,1080,1084,1085,1086,1087,1088,1092,1093,1094,1095,1097,1098,1100,1101,1102,1103,1104,1108,1110,1111,1113,1114,1115,1117,1118,1119,1120,1121,1122,1124,1125,1126,1127,1128,1129,1130,1131,1132,1133,1135,1136,1137,1138,1139,1140,1141,1144,1145,1146,1148,1149,1150,1151,1152,1153,1154,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1272,1273,1276,1277,1279,1280,1281,1282,1283,1285,1286,1287,1289,1291,1292,1293,1295,1295,1296,1298,1299,1300,1302,1303,1304,1305,1306,1307,1308,1311,1312,1314,1315,1316,1317,1318,1319,1320,1322,1323,1324,1327,1329,1330,1331,1332,1333,1336,1337,1338,1339,1340,1341,1342,1344,1345,1346,1348,1349,1350,1351,1352,1353,1355,1356,1357,1358,1359,1360,1361,1363,1364,1365,1367,1368,1369,1370,1371,1373,1374,1375,1376,1377,1378,1379,1380,1381,1383,1385,1386,1387,1388,1390,1391,1392,1393,1294,1395,1396,1399,1400,1401,1402,1403,1404,1405,1407,1408,1411,1412,1413,1414,1417,1418,1419,1421,1423,1424,1425,1427,1428,1429,1430,1431,1432,1434,1435,1437,1439,1440,1441,1442]
    DF_protein = DF_protein[col_list1]

    #merge two dataframes into one
    DF_srna['key'] = 1
    DF_protein['key'] = 1
    DF_result= pd.merge(DF_srna, DF_protein, on='key')
    DF_result= DF_result.drop('key', axis=1)
    new_cols = ['1', '10', '13',"19","45","85","86","87","99","100","102","104","140","149","155","182","183","185","186","187","188","190","191","193","194","195","197","199","200","202","204","205","208","209","210","211","212","213","214","215","216","235","237","241","242","243","246","253","254","257","258","259","260","261","263","264","265","266","267","275","280","282","287","289","291","292","299","305","306","309","311","312","313","314","317","319","322","324","325","326","327","329","332","333","334","335","337","339","344","349","352","355","359","360","361","362","369","374","376","378","379","380","384","385","387","390","391","394","398","400","402","403","404","407","409","410","411","413","414","415","417","421","422","423","425","427","448","457","458","459","460","461","462","463","464","465","466","467","468","471","472","475","476","478","480","481","482","483","486","487","488","489","490","491","492","493","494","496","497","498","500","502","503","504","505","506","508","510","511","512","513","514","515","516","518","519","521","522","523","524","525","526","527","528","529","530","531","534","535","536","537","538","543","545","546","550","551","560","561","563","568","569","570","572","573","575","577","578","579","580","582","586","588","591","592","593","595","596","597","598","599","600","602","605","606","608","609","612","613","614","615","616","617","618","619","621","622","623","624","626","627","629","631","632","633","634","635","636","637","638","639","640","641","642","643","644","645","646","647","648","651","652","653","654","655","658","662","663","664","665","669","670","671","672","673","677","678","679","680","682","683","685","686","687","688","689","693","695","696","698","699","700","702","703","704","705","706","707","709","710","711","712","713","714","715","716","717","718","720","721","722","723","724","725","726","729","730","731","733","734","735","736","737","738","739","741","742","743","744","745","746","747","748","749","750","751","752","753","754","755","758","759","761","762","763","764","765","767","768","769","771","773","774","775","776","777","778","780","781","782","784","785","786","787","788","789","790","793","794","795","796","797","798","799","800","801","802","804","805","806","809","811","812","813","815","818","819","820","821","822","823","824","826","827","828","830","831","832","833","834","835","837","838","839","840","841","842","843","845","846","847","849","850","851","852","853","855","856","857","858","859","860","861","862","863","865","867","868","869","870","872","873","874","875","876","877","878","881","882","883","884","885","886","887","889","890","893","894","895","896","899","900","901","903","905","906","907","909","910","911","912","913","914","916","917","919","921","923","924"]
    DF_result.columns = new_cols
    
    #get number of sequences in protein and sRNA fasta files
    srna_rows = len(DF_srna)
    srna_rows=srna_identifier.shape[0]
    protein_rows = pro_identifier.shape[0]
    sequences=srna_rows*protein_rows
    srna_repeat=sequences / srna_rows
    protein_repeat=sequences / protein_rows
    srna_result = []
    protein_result = []
    for i in range(int(srna_repeat)):
        for j in range(1, srna_rows+1):
            srna_result.append(srna_identifier.iloc[j-1])
   
    for l in range(1, protein_rows+1):
        for k in range(1, int(protein_repeat+1)):
            protein_result.append(pro_identifier.iloc[l-1])
      
    #using trained model for prediction
    loaded_model = joblib.load("Completed_model.joblib")
    DF_result = DF_result.astype({col: 'float64' for col in DF_result.columns})
    ynew=loaded_model.predict_proba(DF_result)
    ynew_df = pd.DataFrame(ynew)
    df_protein_result=pd.DataFrame(protein_result)
    df_srna_result=pd.DataFrame(srna_result)
    df_srna_result = df_srna_result.reset_index(drop=True)
    df_protein_result = df_protein_result.reset_index(drop=True)
    ynew_df = ynew_df.reset_index(drop=True)

    for i in range(len(ynew_df)):
        result = pd.concat([df_srna_result, df_protein_result,ynew_df.iloc[:, 1],ynew_df.iloc[:, 0]],axis=1)
        result.columns=["sRNA ID","protein ID","Propability of interaction","Propability of non-interaction"]
        result_sorted = result.sort_values(by="Propability of interaction",ascending=False)
    result_sorted.to_csv('RPI369.csv', index=False)
    print(result_sorted) 
   
  
    


# In[7]:


pip list


# In[ ]:




