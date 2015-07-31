#!/usr/bin/env python3

import re
import os

class Arguments(object):
    
    def __init__(self):
        self.inputFileName = "karinasample.vcf"
        self.clobber = True
        self.useUnidentifiableGroup = True


     
def yesanswer(question):  #asks the question passed in and returns True if the answer is yes, False if the answer is no, and keeps the user in a loop until one of those is given.  Also useful for walking students through basic logical python functions
    answer = False  #initializes the answer variable to false.  Not absolutely necessary, since it should be undefined at this point and test to false, but explicit is always better than implicit
    while not answer:  #enters the loop and stays in it until answer is equal to True
        print (question + ' (Y/N)')  #Asks the question contained in the argument passed into this subroutine
        answer = input('>>') #sets answer equal to some value input by the user
        if str(answer) == 'y' or str(answer) == 'Y':  #checks if the answer is a valid yes answer
            return True  #sends back a value of True because of the yes answer
        elif str(answer) == 'n' or str(answer) == 'N': #checks to see if the answer is a valid form of no
            return False  #sends back a value of False because it was not a yes answer
        else: #if the answer is not a value indicating a yes or no
            print ('Invalid response.')
            answer = False #set ansewr to false so the loop will continue until a satisfactory answer is given
            



class MapReduce(object):
    
    def __init__(self, data):
        self.data = [str(datum) for datum in data]
        self.generate()
        
    def generate(self):
        self.hashTable = {}
        for datum in self.data:
            try:
                self.hashTable[datum] += 1
            except KeyError:
                self.hashTable[datum] = 0
        keyList = list(self.hashTable.keys())
        keyList.sort()
        self.map = keyList
        self.reduced = []
        for key in self.map:
            self.reduced.append([key,self.hashTable[key]])
        return True


    

class VCFLine(object):
    
    def __init__(self, data, delimiter = "\t"):
        self.line = data.split(delimiter)
        self.variantColumns = self.line[0:9]
        self.sampleColumns = self.line[9:len(self.line)]

class Header(VCFLine):
    
    def generateLists(self, useUnidentifiable):
        self.columnGroupIDs = []
        for item in self.sampleColumns:
            population = ""
            regex = re.search(r'^(\D+)',item)
            if not regex:
                if useUnidentifiable:
                    population = "Unidentifiable"
                else:
                    raise RuntimeError('Unable to extract population from ' + item + '.  Not set to use unidentifiable populations.')
            if not population == "Unidentifiable":
                population = regex.group(0)
            self.columnGroupIDs.append(population)
        mapreduce = MapReduce(self.columnGroupIDs)
        self.groupCountTable = mapreduce.reduced
        self.outputGroupColumns = mapreduce.map
        self.groupHash = {}
        for item in self.outputGroupColumns:
            self.groupHash[item] = [0,0]
        return True
    
class Data(VCFLine):
    
    def integrityCheck(self, headerColumns):
        if len(self.sampleColumns) == len(headerColumns):
            return True
        else:
            return False
    
    def isBiallelic(self):
        self.locus = self.variantColumns[0] + ":" + self.variantColumns[1]
        if "," in self.variantColumns[4]:
            return False
        else:
            return True
        
    def createOutputs(self, columnGroupIDs, groupHash, headerColumns, delimiter = "\t"):
        self.createAlleleCountOutputs(columnGroupIDs, groupHash, headerColumns, delimiter)
        self.createLocusInfoOutputs(delimiter)
        
    def createAlleleCountOutputs(self, columnGroupIDs, groupHash, headerColumns, delimiter):
        position = 0
        for item in self.sampleColumns:
            genotypes = [item[0], item[2]]
            if "." in genotypes:
                continue
            currentGroup = columnGroupIDs[position]
            for genotype in genotypes:
                groupHash[currentGroup][int(genotype)] += 1
            position += 1
        self.refCountsOutput = ""
        self.altCountsOutput = ""
        for column in headerColumns:
            self.refCountsOutput += str(groupHash[column][0]) + delimiter
            self.altCountsOutput += str(groupHash[column][1]) + delimiter
        return True
    
    def createLocusInfoOutputs(self, delimiter):
        self.contig = self.variantColumns[0]
        self.position = self.variantColumns[1]
        self.refAllele = self.variantColumns[3]
        self.altAllele = self.variantColumns[4]
        self.locusInfoRef = self.contig + delimiter + self.position + delimiter + self.refAllele
        self.locusInfoAlt = self.contig + delimiter + self.position + delimiter + self.altAllele
        return True
            


def main(delimiter = "\t"):
    args = Arguments()
    if not os.path.isfile(args.inputFileName):
        raise FileNotFoundError('Unable to find input file: ' + args.inputFileName)
    vcf = open(args.inputFileName,'r')
    locusListFileName = args.inputFileName + '.loci'
    frequencyMatrixFileName = args.inputFileName + '.counts'
    if os.path.isfile(locusListFileName) or os.path.isfile(frequencyMatrixFileName):
        if args.clobber:
            print('Overwriting existing outputs.')
        elif not yesanswer('Output file (s) for ' + args.inputFileName + ' already exist.\nDo you wish to overwrite?'):
            quit('OK. Exiting.')
    locusList = open(locusListFileName,'w')
    frequencyMatrix = open(frequencyMatrixFileName,'w')
    line = vcf.readline()
    while line:
        if line[0] == '#':
            if line [1] == '#':
                line = vcf.readline()
                continue
            else:
                header = Header(line)
                header.generateLists(args.useUnidentifiableGroup)
                groupColumns = delimiter.join(header.outputGroupColumns)
                locusList.write("contig" + delimiter + "position" + delimiter + "allele" + delimiter + groupColumns + "\n")
                line = vcf.readline()
                continue
        else:
            data = Data(line)
            if not data.integrityCheck(header.columnGroupIDs):
                print("Warning: Incorrect number of columns found for locus " + data.locus + ".  Skipping this locus.")
                line = vcf.readline()
                continue
            if not data.isBiallelic():
                print("Warning: Multiple alternative alleles found for locus " + data.locus + ".  Skipping this locus.")
                line = vcf.readline()
                continue
            data.createOutputs(header.columnGroupIDs, header.groupHash, header.outputGroupColumns, delimiter)
            frequencyMatrix.write(data.refCountsOutput + "\n")
            frequencyMatrix.write(data.altCountsOutput + "\n")
            locusList.write(data.locusInfoRef + "\n")
            locusList.write(data.locusInfoAlt + "\n")
            line = vcf.readline()
    vcf.close()
    locusList.close()
    frequencyMatrix.close()
    quit("Done!")
    
main()
            
            