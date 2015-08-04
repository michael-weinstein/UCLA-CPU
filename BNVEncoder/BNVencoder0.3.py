#!/usr/bin/env python3

import re  #import the regex library
import os  #import the OS calling library

class CheckArgs(object):
    def __init__(self):
        import argparse #loads the required library for reading the commandline
        import os  #imports the library we will need to check if the file exists
        parser = argparse.ArgumentParser()
        parser.add_argument ("-f", "--VCFinput", help = "VCF file for collecting data")
        parser.add_argument ("-9", "--clobber", help = "Ignore potential file overwrites (use with caution).", action = "store_true")
        parser.add_argument ("-u", "--useUnidentifiableGroup", help = "Allows the program to use a VCF even if some group IDs can't be identified.  They will be grouped into the group 'Unidentifiable' in the output.", action = "store_true")
        args = parser.parse_args()  #puts the arguments into the args object
        self.VCF = args.VCFinput
        self.useUnidentifiableGroup = args.useUnidentifiableGroup
        self.delimiter = "\t"  #Putting this here for now in case we ever need to use a different delimiter
        if not self.VCF:
            quit('No input VCF specified.')
        if not os.path.isfile(self.VCF):
            quit('Unable to find input VCF: ' + self.VCF)
        self.outputMatrix = self.VCF + ".counts"
        self.outputLoci = self.VCF + ".loci"
        if os.path.isfile(self.outputMatrix) or os.path.isfile(self.outputLoci):  #checks to see if the user set to clobber existing files automatically
            if args.clobber:  #if so, let them know we are overwriting existing files by their command
                print('Outputs already exist.  Set to overwrite in command line arguments.')
            else:  #if not, check if they want to proceed anyway
                if not yesanswer('Output files already exist for this input.  Overwrite them?'):
                    quit('OK!')


     
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
    
    def __init__(self, data): #takes in a list to initialize
        self.data = [str(datum) for datum in data] #makes values in self.data strings
        self.generate() #run the generate function as part of initializing
        
    def generate(self): #generate some useful properties
        self.hashTable = {} #initialize an empty dictionary
        for datum in self.data:  #iterate over data set
            try: #check if the key exists
                self.hashTable[datum] += 1  #try incrementing the correspoding counter
            except KeyError:  #if it didn't exist
                self.hashTable[datum] = 1  #initialize it to a zero
        keyList = list(self.hashTable.keys())  #turn keys into a list 
        keyList.sort() #sort the list
        self.map = keyList  #make the sorted list into an attribute
        self.reduced = [] #initialize an empty list
        for key in self.map:  #iterate over the mapped list
            self.reduced.append([key,self.hashTable[key]]) #use the list of keys and the hashTable to create a two column table of entry and value
        return True


    

class VCFLine(object):
    
    def __init__(self, data, delimiter = "\t"): #take an argument of data and optional specified delimiter
        self.line = data.split(delimiter) 
        self.variantColumns = self.line[0:9] #Takes the first set of columns and turns it into a property
        self.sampleColumns = self.line[9:len(self.line)] #does the same for the second set of columns with sample data

class Header(VCFLine):  #defines Header as an extension class of VCFLine
    
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
            self.columnGroupIDs.append(population)  #add the population to an ordered list of which population the column belongs to
        mapreduce = MapReduce(self.columnGroupIDs)  #call the MapReduce class
        self.groupCountTable = mapreduce.reduced  #initialize an attribute with the number of samples per group
        self.outputGroupColumns = mapreduce.map  #make an attribute with an ordered list of groups (one entry per group)
        self.groupHash = {} #initialize an empty dictionary
        for item in self.outputGroupColumns:  #iterate over our list of groups
            self.groupHash[item] = [0,0]  #for each group as a key in a dictionary, initialize an empty list of 0,0.  Conveniently enough, 0 is the genotype number for the ref allele being called and the index in our list for the reference allele count.  Alt alleles are indicated by 1, which is also the index of their counter in the list.
        return True
    
class Data(VCFLine):
    
    def integrityCheck(self, headerColumns):  #simple function to make sure that we have at least as many columns of data as we have headers
        if len(self.sampleColumns) == len(headerColumns):
            return True
        else:
            return False
    
    def isBiallelic(self):  #this function checks for biallelic variants by looking for a comma in the alt allele column (if present, indicates more than one allele)
        self.locus = self.variantColumns[0] + ":" + self.variantColumns[1]  #also captures the locus as an attribute... this is for the purpose of reporting the issue back to the user
        if "," in self.variantColumns[4]:  #does the actual check
            return False
        else:
            return True
        
    def createOutputs(self, columnGroupIDs, groupHash, headerColumns, delimiter = "\t"):  #supervisor function that tells the object to create its outputs
        self.createAlleleCountOutputs(columnGroupIDs, groupHash, headerColumns, delimiter)
        self.createLocusInfoOutputs(delimiter)
        
    def createAlleleCountOutputs(self, columnGroupIDs, groupHash, headerColumns, delimiter):  #function that creates the output lines for the ref and alt alleles
        position = 0  #initializes an integer to 0 to mark our position within the data columns
        for item in self.sampleColumns:  #iterates over the columns containing sample data
            genotypes = [item[0], item[2]]  #captures the first and third character from the column.  These should always correspond to the first and second genotype called for the locus
            if "." in genotypes:  #if the genotype contains a period, indicating that it was not called
                position += 1  #jump to the next column
                continue  #and then restart the loop
            currentGroup = columnGroupIDs[position] #if we reach this point, we have a called genotype.  Using our position tracker, we will look for the corresponding group ID to know which group to add the current data to
            for genotype in genotypes: #we iterate over the genotypes for this sample
                groupHash[currentGroup][int(genotype)] += 1  #and increment the appropriate genotype (0 for ref, 1 for alt) in our collection of genotypes
            position += 1  #and then increment our position counter
        self.refCountsOutput = ""  #initializing empty strings for ref and alt count output lines
        self.altCountsOutput = ""
        for column in headerColumns: #iterate over our header columns (headers for the output file, that is)
            self.refCountsOutput += str(groupHash[column][0]) + delimiter #and we build our output line by adding on the appropriate count from our hash of values plus a delimiter
            self.altCountsOutput += str(groupHash[column][1]) + delimiter #do the same for alternate reads (all that has to change is looking at position 1 instead of position 0) 
        return True
    
    def createLocusInfoOutputs(self, delimiter):  #creates outputs for our file containing the locus information.  This just defines several attributes based off of their column position in the VCF
        self.contig = self.variantColumns[0]  
        self.position = self.variantColumns[1]
        self.refAllele = self.variantColumns[3]
        self.altAllele = self.variantColumns[4]
        self.locusInfoRef = self.contig + delimiter + self.position + delimiter + self.refAllele  #this line and the next combine all the pertinent information about the ref and alt allele loci and have it ready to output to a file
        self.locusInfoAlt = self.contig + delimiter + self.position + delimiter + self.altAllele
        return True  #mostly useless return, but can be useful if we need some indication that this has actually been run and also provides a clear marker for the end of this function
            


def main():
    args = CheckArgs()  #create an object holding our VALIDATED commandline arguments  (bogus arguments would have made this program quit)
    vcf = open(args.VCF,'r')  #open the VCF for reading
    delimiter = args.delimiter
    locusListFileName = args.outputLoci 
    frequencyMatrixFileName = args.outputMatrix  
    locusList = open(locusListFileName,'w')   #this and the next line open our two output files
    frequencyMatrix = open(frequencyMatrixFileName,'w')
    line = vcf.readline()  #reads a line from the VCF
    counter = 0 #initializes our counter variable to indicate progress
    while line:  #iterates so long as we have an input line.  We must remember to update the line at the end of each iteration, otherwise we keep working over the same line and have an infinite loop
        print("Processed " + str(counter) + " lines.", end="\r")  #display the updated counter to the user.  Does not allow the print to start a new line (so the next update will overwrite the current counter)
        counter += 1  #increment the counter
        if line[0] == '#':  #checks if the line starts with a #
            if line[1] == '#':  #if it does and the second character is also #
                line = vcf.readline()  #read in a new line
                continue #and start the loop over
            else:
                header = Header(line)  #if the line starts with only a single #, it must be the headers, so we initialize our header object (this will be available for reference throughout the run of this program)
                header.generateLists(args.useUnidentifiableGroup)  #tells the header to make its lists.  Just needs to know (in the arguments passed) if it has permission to use groups that it cannot identify
                groupColumns = delimiter.join(header.outputGroupColumns)  #creates a delimited string of our groups for the locus file
                locusList.write("contig" + delimiter + "position" + delimiter + "allele" + delimiter + groupColumns + "\n") #writes the column names (including group IDs) to the locus file
                line = vcf.readline()  #reads next line
                continue  #starts the loop again
        else:  #if the line did not start with a #, it must be a data line
            data = Data(line)  #initialize an object to handle the data line
            if not data.integrityCheck(header.columnGroupIDs):  #if the line fails integrity check (wrong number of columns, probably due to a corruption of the file)
                print("Warning: Incorrect number of columns found for locus " + data.locus + ".  Skipping this locus.")  #warn the user
                line = vcf.readline()  #read the next line
                continue  #and continue to the next iteration of the loop
            if not data.isBiallelic():  #if the line is not for a biallelic locus
                print("Warning: Multiple alternative alleles found for locus " + data.locus + ".  Skipping this locus.")  #warn the user
                line = vcf.readline()  #read the next line
                continue  #and continue to the next iteration of the loop
            data.createOutputs(header.columnGroupIDs, header.groupHash, header.outputGroupColumns, delimiter)  #run the supervisor function for a data line that makes the outputs
            frequencyMatrix.write(data.refCountsOutput + "\n")  #write the newly-created reference output to the matrix file
            frequencyMatrix.write(data.altCountsOutput + "\n")  #do the same on the next line for the alt output
            locusList.write(data.locusInfoRef + "\n")  #then write the locus data to the appropriate file in the same order on this and the next line
            locusList.write(data.locusInfoAlt + "\n")
            line = vcf.readline() #and read the next line of the input file before starting the loop again
    print("Processed " + str(counter) + " lines.")
    vcf.close() #close all the files we were working on
    locusList.close() 
    frequencyMatrix.close()
    quit("Done!")
    
main()
            
            