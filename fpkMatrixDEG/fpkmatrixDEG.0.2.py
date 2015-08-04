#!/usr/bin/env python3

'''
Python3 program for extracting fpkm data from cuffdiff outputs and writing it to a file in a two dimensional matrix of values with one gene per line
and each condition clustered.  Output will preserve the exact order of the input.  Written by Michael Weinstein, UCLA Cohn Lab and Collaboratory, 2015
'''

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

class checkArgs(object):
    def __init__(self):
        import argparse #loads the required library for reading the commandline
        import os  #imports the library we will need to check if the file exists
        parser = argparse.ArgumentParser()
        parser.add_argument ("-c", "--cuffDiffOutput", help = "CuffDiff Output File")
        parser.add_argument ("-g", "--geneList", help = "File containing a list of the genes of interest (1 per line).")
        parser.add_argument ("-9", "--clobber", help = "Ignore potential file overwrites (use with caution).", action = "store_true")
        args = parser.parse_args()  #puts the arguments into the args object
        self.geneList = args.geneList
        self.cuffDiffOutput = args.cuffDiffOutput
        self.outputMatrix = self.cuffDiffOutput + ".matrix"
        self.outputKey = self.cuffDiffOutput + ".key"
        if not args.geneList:
            print("No gene of interest list set.")
            self.geneList = False
        if not self.cuffDiffOutput:
            quit('No CuffDiff file specified.')
        if not os.path.isfile(self.cuffDiffOutput):
            quit('Unable to find CuffDiffOutput file: ' + self.cuffDiffOutput)
        if self.geneList and not os.path.isfile(self.geneList):
            quit('Unable to find gene list file: ' + self.geneList)
        if os.path.isfile(self.outputMatrix) or os.path.isfile(self.outputKey):
            if args.clobber:
                print('Outputs already exist.  Set to overwrite in command line arguments.')
            else:
                if not yesanswer('Output files already exist for this input.  Overwrite them?'):
                    quit('OK!')

def getGenesOfInterest(geneListFile):
    if not geneListFile:
        return False
    file = open(geneListFile, 'r')
    geneListLine = file.readline()
    geneList = []
    while(geneListLine):
        geneListLine = geneListLine.strip()
        geneList.append(geneListLine)
        geneListLine = file.readline()
    file.close()
    return geneList

class CuffDiffDeDataLine(object):
    def __init__(self, line):  #initializes this instance of the object, takes in the line from the cuff diff file
        self.line = line #initializes the line value to the raw input
        self.data = line.split("\t")  #initializes the data value to the split line
        if self.integrityCheck():  #checks to make sure we have enough values on the line
            self.generate()  #runs the subroutine to generate all the attributes
        else:
            raise IndexError("Line: " + self.line + " has less than 9 values.") #stops the program and reports the error
            
    def integrityCheck(self):
        if len(self.data) >= 9:  #if there are 9 or more values on the line
            return True  #returns true to indicate passing the check
        else:
            return False  #otherwise returns False
        
    def generate(self):  #takes no arguments and is called by the initializer.  Initializes all our attributes.
        self.tracking_id = self.data[0]  #initializes the tracking_id value with the first element of the line
        self.condition = self.data[1]
        self.replicate = self.data[2]
        self.raw_frags = self.data[3]
        self.internal_scaled_frags = self.data[4]
        self.external_scaled_frags = self.data[5]
        self.fpkm = self.data[6]
        self.effective_length = self.data[7]
        self.status = self.data[8]
        
def createFpkmDict(deFileName, geneList):
    deFile = open(deFileName, 'r')
    counter = 0
    line = deFile.readline()  #initialize line to the first line of the deFile
    interestDict = {} #intiialize an empty dictionary for our data
    orderList = [] #initializing an empty list for the order
    while line: #as long as there is some value in the line (so not end of file)
        print("Processed " + str(counter) + " lines.", end = "\r")
        if "tracking_id" in line:  #check if the literal string "tracking_id" is in the line
            line = deFile.readline() #if so, we have a header line, we read the next line
            counter += 1
            continue  #and go back to the start of the loop
        currentLine =  CuffDiffDeDataLine(line) #initialize the object for the line.  If we have this object left over from a previous line, we reinitalize it and get rid of the old data
        if geneList and (currentLine.tracking_id not in geneList):    #if there is a geneList and the gene for this line is not one of the genes of interest
            line = deFile.readline()  #read the next line
            counter += 1
            continue #and start the loop over again for it
        else:  #if the gene is in our list of interest
            if [currentLine.tracking_id,currentLine.condition] not in orderList:
                orderList.append([currentLine.tracking_id,currentLine.condition])
            try:  #try to directly drop the fpkm value directly into the appropriate dictionary.  If that fails due to not having a key for it, keep moving up a level and defining the key
                interestDict[currentLine.tracking_id][currentLine.condition][int(currentLine.replicate)] = float(currentLine.fpkm)
            except KeyError:
                try:
                    interestDict[currentLine.tracking_id][currentLine.condition] = {}
                    interestDict[currentLine.tracking_id][currentLine.condition][int(currentLine.replicate)] = float(currentLine.fpkm)
                except KeyError:
                    interestDict[currentLine.tracking_id] = {}
                    interestDict[currentLine.tracking_id][currentLine.condition] = {}
                    interestDict[currentLine.tracking_id][currentLine.condition][int(currentLine.replicate)] = float(currentLine.fpkm)
            counter += 1
            line = deFile.readline()  #read the next line to avoid an infinite loop
    print("Processed " + str(counter) + " lines.")
    deFile.close()
    return (interestDict, orderList)

def matrixOutput(matrixOutputFileName, keyOutputFileName, fpkmDict, fpkmOrder):
    counter = 0
    matrixOutput = open(matrixOutputFileName, 'w')  #open the file we plan to write the matrix to
    keyOutput = open(keyOutputFileName, 'w') #open the file we will write the key output to
    i = 0  #initialize our counter
    while i < len(fpkmOrder):  #go over the fpkmOrder list (each line being a list of gene, condition in the order of the original file)
        print("Wrote " + str(counter) + " lines.", end="\r")
        counter += 1
        gene = fpkmOrder[i][0]  #use the orderList to get our current gene name
        conditions = []  #initialize an empty list for the different conditions for the gene
        conditions.append(fpkmOrder[i][1])  #add the condition for the current line of the order list to our list of conditions for the current gene
        lookAhead = i + 1  #move the index ahead by one
        try:  #setting this try/except block in case we run off the end of the list (which will happen when we reach the last gene)
            while fpkmOrder[lookAhead][0] == gene:  #as long as the next line has the same gene as the gene we are working on at the moment
                conditions.append(fpkmOrder[lookAhead][1])  #add the condition to the list of conditions for the current gene
                lookAhead += 1  #and then move the index up one more, then repeat the loop
        except IndexError: #if we hit an index error (because we ran off the end of the list and finished up)
            i = len(fpkmOrder) + 1  #set the index to one more than the length of the list, ensuring we exit this loop
        i = lookAhead #after getting all the genes for the condition, our lookAhead value will be pointing at the first line for the next gene.  Set our index to that for the next iteration of the loop
        outputList = []  #initialize an empty list for collecting data to put in each line of the matrix output
        conditionCount = []  #initialize an empty list for collecting data to put in each line of the key output (indicating what each matrix line contains)
        for condition in conditions:  #iterate over the conditinos for the gene we are looking at
            replicates = list(fpkmDict[gene][condition].keys())  #get a list of keys
            replicates.sort()  #ensure that the list of replicates (should be integers) is sorted in numerical order.  This is probably not absolutely necessary, but makes it explicit that it will be ordered
            for replicate in replicates:  #iterate over our list of replicates (which we just ensured is in order)
                outputList.append(str(fpkmDict[gene][condition][replicate])) #append the string of the numerical FPKM value to our list of outputs for the matrix.  This will be ordered by condition, then by replicate
            conditionCount.append(condition + "(" + str(len(replicates)) + ")")  #append a string to our condition counts list (for our key output) that will look like condition(N) where N is the number of replicates in the matrix
        matrixOutputLine = gene + "\t" + "\t".join(outputList) + "\n"  #create the output line to the matrix file by joining the list of values with delimiters
        keyOutputLine = gene + "\t" + "\t".join(conditionCount) + "\n"  #do the same for our key output
        matrixOutput.write(matrixOutputLine) #write the appropriate output line to the matrix file
        keyOutput.write(keyOutputLine) #do the same for the key
    print("Wrote " + str(counter) + " lines.")    
    matrixOutput.close() 
    keyOutput.close()
    return True      
    
def main():
    args = checkArgs()
    geneList = getGenesOfInterest(args.geneList)
    fpkmData = createFpkmDict(args.cuffDiffOutput, geneList)
    fpkmDict = fpkmData[0]
    fpkmOrder = fpkmData[1]
    matrixOutput(args.outputMatrix, args.outputKey, fpkmDict, fpkmOrder)
    print('Done!')
    
main()