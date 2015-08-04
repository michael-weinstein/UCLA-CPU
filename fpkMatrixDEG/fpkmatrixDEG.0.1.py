#!/usr/bin/env python3


def getGenesOfInterest(geneListFile):
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
        
def createFpkmDict(deFile, geneList):
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
        if currentLine.tracking_id not in geneList:    #if the gene for this line is not one of the genes of interest
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
    return (interestDict, orderList)

def matrixOutput(filename, fpkmDict, fpkmOrder):
    output = open(filename, 'w')  #open the file we plan to write
    i = 0
    while i < len(fpkmOrder):  #go over the fpkmOrder list (each line being a list of gene, condition in the order of the original file)
        gene = fpkmOrder[i][0]
        conditions = []
        conditions.append(fpkmOrder[i][1])
        lookAhead = i + 1
        try:
            while fpkmOrder[lookAhead][0] == gene:
                conditions.append(fpkmOrder[lookAhead][1])
                lookAhead += 1
        except IndexError:
            i = len(fpkmOrder) + 1
        i = lookAhead
        outputList = []
        for condition in conditions:
            replicates = list(fpkmDict[gene][condition].keys())
            replicates.sort()
            for replicate in replicates:
                outputList.append(str(fpkmDict[gene][condition][replicate]))
        outputLine = gene + "\t" + "\t".join(outputList) + "\n"
        output.write(outputLine)
    output.close()
    return True      
    
def main():
    geneList = getGenesOfInterest('SIG2_GENE_LIST')
    deFile = open('genes.read_group_tracking', 'r')
    fpkmData = createFpkmDict(deFile, geneList)
    fpkmDict = fpkmData[0]
    fpkmOrder = fpkmData[1]
    deFile.close()
    matrixOutput('DEGmatrix', fpkmDict, fpkmOrder)
    print('Done!')
    
main()
                
        
        
    
    
    