#!/usr/bin/env python3
# Identify populations
# Iterate over VCF data lines
# 	Drop anything not bi-allelic
# 	Count 0 and 1 allele calls
# 	Output
# 		Counts of each allele per population
# 			[Ref]	Pop1	Pop2	Pop3...
# 			[Alt]	Pop1	Pop2	Pop3...
# 		Separate file listing accepted loci
# 			Chrom	Pos
# 
# 
# Pseudocode
# 
# vcf = open('file','r')
import re  #will use a regular expression later, so we import the appropriate library now
inputFile = '/Volumes/SlowStore/Collaboratory/pythonUsers/final.recode.vcf'  #made the input file path a variable to make naming output files easier
vcf = open(inputFile,'r')  #opening our input file.  In a later example, we can validate this as an input by ensuring the file actually exists.  We could also allow for it to be passed as a commandline parameter for use in a bash script.
#next two lines open our output files.  See below (OUTPUTS section) to find out why we open them here before entering the loop
locusList = open(inputFile + '.loci','w')  #open an output file for listing accepted loci.  Name the file with the same name as the original input plus .loci .  Lines here will match with lines in the frequency matrix
frequencyMatrix = open(inputFile + '.counts','w') #open an output file for our frequency matrix, use a similar naming scheme to the previous file

# skip over lines starting with '##'
populationsList = [] #initialize an empty list to hold our populations (this will just list any possible populations)
columnPopulations = [] #this will be an ordered list of which population a given column (matched by the order) belongs to
line = vcf.readline()
while line:
    if line[0:2] == '##':
        line = vcf.readline()
        continue
    
# read line that starts with '#' to get our population IDs

    if line[0] == '#':
        line = line.split("\t")

# split line on tabs
        
        
        for field in line[9:len(line)]:  #iterate over columns 9 through the end
        
# ignore constant VCF columns
# extract group ID from column headers

            population = re.search(r'^(\D+)',field)
            if not population: #if our regex returned nothing
                population = False #declare population to be false so we can handle it later
            if population: #but if it did find something
                population = population.group(0) #declare it to tbe the value the regex returned

# 	Regex for ^\D+

            if population:  #this will only proceed if the column returns a population (validating our inputs)
                if population not in populationsList:  #if the population we are looking at is not in the list of possible populations, we add it to the list
                    populationsList.append(population)
                columnPopulations.append(population) #add the population to the ordered list containing the data on which population belongs to which column
            else: #this will deal with any column where the population returned as null (due to a bad header)
                columnPopulations.append(False) #we set the value to a simple false so that we can later quickly identify the column as bogus
        populationsList.sort() #we will sort the populations list (just to have some consistent order)
        line = vcf.readline() #read in the next line of the file (this should be the first dataline)
    
# 	Add returned regex value to a list identifying which group the column belongs to

        continue #and return to the beginning of the loop, which will start us iterating over the datalines

# iterate over datalines
        
    else:
        line = line.split("\t") #split the line on tabs
       
# 	split the line on tabs

        if "," in line[4]: #if the line has a comma in the alternate allele (because more than one alternate allele was identified)
            line = vcf.readline()  #read the next line into our line variable
            continue  #and move on to the next iteration of this loop, skipping any other actions with this line

# 	drop lines that are not bi-allelic
# 		if in column4 there is a ',' skip the line
# 	   	else
        
        freqData = {} #initialize a blank dictionary or ensure that this dictionary is blank (meaning no possible hold-overs from previous lines to cause errors)
        for population in populationsList: #iterate over the list of populations
            freqData[population] = {} #initialize an empty dictionary at the population level.  Unlike PERL or some other languages, you can't directly initialize multiple levels at once.
            for genotype in ['0','1']: #iterate over possible genotype values
                freqData[population][genotype] = 0 #initialize a dictionary of dictionaries indexed by population and genotype with all zero values
        
# 	initialize a dictionary of dictionaries 
# 		data{population}{'0'} = 0
# 		data{population}{'1'} = 0
# 		for each population

        dataColumns = line[9:len(line)] #cheating a bit and creating a new dataset containing only the columns with individual sample data values.  This would be a problem if the lines were especially long, as it would almost double our memory needs, but they are unlikely to be that long ever, and now the indexing of our new set matches that of our list of columnPopulations
        if len(dataColumns) != len(columnPopulations):  #INPUT VALIDATION... if the number of columns we found populations for is not the same as the number of columns we have data for, something is seriously wrong (we lost or gained a column or two)
            raise RuntimeError("Error on variant at " + line[0] + ":" + line[1] + " incorrect number of columns observed.") #We will raise an exception because this is likely a sign that the VCF has been somehow corrupted and needs to be carefully examined
        
        for columnNumber in range(0,len(dataColumns)):  #iterate over our data columns (in this case, maintaining columnNumber as our index counter so we can compare across lists)
            if dataColumns[columnNumber][0] not in ['0','1'] and dataColumns[columnNumber][2] not in ['0','1']: #INPUT VALIDATION... check if both genotypes are valid (only 0 or 1 should be valid).  Uncalled genotypes will be read as a '.' and anything else will be trapped here as well
                continue #and if no genotype was called for the sample, iterate on to the next column
            currentPopulation = columnPopulations[columnNumber]  #this is why we broke off the data columns, we have matching index values and don't have to worry about subtracting
            if not currentPopulation: #make sure that we have a current population value.  We generally should, but if a column had a bad header that didn't contain a valid population identifier, this will catch it
                continue #and we iterate on to the next column without doing anything further
            genotype = dataColumns[columnNumber][0]  #we grab the first genotype from the data
            freqData[currentPopulation][genotype] += 1 #and use it and the current population value to increment the appropriate count by 1 (the += operation says to add 1 to the current value of the variable) 
            genotype = dataColumns[columnNumber][2]  #we grab the second genotype from the data
            freqData[currentPopulation][genotype] += 1 #and use it and the current population value to increment the appropriate count by 1 (notice how repetative this is, maybe in future versions this could become a function somehow)
# this could also be written as    freqData[columnPopulations[columnNumber]][dataColumns[columnNumber][0]] += 1
# and                              freqData[columnPopulations[columnNumber]][dataColumns[columnNumber][2]] += 1   but that would be very confusing.  Try to puzzle through it if you want to better understand indexing.
# 	for column 9 through last column
# 		get column[0]
# 		if column[0] is a '.' (indicates missing data
# 			go to the next column
# 			else
# 		data{population}{column[0]} += 1
# 		do the same for column[2]
# 	OUTPUTS
#at this point, you should realize that we don't want to be constantly opening our output files, which will happen if we put them inside the loop.  We can either open them within an if statement, or open them before starting to loop.  For simplicity, we will put them before the loop.

        delimiter = "\t"  #next group of lines are just saving values from the dataline to variables.  This is slightly less efficient than calling them directly, but simpler to read at this level
        contig = line[0]
        position = line[1]
        refAllele = line[3]
        altAllele = line[4]
        for allele in [refAllele, altAllele]: #iterate over the genotypes (in the order listed)
            locusList.write(contig + delimiter + position + delimiter + allele + "\n") #and write the locus info line for each one
        
        for allele in ['0','1']:
            for population in populationsList:
                currentCount = freqData[population][allele]  #note that this value is an integer, not a string at this point, we risk errors if we try to write it directly
                currentCount = str(currentCount) #so we force it to become a string
                frequencyMatrix.write(currentCount + delimiter)  #and we write it plus our chosen delimiter (probably a tab) to the output file (but do not yet create a new line, as we might have more values to write)
            frequencyMatrix.write("\n") #once we reach the end of our populations, we terminate the line and write data for the next allele if we didn't just do so
        
        line = vcf.readline() #and we get a new line of data from the input file before starting the loop again (note that if we reach the end of the file and have no data, the loop will exit)  
# 		accepted loci file
# 						Chr	Pos	Ref
# 						Chr	Pos	Alt
# 			write column[0]\t column[1]\t column[3]
# 			write column[0]\t column[1]\t column[4]
# 	
# 		frequency data file
# 			iterate over alleles (0 and 1)
# 				iterate over populations IN ORDER
# 					write data{population}{allele}\t
# 				write \n
        
locusList.write('#Population order:' + delimiter) #just to be sure that the order of populations in the matrix is explicit, we add a last line to the locus list (we add it after the end so that the locus lines will still match to the data)
for population in populationsList:  #iterate over the list of populations in the order we were sorted them
    locusList.write(population + delimiter) #and write them to the file, separated by our delimiter of choice (generally a tab)
            
vcf.close()  #and we tidy up by closing our files
locusList.close()
frequencyMatrix.close()
print('Done!') #and saying we have finished before we
quit()
        