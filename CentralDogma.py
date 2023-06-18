from sys import *

#read the in the nucleotide sequence as a file object
inputFile = open(argv[1],'r')

#read in conversion table
conversionTable = open(argv[2],'r')

#getCodons function: splits a nucleotide string into a list of codons
def getCodons(nucleotideString):
    #create a string of nucleotides with the codons separated by spaces
    spacedNucleotides = ''

    #add a space between every 3 nucleotides
    for i in range(0,len(nucleotideString),3): 
        spacedNucleotides += nucleotideString[i:i+3]
        spacedNucleotides += ' '

    #create a list of codons
    codonList = spacedNucleotides.strip().split()
    return codonList

#getCodonLists function: returns 3 lists of codons for the 3 reading frames of a nucleotide string input 
def getCodonLists(nucleotideString):
    #define codons for reading frame 1
    frame1Codons = getCodons(nucleotideString)

    #define codons for reading frames 2 & 3
    nucleotideString2 = nucleotideString[1:]
    nucleotideString3 = nucleotideString[2:]

    frame2Codons = getCodons(nucleotideString2)
    frame3Codons = getCodons(nucleotideString3)

    #delete the last element of the list with < 3 nucleotides
    ele = frame2Codons.pop()
    ele = frame3Codons.pop()

    listOfLists = [frame1Codons, frame2Codons, frame3Codons]
    return listOfLists

#getStartIndexes: takes in a list of codons & returns the indexes of start codons
def getStartIndexes(readingFrame):
    indexes = []
    for i in range(0,len(readingFrame)):
        if readingFrame[i] == 'atg':
            indexes.append(str(i))
    return indexes

#getStopIndexes: takes in a list of codons & returns the indexes of stop codons
def getStopIndexes(readingFrame):
    indexes = []
    for i in range(0,len(readingFrame)):
        if readingFrame[i] == 'taa':
            indexes.append(str(i))
        elif readingFrame[i] == 'tag':
            indexes.append(str(i))
        elif readingFrame[i] == 'tga':
            indexes.append(str(i))
    return indexes

# function to create each possible codon sequence to be translated for each reading frame 
def getTranslationCodons(readingFrame):
    translationCodons = []
    #get start & stop indexes
    startIndexes = []
    stopIndexes = []
    startIndexes = getStartIndexes(readingFrame)
    if startIndexes != []: # if there are start codon(s)
        stopIndexes = getStopIndexes(readingFrame)
        #create a codon sequence starting @ each start codon
        for i in range(0, len(startIndexes)): #i = start codon
            startCodon = int(startIndexes[i]) 
            #identify the nearest stop codon after the start codon
            if stopIndexes != []:
                for n in range(0,len(stopIndexes)):
                    if int(stopIndexes[n]) > startCodon:
                        stopCodon = int(stopIndexes[n]) 
                        break
                    elif int(stopIndexes[-1]) <= startCodon:
                        stopCodon = len(readingFrame)
                        break
            #continue til the end of the readingFrame codons if there's no stop codon
            else:
                stopCodon = len(readingFrame)
                
            translationCodons.append(readingFrame[startCodon:stopCodon+1])
    else: #if there are no start codons, return an empty list for translation codons for that reading frame
        translationCodons = []
    return translationCodons

#getProteins function: takes a list of codons as input & returns the protein abbreviation output
def getProteins(codonList):
    output = '' #create an output string
    
    for codon in codonList: #iterate through each codon
        output += codonDict[codon]
        
    return output

#function that gives all possible protein output strings for a reading frame input
def getOutput(readingFrame):
    output = ''
    translationCodons = getTranslationCodons(readingFrame)
    for n in range(0,len(translationCodons)):
        output += getProteins(translationCodons[n])+'\n'
    return output

#create codon dictionary from conversion table

lineList = [] 

for line in conversionTable:
    if line != []:
        lineList.append(line.strip().split())
#print(lineList)

cleanedList = [ele for ele in lineList if ele != []]
#print(cleanedList)


codonDict = {} #create empty dictionary

for ele in cleanedList:
    nucleotides = ele[0].lower() #make lowercase to match input nucleotides
    aminoAcid = ele[1]
    codonDict[nucleotides] = aminoAcid #define the dictionary

#create a string for the nucleotides
nucleotideString = ''

for currentLine in inputFile:
    if '>' in currentLine: 
        continue #skip the line(s) that start with '>'
    else:
        nucleotideString += currentLine.strip().lower() #add the nucleotides to the string & make sure it's all lowercase

#create nucleotide string for reverse strand
reversedString = nucleotideString[::-1] #reverse the order

#switch a/t and c/g
backwardString = ''
for i in range(0,len(reversedString)):
    if reversedString[i] == 'a':
        backwardString += 't'
    elif reversedString[i] == 't':
        backwardString += 'a'
    elif reversedString[i] == 'c':
        backwardString += 'g'
    else:
        backwardString += 'c'

#define each reading frame
forwardStrandCodonFrames = getCodonLists(nucleotideString)
forwardRF1 = forwardStrandCodonFrames[0]
forwardRF2 = forwardStrandCodonFrames[1]
forwardRF3 = forwardStrandCodonFrames[2]

backwardStrandCodonFrames = getCodonLists(backwardString)
backwardRF1 = backwardStrandCodonFrames[0]
backwardRF2 = backwardStrandCodonFrames[1]
backwardRF3 = backwardStrandCodonFrames[2]

output = 'possible proteins are:\n'

#get & print all the possible proteins for each reading frame & display them

output = 'possible proteins are:\n' + getOutput(forwardRF1) + getOutput(forwardRF2) + getOutput(forwardRF3) + getOutput(backwardRF1) + getOutput(backwardRF2) + getOutput(backwardRF3)
    
print(output)
