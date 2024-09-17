'''This was a group project for the class "Practical Computer Concepts for Bioinformatics" where
our group of three Alexander A., Christine W. and Myself were tasked with creating a python program
that would find open reading frames (ORFs) within sequences of a given FASTA file.'''



#Part 1 group project, Input Handling and Parsing (read and parse fasta file)
import re

# Main Author: Christine 
def fasta(seeit):    
#set variables to store the header and sequence and place them in a dictionary
    fasta_dict = {}
    header = ""
    sequence = ""
#open the file and read line by line searching for the header
    with open(seeit, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
#if the header line is found make it a key for the following sequence
                if header:
                    fasta_dict[header] = sequence
                    sequence = "" 
                header = line[1:].strip().upper() 
#if the header line is not found add the line to the sequence
            else:
                sequence += line.strip().upper() 
#if another header line is found set the previous header as the key and the sequence that followed it as the value
        if header:
            fasta_dict[header] = sequence    
    return fasta_dict




#Part 2 ORF identification

# Main Author: Alex
def make_antisense(inputStrand):  # Makes corresponding antisense DNA strand  
    dna_complement = ''  # Empty string to fill with complementary bases

    # Go nt by nt for the entire length of inputStrand
    for nt in inputStrand:

        # Add complementary bases to nascent dna_complement strand
        if nt == 'A':
            dna_complement += 'T'  # A pairs with T
        elif nt == 'T':
            dna_complement += 'A'
        elif nt == 'C':
            dna_complement += 'G'  # C pairs with G
        elif nt == 'G':
            dna_complement += 'C'

    # Return the string, but in reverse so its in the same orientation as the sense strand
    return dna_complement[::-1]

# Main Author: Alex
def orf_finder(inputStrand, fastaHeader):  # Find all ORFs in all six reading frames of the input strand

    # Make substrings ('frames') from sense strand by slicing the inputStrand
    senseFrames = [inputStrand, inputStrand[1:], inputStrand[2:]]

    # Now do the same for the antisense strand using make_antisense()
    antisenseStrand = make_antisense(inputStrand)  # Generating the antisense strand
    antisenseFrames = [antisenseStrand, antisenseStrand[1:], antisenseStrand[2:]]

    totalFrames = senseFrames + antisenseFrames  # Group all the framed substrings in one list

    # Hard-code start and stop codons, denoting start/stop signals for orf_finder
    startCodon = 'ATG'
    stopCodon = ['TAA', 'TAG', 'TGA']

    # Make a dictionary to store ORF information
    orfDictionary = {}

    # Make a counter to keep track of which frame is being scanned, for frame/position negation purposes
    frameNum = 1  # Starts by scanning the first item in the totalFrames list
    
    # Main for loop, iterate through all the frames in search of ORFs
    for frame in totalFrames:
        length = len(frame)  # Get the length of the substring
        startIndex = 0  # Starting at index 0, or the beginning of the sequence
        while startIndex <= length - 3:
            if frame[startIndex:startIndex + 3] == startCodon: # Scan nts in groups of three (codon), check for ATG
                for currentIndex in range(startIndex + 3, length, 3): # Scan subsequent codons for one of the valid stop codons
                     if frame[currentIndex:currentIndex + 3] in stopCodon:# If/when in-frame stop codon is found
                        # Compute ORF sequence & positional data based on start/stop codon indices
                        start = startIndex  # Position of start codon
                        stop = currentIndex + 3  # Position of stop codon
                        orf_sequence = frame[start:stop]  # Slice the substring using the start/stop codon indices
                        orf_length = len(orf_sequence)  # The length of the resultant sliced string

                        if frameNum > 3:  # Frames 4, 5, and 6 refer to the antisense strand
                            start = -(length - stop)

                             # Arrange the ORF data as a dictionary entry
                        orf_data = {
                            'Header': fastaHeader,  
                            'Frame': frameNum,
                            'Position': (start + 1),  # First nt in start codon, + 1 to account for zero-index offset
                            'Length': orf_length
                            }

                        orfDictionary[orf_sequence] = orf_data
                        startIndex=currentIndex +3
                        break
                else:
                    startIndex+=3
                    continue
            else:
                startIndex +=3
        frameNum += 1
    return orfDictionary


# Main Author: Matthew
def orf_format(orf_sequence):
        #create codons (triplets from all of the ORF sequence)
        codons = (orf_sequence[i:i+3] for i in range(0, len(orf_sequence), 3))
        #format codons by inputting spaces between codons
        format_codons = ' '.join(codons)
        #format codons further by creating a new line every 15 codons accounting for spaces
        format_codons = '\n'.join(format_codons[i:i+60] for i in range(0, len(format_codons), 60))
        #store formatted orf
        return format_codons


# Main Author: Matthew
def dict_sort(key, orfDictionary, orf_lenUser):
    #use the ORF sequence as a key to get its corresponding data from part2
    data = orfDictionary.get(key)
    if data:
    #if data is returned from the ORF format its sequence for output
        proper_orf = orf_format(key)
        final_formatting(data, proper_orf, orf_lenUser)
    
    return None
             


# Main Author: Matthew, Christine
def final_formatting(data, proper_orf, orf_lenUser):
    #make sure the orf is the proper length desired by the user
    if int(data['Length']) > orf_lenUser:
        #format and output header line of final output
        results = f">{data['Header']} | FRAME = {data['Frame']} POS = {data['Position']} LEN = {data['Length']}\n"
        #add sequence to final output
        results += proper_orf
        return print(results)
    else:
        return ''

# Main Author: Matthew
def main():
    print(
          """This program will find ORFs for sequences. When prompted, please enter a FASTA file with 1 or more sequences."""
          )
    
    seeit = input("Enter FASTA file: ")
    fasta_dict = fasta(seeit)
    main_orfDictionary={} 

    #from part 1, organize outputs: keys(headers) and values(sequences)
    header_keyList = fasta_dict.keys()

    #Use keys(headers) to generate a sequence and plug it into part 2
    for key in header_keyList:
        fasta_seq = fasta_dict[key]
        orfDictionary = orf_finder(fasta_seq, key)
        main_orfDictionary.update(orfDictionary)

    #output from part2 is 'orfDictionary'
    #set up part 2 to iterate through and read all values(sequences) and return them with the proper keys(headers)
    #part 2 outputs: keys(orf sequences) and values(     *another dictionary*  *keys(Header, Frame, Position, Length)  *values(headerValue, frameValue, posValue, lenValue)   )

    global orf_lenUser
    orf_lenUser = int(input("Enter minimum length in bp for ORFs: "))

     #remove all orfs not at proper length
    for key in main_orfDictionary:
        dict_sort(key, main_orfDictionary, orf_lenUser)


#test file EXsequence2.fasta
main()

