#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3
# Edit the shebang line for your OS
# collate_seq_files.py
# Takes a set of .seq files and collates their contents into a single .fasta
# file

import os, Bio, re, sys, datetime # Dependencies
from Bio.Seq import Seq
import Bio.Alphabet
from Bio.Alphabet import IUPAC

# Change directory to where seq files are stored
print('''Enter the ABSOLUTE path to the directory containing your .seq files.
If you are already in that directory, hit 'Enter':''')
directory = input()
if directory == '':
    directory = '.'
os.chdir(directory)

### Look in directory for seq files
# Regex to get only the names of .seq files in the directory
seq_regex = re.compile(r'([a-zA-Z0-9_-]*)(\.seq)')

# convert the directory contents to a list
file_list = os.listdir()
file_str = ', '.join(file_list)

# find all the .seq files
seq_matches = seq_regex.findall(file_str)
if seq_matches == []:
    print('There are no .seq files in ' + str(os.getcwd()))
    sys.exit()
    
# pass the names of each file to a new list, so they can be used for headers
# in the collated .fasta file
name_list = []
for i in range(len(seq_matches)):
    name_list.append(seq_matches[i][0])

# Replaces dash characters with underscores
# This code can be adapted to replace any other characters in file headers that
# are messing with the BLAST algorithm
name_list = [w.replace('-','_') for w in name_list]

### Loop through each seq file, read its contents and append it to a single
### fasta file with an appropriate header
print('Enter a name for the collated .fasta file (remember to escape any special characters with \):')
fasta_name = input()
now = datetime.datetime.now()
now = (str(now.day) + '-' + str(now.month) + '-' + str(now.year))
fasta_file = open('%s_%s.fasta' % (fasta_name, now), 'w')

sequence = '' # initialise a string to hold the nucleotide sequence
i = 0         # initialise i at 0 to act as a counter 

# loop through the cwd           
for file in os.listdir(directory):

    if file.endswith('.seq'):
        fasta_file.write('>%s\n' % name_list[i])       #Adds the filename as a FASTA header
        sequence = open('%s' % file)                    #Opens the file
        sequence = ''.join(sequence.readlines())        #Lines are read from the file, and joined into a string
                                                        #to make it useable
        sequence = Seq(sequence, IUPAC.ambiguous_dna)   #Converts the file contents to an IUPAC format
        fasta_file.write(str(sequence))                 #The sequence is written under the header
        i += 1                                          #Increases i by 1, so that name_list[i] points at the
                                                        #correct index
            


    

