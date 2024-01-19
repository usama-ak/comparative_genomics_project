#!/usr/bin/python3

import os

# Function to read files and merge data
def merge_files(file1_path, file2_path, output_path, a, b):
    file2_data = {}
    
    with open(file2_path, 'r') as file2:
        for line in file2:
            columns = line.strip().split()
            if len(columns) >= 3:
                key = columns[0]
                value2 = columns[1]
                value3 = columns[2]
                file2_data[key] = (value2, value3)

    with open(file1_path, 'r') as file1, open(output_path, 'w') as output_file:
        for line in file1:
            columns = line.strip().split()
            if len(columns) >= a:
                key = columns[b]
                if key in file2_data:
                    value2, value3 = file2_data[key]
                    output_line = f"{line.strip()}\t{value2}\t{value3}\n"
                    output_file.write(output_line)
                else:
                    output_file.write(f"{line.strip()} \t \n")
                    

file1_path = 'Medicago_truncatula_Blastp_longIsoforme'
file2_path = 'Medicago_truncatula_protein_length'
file3_path = 'temp'
output_file= 'output.txt'


merge_files(file1_path, file2_path, file3_path, 1,0)
merge_files(file3_path, file2_path, output_file, 2,1)



def filter_hits(file1_path, output_path):
    mydict = {}

    # Read file1 and store bitscore values in a dictionary for both directions
    with open(file1_path, 'r') as file1:
        for line in file1:
            columns = line.strip().split()
            key1 = columns[0]
            key2 = columns[1]
            bitscore = float(columns[11])

            # Consider both directions: (key1, key2) and (key2, key1)
            mydict[(key1, key2)] = max(mydict.get((key1, key2), -float('inf')), bitscore)
            mydict[(key2, key1)] = max(mydict.get((key2, key1), -float('inf')), bitscore)

    # Read file1 again, filter data based on bitscore, and write the filtered data to the output file
    with open(file1_path, 'r') as file1, open(output_path, 'w') as output_file:
        for line in file1:
            columns = line.strip().split()
            key = (columns[0], columns[1])

            # Retrieve the maximum bitscore for both directions
            max_bitscore = max(mydict.get((key[0], key[1]), -float('inf')),
                              mydict.get((key[1], key[0]), -float('inf')))

            if float(columns[11]) == max_bitscore:  # Only write if it's the maximum bitscore
                output_file.write(f"{line.strip()} \t \n")
	    

filter_hits('output.txt', 'filtered_data.txt')

os.remove(file3_path)
os.remove("output.txt")
print("The output file is saved in 'filtered_data.txt'.")
			
			
			


