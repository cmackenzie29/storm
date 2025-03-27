import os
import sys
import numpy as np
import pandas as pd

### Lists the number corresponding to the colors in the txt storm files

def binary_input(message, input1, input2):
    while True:
        user_input = input(message)
        if (user_input == input1) or (user_input == input2):
            return user_input

def get_storm_colors(folder, f):
    directory = os.path.join(folder, f)
    with open(directory, 'rb') as f:
        contents = f.read()
    rows = str(contents,'utf-16').split('\r\n')[1:]
    storm_data = pd.DataFrame(data = [row.split('\t') for row in rows[1:(len(rows)-1)]], columns=rows[0].split('\t'))
    
    return np.unique(np.array(storm_data.Color).astype(int))


os.system("clear")

while True:
    folder = input("Enter name of folder or path to directory containing .txt files: ")
    os.system("clear")
    if os.path.isdir(folder):
        break
    else:
        print(f"Could not find \"{folder}\". Current directory: \"{os.getcwd()}\"\n")

files = [x for x in os.listdir(folder) if ".txt" in x]
if len(files) < 1:
    print("No .txt files found. Try again with a different directory containing .txt files.")
    sys.exit()

while True:
    print(f"Files in \"{folder}\":\n" + "\n".join(files) + "\n")
    condition_prefix = input(f"Enter a prefix unique to the files for a single condition: ")
    os.system("clear")

    selected_files = [x for x in files if x[:len(condition_prefix)]==condition_prefix]
    if len(selected_files) > 0:
        print(f"Files in \"{folder}\" with prefix \"{condition_prefix}\":\n" + "\n".join(selected_files) + "\n")
        confirmation = binary_input("Proceed with analyzing these files together (type y) or enter a new condition prefix (type n)? ", "y", "n")
        if confirmation == "y":
            break
    else:
        confirmation = binary_input(f"No files in \"{folder}\" found matching prefix \"{condition_prefix}\". Enter a different condition prefix (type y) or quit (type q)? ", "y", "q")
        os.system("clear")
        if confirmation == "q":
            sys.exit()

os.system("clear")

all_colors = np.zeros(0)
for f in selected_files:
	all_colors = np.append(all_colors, get_storm_colors(folder, f))

print(f"Colors found in files: {np.unique(all_colors).astype(int)}")

