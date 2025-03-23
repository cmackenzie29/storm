import sys
import random
from storm_analysis import *

def binary_input(message, input1, input2):
    while True:
        user_input = input(message)
        if (user_input == input1) or (user_input == input2):
            return user_input

def float_input(message):
    while True:
        user_input = input(message)
        try:
            user_input = float(user_input)
            if user_input <= 0:
                continue
            return user_input
        except:
            continue

def get_user_inputs():
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

    # Create results folder if it does not exist
    if not os.path.isdir(os.path.join(folder, "results")):
        os.mkdir(os.path.join(folder, "results"))
    
    os.system("clear")

    dim = int(binary_input("Is this data 2D STORM or 3D STORM (type 2/3): ", "2", "3"))

    while True:
        os.system("clear")
        red_name = input("Enter the name of the molecule tagged by red (or press enter to exclude): ")
        green_name = input("Enter the name of the molecule tagged by green (or press enter to exclude): ")
        magenta_name = input("Enter the name of the molecule tagged by magenta/far-red (or press enter to exclude): ")
        os.system("clear")
        print("\nRed: " + ("-SKIP-" if red_name=="" else "\""+red_name+"\""))
        print("Green: " + ("-SKIP-" if green_name=="" else "\""+green_name+"\""))
        print("Magenta: " + ("-SKIP-" if magenta_name=="" else "\""+magenta_name+"\""))
        confirmation = binary_input("\nIs this correct (enter y/n): ", "y", "n")
        if confirmation == "y":
            break

    colors = {}
    if red_name != "":
        colors[255] = red_name
    if green_name != "":
        colors[65280] = green_name
    if magenta_name != "":
        colors[16711935] = magenta_name

    os.system("clear")
    print("How would you like your files analyzed?\n")
    print("Clustering mode looks for clusters of molecules using DBSCAN and provides a more thorough analysis.\n")
    print("Basic mode provides nearest neighbor distances between molecules, regardless of if your files have clusters.\n")
    mode = binary_input("Proceed in clustering mode or basic mode (type c/b): ", "c", "b")

    os.system("clear")
    output_mode = binary_input("Would you like to output only summary statistics or raw data entries as well (type s/r): ", "s", "r")
    os.system("clear")

    if mode == "b":
        # Basic mode
        print(f"Analyzing files with prefix \"{condition_prefix}\" in basic mode...\n")
        analyze_basic(folder, selected_files, condition_prefix, colors, output_mode)
        sys.exit()

    # Clustering mode
    clustering_mode = binary_input("Do the clusters contain all molecule colors together or just one color (type a/o): ", "a", "o")

    eps = 0.5
    min_samples = 10

    positions, storm_colors, _ = load_storm_file(folder, random.choice(selected_files), colors)

    while True:
        os.system("clear")
        print(f"DBSCAN parameters: eps = {eps}, min_samples = {min_samples}. Molecules are shown in grey, with clusters shown in colors. Close image to continue...\n")
        if clustering_mode == "a":
            display_clusters(positions, eps, min_samples)
        else:
            display_clusters(positions, eps, min_samples, storm_colors)
        confirmation = binary_input("Does the clustering look good (type y/n): ", "y", "n")
        if confirmation == "y":
            break
        else:
            print("\nEnter new DBSCAN parameters to try: ")
            eps = float_input("eps = ")
            min_samples = round(float_input("min_samples = "))

    os.system("clear")
    print(f"Analyzing files with DBSCAN parameters: eps = {eps}, min_samples = {min_samples}...")
    if clustering_mode == "a":
        analyze_clustering_colors_together(folder, files, condition_prefix, dim, colors, eps, min_samples, output_mode)
    else:
        analyze_clustering_colors_individual(folder, files, condition_prefix, dim, colors, eps, min_samples, output_mode)

get_user_inputs()
