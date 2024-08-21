from collections import Counter

# Read the isolate names from the file
with open('isolates.txt', 'r') as infile:
    lines = infile.readlines()

# Extract the sample names (ignoring the country names)
sample_names = [line.strip().split('\t')[1] for line in lines]

# Count occurrences of each sample name
name_counts = Counter(sample_names)

# Identify duplicates
duplicates = {name: count for name, count in name_counts.items() if count > 1}

# Output the results
if duplicates:
    print(f"Found {len(duplicates)} duplicate sample names.")
    for name, count in duplicates.items():
        print(f"Sample name '{name}' appears {count} times.")
else:
    print("No duplicate sample names found.")
