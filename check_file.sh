#!/bin/bash

# Path to the folder
folder_path="/home/ma/dk/UCEC/raw/"

# Read names from the fourth column of the file and search for them in files within the folder
awk '{print $2}' ~/dk/UCEC/gdc_manifest.2023-08-13.txt | while IFS= read -r name; do
    # Search for the name in the contents of files within the folder
    if [ -e "$folder_path"/"$name" ]; then
        echo "$folder_path"/"$name: yes"
    else
        echo "$folder_path"/"$name: no"
    fi
done

