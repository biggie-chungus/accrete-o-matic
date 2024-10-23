#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:15:35 2024

@author: jonw
"""
import os
import pandas as pd
import re

def get_file_number(filename):
    # Extract the number from the filename
    match = re.search(r'_(\d+)\.csv$', filename)
    if match:
        return int(match.group(1))
    return None

def merge_csv_files(directory_path, output_file):
    # List to store individual dataframes and their file numbers
    dataframes = []

    # Iterate through all files in the directory
    for filename in os.listdir(directory_path):
        if filename == 'merged_file.csv':
            continue
        if filename.endswith(".csv"):
            file_number = get_file_number(filename)
            print (filename)
            if file_number is not None:
                file_path = os.path.join(directory_path, filename)
                # Read each CSV file into a dataframe
                df = pd.read_csv(file_path, low_memory=False)
                dataframes.append((file_number, df))

    # Sort dataframes based on the extracted file number
    dataframes.sort(key=lambda x: x[0])

    # Extract dataframes from the sorted list
    sorted_dataframes = [df for _, df in dataframes]

    # Concatenate all dataframes
    merged_df = pd.concat(sorted_dataframes, ignore_index=True)

    # Save the merged dataframe to a new CSV file
    merged_df.to_csv(output_file, index=False)




# Specify the directory containing the CSV files and the output file path
directory_path = '/Users/jonw/Desktop/earth_core/output/'
output_file = '/Users/jonw/Desktop/earth_core/output/merged_file.csv'

# Merge the CSV files
merge_csv_files(directory_path, output_file)
