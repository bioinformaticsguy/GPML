import os
from multiprocessing import Pool
from pathlib import Path

from src.parallel_script import run_command

input_folder = Path('Data/dbNSFP_input_files')  # Replace with your actual folder path
print(input_folder)
input_files = list(input_folder.glob('*'))
all_files = [file for file in input_folder.iterdir() if file.is_file()]

print(all_files)

num_processes = min(13, os.cpu_count())

with Pool(processes=num_processes) as pool:
    pool.map(run_command, input_files)
