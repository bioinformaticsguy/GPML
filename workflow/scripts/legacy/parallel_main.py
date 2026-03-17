import os
from multiprocessing import Pool

from pathlib import Path
import subprocess

OUTPUT_FOLDER = Path("dbNSFP_output_folder")


def run_command(input_file):
    base_name = input_file.stem
    output_file = OUTPUT_FOLDER / f"{base_name}_output.csv"

    command = ["java", "search_dbNSFP45a", "-i", str(input_file), "-o", output_file]
    # command = ['bash', '-c', f'echo "{content}" > {output_file}']
    subprocess.run(command)


input_folder = Path('dbNSFP_input_files')  # Replace with your actual folder path
print(input_folder)
input_files = list(input_folder.glob('*'))
print(input_files)
# all_files = [file for file in input_folder.iterdir() if file.is_file()]

# print(all_files)

num_processes = min(3, os.cpu_count())




with Pool(processes=num_processes) as pool:
    pool.map(run_command, input_files)
