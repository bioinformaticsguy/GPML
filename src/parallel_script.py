import os
from multiprocessing import Pool
from pathlib import Path
import subprocess

OUTPUT_FOLDER = Path("Data/dbNSFP_output_folder")


def run_command(input_file):
    base_name = input_file.stem
    output_file = OUTPUT_FOLDER / f"{base_name}_output.csv"
    print(output_file)

    # command = ["java", "search_dbNSFP44a", "-i", str(input_file), "-o", output_file]
    # command = ["touch", output_file]
    # subprocess.run(command)
