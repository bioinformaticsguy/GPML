from pathlib import Path

def get_updated_path(original_path, new_subdirectory, new_extension):
    new_subdirectory = original_path.parent / new_subdirectory
    new_path_with_subdir = new_subdirectory / original_path.name
    new_path_with_extension = new_path_with_subdir.with_suffix(new_extension)

    return  new_path_with_extension

# MultiPreprocessing
import subprocess
import multiprocessing

# Function to run the command for a given input and output file
def run_command(input_file, output_file):
    command = f"./run_mutpred2.sh -i {input_file} -p 1 -c 1 -b 0 -t 0.1 -f 4 -o {output_file}"
    subprocess.call(command, shell=True)

if __name__ == "__main__":
    input_files = [file for file in multiprocessong_folder_path.rglob('*') if file.is_file()] * 36  # Assuming you want to use the same input file for all cores
    output_files = [get_updated_path(Path(orignal_path), "output_folder", ".out") for orignal_path in input_files] * 36  # Assuming you want to use the same output file for all cores

    pool = multiprocessing.Pool(processes=36)  # Number of cores
    pool.starmap(run_command, zip(input_files, output_files))
    pool.close()
    pool.join()
