import subprocess
import multiprocessing

# Function to run the command for a given input and output file
def run_command(input_file, output_file):
    command = f"./run_mutpred2.sh -i {input_file} -p 1 -c 1 -b 0 -t 0.1 -f 4 -o {output_file}"
    subprocess.call(command, shell=True)

if __name__ == "__main__":
    input_files = ["t_36_seq.fasta"] * 36  # Assuming you want to use the same input file for all cores
    output_files = ["t_36_seq.out"] * 36  # Assuming you want to use the same output file for all cores

    pool = multiprocessing.Pool(processes=36)  # Number of cores
    pool.starmap(run_command, zip(input_files, output_files))
    pool.close()
    pool.join()

