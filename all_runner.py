import subprocess

# List of your scripts
scripts = ["main_dataframe_preprocessor.py",
           "main_pssm_baseline.py",
           "main_corelation_calculator.py",
           "main_plot_graphs.py"]

for script in scripts:
    subprocess.call(["python", script])


print("Debug Pause")