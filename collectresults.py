"""
Code to collect results from multiple experiments and generate figures and tables.

Should NOT be run with a dirty working directory! Will travel to multiple git commits.
"""
import os
import subprocess
import tempfile
import shutil
import csv
# Don't have git, so use subprocess to call git commands.


def run_experiment(commit_hash, experiment_args, repeats = 3):
    """
    Run the experiment at the given commit hash with specified arguments and number of repeats.
    Returns the array of results.
    """
    # Results are scraped from the line with "Total time: <time_in_seconds_float> seconds".
    # checkout, make clean, make, "./main.exe <args>", get results.
    subprocess.run(['git', 'checkout', commit_hash], check=True)
    subprocess.run(['make', 'clean'], check=True)
    subprocess.run(['make'], check=True)
    results = []
    print("Running repeats...", end=" ")
    for i in range(repeats):
        print(i + 1, end=" ", flush=True)
        completed_process = subprocess.run(['./main.exe'] + experiment_args, capture_output=True, text=True, check=True)
        for line in completed_process.stdout.splitlines():
            line = line.strip()
            if line.startswith("Total time:"):
                time_str = line.split("Total time:")[1].strip().split(" ")[0]
                results.append(float(time_str))
                break
    print()
    return results

commits = [
    # ("Initial", "3d0550609610e188a2ea1bbfd10b9345b8e73dc0"),
    ("Predict dE", "83305cea551edf6f0b791d77ab53b8a8fff7cbef"),
    ("Reduce Storage", "fd8e2e7b76978e594edb6703e5cc566c87bb2f45"),
    ("Return dE Directly", "88dd48dc02ae052be21034ce08460ac560b19518"),
    ("Naive Parallelism", "8c6b7940ed3d1c518fc6dc71622d6ee4b9c967cd"),
    ("Parallel Rand (Without Balance)", "5b09082cc7e0cd91cf1b01c31f14f123fa596b69"),
    ("Parallel Rand (Load Balancing)", "fc4d7eeafd514e4c15a9bef1ffcbbb0e293ccf3e")
]

def avg(data):
    return sum(data) / len(data)



if __name__ == "__main__":
    all_results = []
    # Check we are in the "host" directory.
    cwd = os.getcwd()
    if not os.path.exists(os.path.join(cwd, "main.cpp")):
        print("WARNING: this script should be run from within /host of the docker image!")
    else:
        subprocess.run(['git', 'config', '--global', '--add', 'safe.directory', cwd], check=True)

    for name, commit_hash in commits:
        avg_small = avg(run_experiment(commit_hash, ["100", "2.269", "1000000", "1"]))
        
        n_its_large = 100000
        n_its_large_long = 1000000000
        if name == "Initial":
            # Too long otherwise.
            n_its_large = n_its_large / 100
            n_its_large_long = n_its_large_long / 100000
        avg_large = avg(run_experiment(commit_hash, ["1024", "2.269", str(int(n_its_large)), "1"]))
        avg_large_long = avg(run_experiment(commit_hash, ["1024", "2.269", str(int(n_its_large_long)), "1"]))
        if name == "Initial":
            # Scale back up to be comparable.
            avg_large = avg_large * 100
            avg_large_long = avg_large_long * 100000

        print(f"Results for {name}:")
        print(f"  Small (1M its): {avg_small:.4f} seconds")
        print(f"  Large (100M its): {avg_large:.4f} seconds")
        print(f"  Large Long (1B its): {avg_large_long:.4f} seconds")
        all_results.append((name, avg_small, avg_large, avg_large_long))


        # Create a summary table
        baseline = all_results[0][1:]
        previous = baseline
        with open("results_summary.csv", "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Name", "Small (1M its)", "Cumulative Speedup", "Individual Speedup", "Large (100M its)", "Cumulative Speedup", "Individual Speedup", "Large Long (1B its)", "Cumulative Speedup", "Individual Speedup"])
            for result in all_results:
                row = [result[0]]
                times = result[1:]
                for i in range(len(times)):
                    time = times[i]
                    row.append(f"{time:.4f}")
                    cum_speedup = baseline[i] / time
                    indiv_speedup = previous[i] / time
                    row.append(f"{cum_speedup:.2f}x")
                    row.append(f"{indiv_speedup:.2f}x")
                previous = times
                writer.writerow(row)