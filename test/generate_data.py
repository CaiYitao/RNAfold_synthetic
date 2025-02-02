import RNA
import random
import pandas as pd
import os


def generate_rna_dataset(
    num_samples=100,
    mean_length=88,
    std_dev_length=66,
    output_file="rna_dataset.csv",
    output_ct=None,
    distribution="normal",
    **kwargs,
):
    """
    Generate a dataset of random RNA sequences and their secondary structures.

    Parameters:
        num_samples (int): Number of RNA samples to generate.
        mean_length (int): Mean of the RNA sequence length distribution.
        std_dev_length (int): Standard deviation of the RNA sequence length distribution.
        output_file (str): Path to the output CSV file.
        output_ct (str): Directory to save .ct files (if None, .ct files are not generated).
        distribution (str): Distribution to sample sequence lengths.
        **kwargs: Additional parameters for custom distributions.

    Returns:
        pd.DataFrame: A DataFrame containing the dataset.
    """

    # Calculate bounds dynamically based on mean_length and std_dev_length
    min_length = max(
        1, mean_length - 2 * std_dev_length
    )  # Ensure min_length is positive
    max_length = mean_length + 2 * std_dev_length

    # Initialize the ViennaRNA fold model
    md = RNA.md()

    # Data storage
    data = []

    # Generate the dataset
    for i in range(num_samples):
        # Generate sequence length
        # n = generate_length(distribution, min_length, max_length, **kwargs)
        # Generate random RNA sequence
        sequence = generate_sequence()

        # Compute secondary structure and MFE
        structure, mfe = compute_structures(sequence, subopt=False)

        # Append to data list
        data.append({"Sequence": sequence, "Structure": structure, "MFE": mfe})

        # If output_ct is specified, save .ct file
        if output_ct:
            if not os.path.exists(output_ct):
                os.makedirs(output_ct)
            ct_filename = os.path.join(output_ct, f"rna_sample_{i + 1}.ct")
            with open(ct_filename, "w") as ct_file:

                ct_data = rna_to_ct(sequence, structure, mfe)
                ct_file.write(ct_data)
                ct_file.close()

    # Create a DataFrame
    df = pd.DataFrame(data)

    # Save to CSV
    path = os.path.join(os.getcwd(), output_file)
    df.to_csv(path, index=False)
    print(f"Dataset saved to {path}")

    if output_ct:
        print(f".ct files saved to {os.path.abspath(output_ct)}")

    return df


def rna_to_ct(rna_sequence, dot_bracket_structure, mfe):
    """
    Convert RNA sequence data and dot-bracket structure into a CT format string.

    Parameters:
        rna_sequence (str): RNA sequence.
        dot_bracket_structure (str): Secondary structure in dot-bracket notation.
        mfe (float): Minimum Free Energy value.

    Returns:
        str: CT file content as a string.
    """
    # Validate inputs
    if len(rna_sequence) != len(dot_bracket_structure):
        raise ValueError(
            "The RNA sequence and dot-bracket structure must have the same length."
        )

    n = len(rna_sequence)
    # Initialize base-pairing information with all zeros (unpaired)
    base_pairs = [0] * n
    stack = []

    # Parse dot-bracket structure to get base pairings
    for i, symbol in enumerate(dot_bracket_structure):
        if symbol == "(":
            stack.append(i)
        elif symbol == ")":
            if not stack:
                raise ValueError(
                    "Unmatched closing parenthesis in dot-bracket structure."
                )
            j = stack.pop()
            base_pairs[i] = j + 1
            base_pairs[j] = i + 1

    if stack:
        raise ValueError("Unmatched opening parenthesis in dot-bracket structure.")

    # Create CT lines
    ct_lines = [f"{n} {mfe:.2f}"]  # Header line: number of nucleotides and MFE
    for i in range(n):
        line = f"{i + 1} {rna_sequence[i]} {i} {i + 2 if i + 1 < n else 0} {base_pairs[i]} {i + 1}"
        ct_lines.append(line)

    return "\n".join(ct_lines)


import RNA
import random
import csv
from multiprocessing import Pool, cpu_count
from functools import partial


def generate_length(distribution, min_length, max_length, **kwargs):
    if distribution == "normal":
        mean = kwargs["mean"]
        std = kwargs["std"]
        length = int(random.gauss(mean, std))
    elif distribution == "uniform":
        length = random.randint(min_length, max_length)
    elif distribution == "exponential":
        lambd = kwargs["lambda"]
        length = int(random.expovariate(lambd))
    elif distribution == "log-normal":
        mu = kwargs["mu"]
        sigma = kwargs["sigma"]
        length = int(random.lognormvariate(mu, sigma))
    else:
        raise ValueError(f"Unsupported distribution: {distribution}")
    # Ensure the length is within min and max
    return max(min_length, min(max_length, length))


def generate_sequence():
    params = distribution_params.get(distribution, {})
    length = generate_length(
        distribution=distribution,
        min_length=min_length,
        max_length=max_length,
        **params,
    )
    return RNA.random_string(length, "ACGU")


def compute_structures(sequence, subopt=True):
    md = RNA.md()
    fc = RNA.fold_compound(sequence, md)

    # Compute suboptimal structures within delta_energy from MFE
    if subopt:
        subopt_solutions = fc.subopt(delta_energy)
        return subopt_solutions
    # return mfe_structure, mfe_energy, subopt_solutions
    else:
        mfe_structure, mfe_energy = fc.mfe()
        return mfe_structure, mfe_energy


# def classify_sequence(sequence, mfe_structure, mfe_energy, subopt_solutions):
def classify_sequence(sequence, subopt_solutions):
    if len(subopt_solutions) < 2:
        return None  # Not enough structures to determine d
    # The first solution is MFE, the second is next best
    mfe = subopt_solutions[0].energy
    mfe_structure = subopt_solutions[0].structure
    second_best = subopt_solutions[1]
    d = abs(mfe - second_best.energy)
    data = {
        "sequence": sequence,
        "mfe_structure": mfe_structure,
        "mfe_energy": mfe,
        "second_best_structure": second_best.structure,
        "second_best_energy": second_best.energy,
        "d": d,
    }
    if d > theta1:
        return ("easy", data)
    elif theta2 <= d <= theta1:
        return ("medium", data)
    else:
        return ("hard", data)


# def process_sequence(_):
#     try:
#         seq = generate_sequence()
#         solutions = compute_structures(seq)
#         result = classify_sequence(solutions)
#         if result:
#             result[1]["sequence"] = seq  # Add raw sequence
#             return result
#     except Exception as e:
#         return None


# def generate_dataset_complex(num_samples):
#     datasets = {"easy": [], "medium": [], "hard": []}

#     with Pool(cpu_count()) as pool:
#         # Process sequences in parallel batches
#         batch_size = 100
#         needed = 3 * num_samples

#         while sum(map(len, datasets.values())) < needed:
#             results = pool.map(process_sequence, range(batch_size))

#             for result in results:
#                 if not result:
#                     continue
#                 class_label, data = result
#                 if len(datasets[class_label]) < num_samples:
#                     datasets[class_label].append(data)

#                 # Early exit if all classes filled
#                 if all(len(v) >= num_samples for v in datasets.values()):
#                     return datasets
#     return datasets


def generate_dataset_complex(num_samples):
    datasets = {"easy": [], "medium": [], "hard": []}
    attempts = 0
    max_attempts = 3 * num_samples**2  # Prevent infinite loops
    while (
        len(datasets["easy"]) < num_samples
        or len(datasets["medium"]) < num_samples
        or len(datasets["hard"]) < num_samples
    ) and attempts < max_attempts:
        sequence = generate_sequence()
        # mfe_structure, mfe_energy, subopt_solutions = compute_structures(sequence)
        subopt_solutions = compute_structures(sequence)

        # print(mfe_structure, mfe_energy, subopt_solutions)
        classification = classify_sequence(sequence, subopt_solutions)
        attempts += 1
        if classification is None:
            continue
        class_label, data = classification
        if len(datasets[class_label]) < num_samples:
            datasets[class_label].append(data)
    print("attempts:", attempts)
    return datasets


def save_to_csv(datasets, prefix="rna_dataset"):
    for class_label, data in datasets.items():
        filename = f"{prefix}_{class_label}.csv"
        with open(filename, "w", newline="") as csvfile:
            if not data:
                continue  # Skip empty datasets
            writer = csv.DictWriter(csvfile, fieldnames=data[0].keys())
            writer.writeheader()
            writer.writerows(data)


# rna_sequence = "CCGGUACGAUACUUAAUAAAUUUUGAUGUGUUUGGGCGUAUGCGGUUGGACUUCUAGGUGCGAAGGGGGAUGGUUCUGUCACCUGGUUCUGAAGAUCAGCUUAGAUCCCCCUGCCCUUGCCUGACACUGCGGAACAUAAUGUUCUUAUCUGAAUCCACUAUAGCGAAGUAGCAGCCCCGAGACCCAAUCAGGUUCAUUGUGGGAUUGUUCACCCGUUAGUGUCGGUCACAAUUAUCUAGACCCACGCUCCGCGCAAUUGAAGAACGGAUCCUCCGCAGAAUACAAAACUUCAACACUGCAAUUACCAUUUGCUUAAUGAUUCAAUACAUUUCACAAUGAUGACGUUACAAAUCCAAUGAAUAAUUCUAGUGC"
# dot_bracket_structure = "..((((.(((((((((......)))).))))).((((....(((.((((....)))).)))..(((((((....((((....((((((.....))))))..))))))))))))))).)))).....(((((((.......((((......))))......(((((..(((((((((((((((.((.....)).)).))).))).))))).)).)))))((((.((((...........)))).))))((((..(.......)..))))...)))))))...............(((((.(((((.((((.(.((..(((((((...((((....)))).))).))))..)).).))))..))))).)))))."
# mfe = -88.0999984741211

# ct_content = rna_to_ct(rna_sequence, dot_bracket_structure, mfe)
# print(ct_content)

if __name__ == "__main__":

    import time

    start_time = time.time()
    num_samples = 1666
    mean_length = 198
    std_dev_length = 66
    output_ct = f"rna_ct_files_{mean_length}"
    distribution_params = {
        "normal": {"mean": mean_length, "std": std_dev_length},
        "uniform": {},
        "exponential": {"lambda": 0.01},
        "log-normal": {"mu": 6, "sigma": 0.5},
        # "custom": {"custom_lengths": [min_length, max_length]}
    }
    distribution = "normal"
    min_length = max(1, mean_length - 2 * std_dev_length)
    max_length = mean_length + 3 * std_dev_length

    # Generate RNA dataset
    df = generate_rna_dataset(
        num_samples=num_samples,
        mean_length=mean_length,
        std_dev_length=std_dev_length,
        output_file=f"rna_dataset_{mean_length}.csv",
        output_ct=output_ct,  # Directory for .ct files
        distribution="normal",
    )
    #     print(df)

    # Configuration parameters
    # theta1 = 1.6  # Threshold for easy datasets
    # theta2 = 0.5  # Threshold for medium datasets
    # delta_energy = 500  # Energy delta for suboptimal structures
    # num_samples_per_class = 100  # Number of samples per class

    # # Sequence length parameters
    # mean_length = 200
    # std_dev_length = 68
    # min_length = max(1, mean_length - 2 * std_dev_length)
    # max_length = mean_length + 3 * std_dev_length
    # distribution = "normal"  # Choose from "normal", "uniform", "exponential", "log-normal", "custom"
    # # Parameters for each distribution type

    # # Generate and save datasets
    # datasets = generate_dataset_complex(num_samples_per_class)
    # save_to_csv(datasets)

    # print(
    #     f"Generated {len(datasets['easy'])} easy, {len(datasets['medium'])} medium, {len(datasets['hard'])} hard samples."
    # )
    # print the time cost
    print("Time cost: ", time.time() - start_time)
