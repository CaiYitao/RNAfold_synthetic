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

    def generate_length():
        if distribution == "normal":
            return max(
                min_length,
                min(
                    max_length,
                    int(
                        random.gauss(
                            kwargs.get("mean", mean_length),
                            kwargs.get("std", std_dev_length),
                        )
                    ),
                ),
            )
        elif distribution == "uniform":
            return random.randint(min_length, max_length)
        elif distribution == "exponential":
            return max(
                min_length,
                min(max_length, int(random.expovariate(kwargs.get("lambda", 0.01)))),
            )
        elif distribution == "log-normal":
            return max(
                min_length,
                min(
                    max_length,
                    int(
                        random.lognormvariate(
                            kwargs.get("mu", 6), kwargs.get("sigma", 0.5)
                        )
                    ),
                ),
            )
        elif distribution == "custom":
            return random.choice(kwargs.get("custom_lengths", [min_length, max_length]))
        else:
            raise ValueError("Unsupported distribution.")

    # Initialize the ViennaRNA fold model
    md = RNA.md()

    # Data storage
    data = []

    # Generate the dataset
    for i in range(num_samples):
        # Generate sequence length
        n = generate_length()
        # Generate random RNA sequence
        sequence = RNA.random_string(n, "ACGU")

        # Compute secondary structure and MFE
        fc = RNA.fold_compound(sequence, md)
        structure, mfe = fc.mfe()

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


# rna_sequence = "CCGGUACGAUACUUAAUAAAUUUUGAUGUGUUUGGGCGUAUGCGGUUGGACUUCUAGGUGCGAAGGGGGAUGGUUCUGUCACCUGGUUCUGAAGAUCAGCUUAGAUCCCCCUGCCCUUGCCUGACACUGCGGAACAUAAUGUUCUUAUCUGAAUCCACUAUAGCGAAGUAGCAGCCCCGAGACCCAAUCAGGUUCAUUGUGGGAUUGUUCACCCGUUAGUGUCGGUCACAAUUAUCUAGACCCACGCUCCGCGCAAUUGAAGAACGGAUCCUCCGCAGAAUACAAAACUUCAACACUGCAAUUACCAUUUGCUUAAUGAUUCAAUACAUUUCACAAUGAUGACGUUACAAAUCCAAUGAAUAAUUCUAGUGC"
# dot_bracket_structure = "..((((.(((((((((......)))).))))).((((....(((.((((....)))).)))..(((((((....((((....((((((.....))))))..))))))))))))))).)))).....(((((((.......((((......))))......(((((..(((((((((((((((.((.....)).)).))).))).))))).)).)))))((((.((((...........)))).))))((((..(.......)..))))...)))))))...............(((((.(((((.((((.(.((..(((((((...((((....)))).))).))))..)).).))))..))))).)))))."
# mfe = -88.0999984741211

# ct_content = rna_to_ct(rna_sequence, dot_bracket_structure, mfe)
# print(ct_content)

if __name__ == "__main__":
    num_samples = 10
    mean_length = 88
    std_dev_length = 66
    output_ct = f"rna_ct_files_{mean_length}"

    df = generate_rna_dataset(
        num_samples=num_samples,
        mean_length=mean_length,
        std_dev_length=std_dev_length,
        output_file=f"rna_dataset_{mean_length}.csv",
        output_ct=output_ct,  # Directory for .ct files
        distribution="normal",
    )
#     print(df)
