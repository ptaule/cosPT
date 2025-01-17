import argparse
import numpy as np

# Helper function for measurements
def zscore(a_val, a_err, b_val, b_err):
    err = np.sqrt(a_err**2 + b_err**2)
    return (a_val - b_val)/(err)

def parse_commandline():
    parser = argparse.ArgumentParser(description="Compare columns in two files with optional error columns.")

    parser.add_argument("--col_A", type=int, default=1, help="Column in file A to compare (1-based index).")
    parser.add_argument("--col_err_A", type=int, help="Column for error in file A (1-based index).")
    parser.add_argument("--col_B", type=int, default=1, help="Column in file B to compare (1-based index).")
    parser.add_argument("--col_err_B", type=int, help="Column for error in file B (1-based index).")
    parser.add_argument("--verbose", "-v", action="store_true", help="Be verbose.")
    parser.add_argument("fileA", type=str, help="Path to file A.")
    parser.add_argument("fileB", type=str, help="Path to file B.")

    return parser.parse_args()

def main():
    args = parse_commandline()

    # Load data from files
    file_a = np.genfromtxt(args.fileA, comments="#")
    file_b = np.genfromtxt(args.fileB, comments="#")

    col_a = args.col_A - 1  # Convert to 0-based index
    col_err_a = args.col_err_A - 1 if args.col_err_A else None
    col_b = args.col_B - 1  # Convert to 0-based index
    col_err_b = args.col_err_B - 1 if args.col_err_B else None

    # Create measurements with or without errors
    a_errors = file_a[:, col_err_a] if col_err_a is not None else np.zeros(file_a.shape[0])
    b_errors = file_b[:, col_err_b] if col_err_b is not None else np.zeros(file_b.shape[0])

    a_vals = file_a[:, col_a]
    b_vals = file_b[:, col_b]

    if len(a_vals) != len(b_vals):
        print("Error: Column dimensions in file A and file B differ.")
        return

    print(f"Comparing column {col_a + 1} from {args.fileA} and column {col_b + 1} from {args.fileB}:")

    if np.array_equal(a_vals, b_vals):
        print("Columns are equal.")
        return

    # Compute z-scores
    zscores = zscore(a_vals, a_errors, b_vals, b_errors)

    for i, z in enumerate(zscores):
        if args.verbose or z > 1:
            print(f"Row {i + 1}: a = {a_vals[i]}, b = {b_vals[i]}, zscore = {z}")

    print("done.")

if __name__ == "__main__":
    main()
