import sys

def generate_ranges(start, label, stop, fname="guide_file_", step=100000):
    start = int(start)  # Ensure start is an integer
    stop = int(stop)  # Ensure stop is an integer
    filename = fname + label + ".txt"
    
    with open(filename, "w") as f:
        while start + step <= stop:
            end = start + step
            range_col = f"{label}:{start}-{end-1}"
            f.write(f"{start}\t{end}\t{label}\t{step}\t{end-1}\t{start}\t{range_col}\n")
            start = end  # Update start for the next row

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <start> <label> <stop>")
    else:
        start_value = sys.argv[1]
        label_value = sys.argv[2]
        stop_value = sys.argv[3]
        generate_ranges(start_value, label_value, stop_value)

