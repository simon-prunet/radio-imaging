import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

lines = None

with open(input_file) as file:
    lines = [line.replace(",", ".") if "USER" not in line else line.replace("kB_rd/s", "0").replace("kB_wr/s", "0") for line in file if 'Linux' not in line and line != '\n']

with open(output_file, 'w') as f:
    for line in lines:
        f.write(f"{line}")
