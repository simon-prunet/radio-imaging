import casatasks
import argparse

CLI = argparse.ArgumentParser()

CLI.add_argument("--idatasets", default=["test.ms"], nargs='+')
CLI.add_argument("--odataset_name", default="out.ms")

args = CLI.parse_args()

input_ms = args.idatasets
output_ms = args.odataset_name

print(input_ms)

casatasks.concat(vis=input_ms, concatvis=output_ms)