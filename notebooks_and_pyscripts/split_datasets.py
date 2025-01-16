from casatools import ms
import argparse

CLI = argparse.ArgumentParser()

CLI.add_argument("--idataset_name", default="in.ms")
CLI.add_argument("--num_pix", type=int, default=1)
CLI.add_argument("--pix_rad", type=float, default=1)
CLI.add_argument("--halfwidth", type=int, default=1)
CLI.add_argument("--centers", nargs="*", type=int, default=[1, 2, 3])
CLI.add_argument("--odataset_name", default="out")

args = CLI.parse_args()

full_ms_filename = args.idataset_name
num_pix = args.num_pix
pix_size_rad = args.pix_rad
hw = args.halfwidth
centers = args.centers
o_ms_name = args.odataset_name

centers.sort()

max_pix = 2 * num_pix

m = ms()
m.open(full_ms_filename)

for i, band in enumerate(centers):
	band_low = 0 if i == 0 else centers[i - 1] - hw
	band_high = band + hw

	band_low_lambda = band_low / (num_pix * pix_size_rad)
	band_high_lambda = band_high / (num_pix * pix_size_rad)

	curr_band_filename = o_ms_name + str(band) + ".ms"

	m.split(outputms=curr_band_filename,uvrange=str(band_low_lambda) + "~" + str(band_high_lambda) + "lambda")


band_low = centers[-1]
band_high = max_pix

band_low_lambda = band_low / (num_pix * pix_size_rad)
band_high_lambda = band_high / (num_pix * pix_size_rad)

curr_band_filename = o_ms_name + str(band_low) + "plus.ms"

m.split(outputms=curr_band_filename,uvrange=str(band_low_lambda) + "~" + str(band_high_lambda) + "lambda")

m.close()