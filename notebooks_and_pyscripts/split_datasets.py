from casatools import ms

from casatools import ms
import sys

full_ms_filename = sys.argv[1]
low_ms_filename = sys.argv[2]
high_ms_filename = sys.argv[3]
lambda_low = sys.argv[4]
lambda_high = sys.argv[5]

m = ms()
m.open(full_ms_filename)
m.split(outputms=low_ms_filename,uvrange="<" + lambda_low + "lambda")
m.split(outputms=high_ms_filename,uvrange=">" + lambda_high + "lambda")
m.close()