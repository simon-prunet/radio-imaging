#!/bin/bash
wget https://cloud.oca.eu/index.php/s/mAERL43dMbGKi2n/download
unzip download
mv results/* .
rm download
rm -r results