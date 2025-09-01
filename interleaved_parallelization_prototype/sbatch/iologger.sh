#!/bin/bash

rm cyg_a_io.dat
pidstat -d 1 -U $USER -C python -H >> cyg_a_io.dat &
processID=$!
sleep 5
python ip_l1.py configs/cyg_a_test.config
sleep 5
kill $processID
python clean_io.py cyg_a_io.dat cyg_a_io_cleaned.dat
#python ioplotter.py cyg_a_io_cleaned.dat cyg_a_io.png

rm hltau_io.dat
pidstat -d 1 -U $USER -C python -H >> hltau_io.dat &
processID=$!
sleep 5
python ip_l1.py configs/hltau_test.config
sleep 5
kill $processID
python clean_io.py hltau_io.dat hltau_io_cleaned.dat
#python ioplotter.py hltau_io_cleaned.dat hltau_io.png