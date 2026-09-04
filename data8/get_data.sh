#!/bin/bash
wget https://cloud.oca.eu/index.php/s/C9GinP22frBbMco/download
unzip download
rm download

mv data8/* .
unzip HLTau.zip
rm HLTau.zip

unzip LOW.zip
rm LOW.zip

unzip MID.zip
rm MID.zip

rm -r data8