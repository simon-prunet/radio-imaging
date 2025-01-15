#!/bin/bash
wget https://cloud.oca.eu/index.php/s/zjdNWwNYA7pYNGd/download
unzip download
mv interleaved_datasets/* .
rm download
rm -r interleaved_datasets