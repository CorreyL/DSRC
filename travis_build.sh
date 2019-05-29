#!/bin/bash
set -e

# for the moment compile only the binary
make bin

# run the basic tests
python test/test_dsrc_bin.py test/files
