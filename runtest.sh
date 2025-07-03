#!/bin/sh
./build/rmsd-conformer -align -remap -pipe xyz -stats -N 8 < "$1" 2>&1 >output
