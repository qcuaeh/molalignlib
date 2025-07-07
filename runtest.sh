#!/bin/sh
./build/rmsd-conformer -align -remap -pipe xyz -stats -N 20 < "$1" 2>&1 >output
