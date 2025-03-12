#!/bin/bash
source /opt/homebrew/opt/root/bin/thisroot.sh
g++ -o Normalized Normalized.cc `root-config --cflags --libs` 