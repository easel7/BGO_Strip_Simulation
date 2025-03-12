#!/bin/bash
exe=/Users/xiongzheng/software/B4/B4c/build/exampleB4c
Tag=Deuteron
# Loop through energy values from 10 to 100 GeV with a step of 10
for energy in $(seq 20 10 90)
do
    ${exe} -m ${Tag}_${energy}GeV.mac
done

# for energy in $(seq 100 100 1000)
# do
#     ${exe} -m ${Tag}_${energy}GeV.mac
# done