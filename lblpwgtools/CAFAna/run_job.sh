#!/bin/bash

stateFiles=("/pnfs/dune/persistent/users/cnguyen/cafana_example/State_")

systSets=("nosyst")

for i in "${!stateFiles[@]}"; do
    echo $i
    stateFile=${stateFiles[$i]}

    for systSet in "${systSets[@]}"; do
        echo $systSet
        stateFileName="${stateFiles##*/}"

        # Loop through integers 0 to 36
        for intVal in {9..9}; do
            # Dynamically create the output file name based on stateFile, systSet, and intVal
            outputFile="${stateFileName%.*}_${systSet}_624ktmwyr_test_v2.root"
            echo "Running: cpv_joint $stateFile $outputFile $systSet $intVal"
            #cpv_joint "$stateFile" "$outputFile" "$systSet" "$intVal"
        done
    done
done

