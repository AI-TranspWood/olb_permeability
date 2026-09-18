#!/bin/bash

# Run a calculation on image.vti and compare the output in permeability.dat with baseline.dat
BASEDIR=$(dirname "$0")

rm -fr permeability.dat

mpirun -np 1 PermPorousWood tests/image.vti ImageFile 1e-5 1e-3 40 1e5 0.8 1e-16 1e-6 1000 1e-5 2 0

if [ ! -f permeability.dat ]; then
    echo "Error: permeability.dat not found. The calculation may have failed."
    exit 1
fi

baseline=$(cat $BASEDIR/baseline.dat)
test=$(cat permeability.dat)
threshold=6e-2

diff=$(awk -v b="$baseline" -v t="$test" 'BEGIN { print (b > t) ? b - t : t - b }')
diff_relative=$(awk -v d="$diff" -v b="$baseline" 'BEGIN { print -d / b }')

echo "Baseline: $baseline"
echo "Test: $test"
echo "Difference: $diff"
echo "Relative Difference: $diff_relative"
if awk -v d="$diff_relative" -v t="$threshold" 'BEGIN { exit (d <= t) ? 0 : 1 }'; then
    echo "Test passed: permeability.dat is within the threshold of baseline.dat"
    exit 0
else
    echo "Test failed: permeability.dat differs from baseline.dat by $diff, which is above the threshold of $threshold"
    exit 1
fi

