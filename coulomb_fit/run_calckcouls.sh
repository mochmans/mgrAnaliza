#!/bin/bash

rvalues=(3.67 3.69 3.80 3.82 3.84 3.85 3.87 3.89 3.91 3.93 3.95 5.0)

for i in ${rvalues[@]};
do
    ./calckcoul 3 $i
done

