#!/bin/bash

if [ "$#" -eq 0 ]; then
	echo "Usage: $0 path/to/YamlDirectory/"
	exit 1
fi

YAMLDIRPATH=$1

for yaml in `ls $YAMLDIRPATH/*.yaml`
do
	echo $yaml
	./GHphase4.py $yaml
done


