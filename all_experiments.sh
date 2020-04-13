#! /bin/bash
for f in data/*.pickle
do
    python art_eval.py $f
done
