#!/bin/bash
../build/tb ../data/Buddha200.pgm --benchmark
for file in "Casting200.pgm" "Dancing200.pgm" "Fertility200.pgm" "Filigree200.pgm" \
            "Genus3200.pgm" "Grayloc200.pgm" "Greek200.pgm" "Hand200.pgm" \
            "Neptune_smaller200.pgm" "Pegasus200.pgm" "Twirl200.pgm"
do
  ../build/tb ../data/${file} --benchmark
done