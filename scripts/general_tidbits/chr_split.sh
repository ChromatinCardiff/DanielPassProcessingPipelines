#!/bin/bash

grep -w 'Chr1' ES16.sam >Chr1_grep.txt &
grep -w 'Chr2' ES16.sam >Chr2_grep.txt &
grep -w 'Chr3' ES16.sam >Chr3_grep.txt &
grep -w 'Chr4' ES16.sam >Chr4_grep.txt &
grep -w 'Chr5' ES16.sam >Chr5_grep.txt &
grep -w 'mitochondria' ES16.sam >mito_grep.txt &
grep -w 'chloroplast' ES16.sam >chloro_grep.txt &
