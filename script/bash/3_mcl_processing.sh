#!/usr/bin/bash

awk ' BEGIN{ print "geneIDA\tgeneIDB\tbitScore"} $3 >=30 && ($8-$7+1)/$14 >=0.5 && ($10-$9+1)/$16 >= 0.5 {print $13"\t"$15"\t"$12}' filtered_data.txt >  mcl_input_30_50.txt;
awk ' BEGIN{ print "geneIDA\tgeneIDB\tbitScore"} $3 >=30 && ($8-$7+1)/$14 >=0.8 && ($10-$9+1)/$16 >= 0.8 {print $13"\t"$15"\t"$12}' filtered_data.txt >  mcl_input_30_80.txt;
awk ' BEGIN{ print "geneIDA\tgeneIDB\tbitScore"} $3 >=50 && ($8-$7+1)/$14 >=0.5 && ($10-$9+1)/$16 >= 0.5 {print $13"\t"$15"\t"$12}' filtered_data.txt >  mcl_input_50_50.txt;
awk ' BEGIN{ print "geneIDA\tgeneIDB\tbitScore"} $3 >=50 && ($8-$7+1)/$14 >=0.8 && ($10-$9+1)/$16 >= 0.8 {print $13"\t"$15"\t"$12}' filtered_data.txt >  mcl_input_50_80.txt;

mcl mcl_input_30_50.txt --abc -o mcl_output_30_50;
mcl mcl_input_30_80.txt --abc -o mcl_output_30_80;
mcl mcl_input_50_50.txt --abc -o mcl_output_50_50;
mcl mcl_input_50_80.txt --abc -o mcl_output_50_80;

awk -F'\t' '{size = gsub(/\t/, ""); if (size > max_size) max_size = size} END {print "Total Families:", NR, "Max Family Size:", max_size}' mcl_output_30_50;   
awk -F'\t' '{size = gsub(/\t/, ""); if (size > max_size) max_size = size} END {print "Total Families:", NR, "Max Family Size:", max_size}' mcl_output_30_80;
awk -F'\t' '{size = gsub(/\t/, ""); if (size > max_size) max_size = size} END {print "Total Families:", NR, "Max Family Size:", max_size}' mcl_output_50_50;
awk -F'\t' '{size = gsub(/\t/, ""); if (size > max_size) max_size = size} END {print "Total Families:", NR, "Max Family Size:", max_size}' mcl_output_50_80;


