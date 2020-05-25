#!/bin/bash

awk -F '\t' 'BEGIN {f=0;} /^#/ {next;} f==0 { OFS="\t"; print $0,"Prop_noncritical_fatal"; f=1; next; } { OFS="\t"; print $0,$4-$3*$7*$5*$6; }' agegroup.infparams.ferguson.original > agegroup.infparams.ferguson
