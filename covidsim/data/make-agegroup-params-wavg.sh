#!/bin/bash

awk -F ',' '
function arr2str(a,  tmpstr,i){
  tmpstr="";
  for (i=1; i<=length(a); ++i)
    tmpstr=tmpstr "\t" a[i];
  return(substr(tmpstr,2));
}

/^#/{next;}

headerline==0 {
  headerline=$0;
  split($0,headerarr,",");
  print "agegroup\t" arr2str(headerarr);
  next;
}

{ agegroup[$1]=$0; }

function calc_weighted_avg(a,b,wa,wb){
  wt=wa+wb;
  wa=wa/wt;
  wb=wb/wt;
  for (i=2; i<length(a); ++i)
    a[i]=a[i]*wa+b[i]*wb;
}

END{
  for (i in agegroup) {
    if (i==80) break;
    split(agegroup[i],a,",");
    a[1]=i-5;
    print (i-10)/5 "\t" arr2str(a);
    a[1]=i;
    print (i-5)/5 "\t" arr2str(a);
  }
  split(agegroup[80],a,",");
  split(agegroup[100],b,",");
  calc_weighted_avg(a,b,3388.488 + 2442.147, 1736.567 + 1077.555 + 490.577 + 130.083 + 15.834);
  a[1]=80-5;
  print (80-10)/5 "\t" arr2str(a);
  a[1]=80;
  print (80-5)/5 "\t" arr2str(a);
}' agegroup.infparams.orig > agegroup.infparams

