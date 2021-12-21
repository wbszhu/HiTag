#!/bin/bash -ue
zcat test2.pairs.gz | grep -v "^#" | awk 'BEGIN{c=0;t=0;OFS="	"}{if($2==$4){c+=1}else{t+=1}}END{print NR,c,t}' > test2.bam2pairs.count
