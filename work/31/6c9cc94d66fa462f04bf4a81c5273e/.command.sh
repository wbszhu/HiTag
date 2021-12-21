#!/bin/bash -ue
cat test2.valid.pairs | grep -v "^#" | awk 'BEGIN{c=0;t=0;OFS="	"}{if($2==$4){c+=1}else{t+=1}}END{print NR,c,t}' > test2.remove_noise.count
