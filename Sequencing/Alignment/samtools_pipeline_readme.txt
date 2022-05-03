##samtools pipeline

#pairend
samtools view -Sb [input dir] | samtools view -b -q 10 -f 2 - | samtools sort - [output dir]; 

#singleend
samtools view -Sb [input dir] | samtools view -b -q 10 - | samtools sort - [output dir]; 