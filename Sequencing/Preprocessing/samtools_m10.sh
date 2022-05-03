#!/bin/bash

echo `samtools view -q 10 -h -o ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_default_m10 ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_default`

echo `samtools view -q 10 -h -o ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_local_default_m10 ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_local`

echo `samtools view -q 10 -h -o ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_mp_10_2_m10 ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_mp_10_2`

echo `samtools view -q 10 -h -o ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_Npenalty_m10 ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_Npenalty`

echo `samtools view -q 10 -h -o ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_rdg_10_6_m10 ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_rdg_10_6`

echo `samtools view -q 10 -h -o ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_rfg_10_6_m10 ~/sequence_data/hg38_aligned/hg38_soft_fullchrom_rfg_10_6`
