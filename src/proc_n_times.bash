#! /bin/bash

export PYTHONPATH=/data/legs/rpete/flight/analysis_functions

python=/usr/bin/python
pdir=./data/pickle/scale_lengths
pdfdir=./plots/scale_lengths

src=mkn421
caldb=qe_N0014_qeu_N0013
b=0
e=2

#rm -f times.out

#for n in 2 3 4 5 6 7 8 9 10 11 12
for n in 5 6 7 8 9 10 11 12
do
  rm data/pickle/scale_lengths/$caldb/mkn421_sl_B${b}_E${e}*
  read min max <<<$( python3 minmax.py $src lengths $caldb $b $e )
  t=$(sh -c "time mpirun --oversubscribe -n $n \
    $python src/scale_lengths.py \
    s l \
    ./rdb/acis/${src}_fits_B${b}.rdb \
    ./rdb/$caldb/${src}_fits_B${b}.rdb \
    -i 0.2 \
    -e $e \
    --minscale $min --maxscale $max \
    --pdir $pdir/$caldb \
    --pdfdir $pdfdir/$caldb \
    --models" 2>&1 | tail -2 | head -1)

  echo $n: $t >> results/proc_n_times.out

done
