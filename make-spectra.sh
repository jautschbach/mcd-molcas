for i in 0 1 2 3; do
  ~/code/mcd-c-molcas/plot-spectrum < mcdspectrum-$i
  ~/code/mcd-c-molcas/graph-it.sh -r 0.001 "ucl6-mcd-$i"
  rm graph.dat impulses.dat
done
