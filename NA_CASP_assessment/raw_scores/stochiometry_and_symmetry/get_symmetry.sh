sym_out='casp16_sym_round1_rna-homomultimers.csv'
echo "model,symmetries,rmsds" > $sym_out

for f in */R1???o/* */R1???v?o/*
do
  f2=$(basename $f)_simple
  out=$(basename $f)_sym
  # copy just one atom per res (P)
  cat $f | grep "         P" > $f2
  sed -i "s/ P     / CA    /g" $f2
  sed -i "s/        P/        C/g" $f2
  sed -i "s/   C / ALA /g" $f2
  sed -i "s/   G / ALA /g" $f2
  sed -i "s/   A / ALA /g" $f2
  sed -i "s/   U / ALA /g" $f2
  ~/ananas $f2 -C 50 > $out
  syms=$(sed -n -e 's/^.*Symmetry group : //p' $out)
  rmsds=$(sed -n -e 's/^.*Average RMSD : //p' $out)
  echo ${f2%.*}"," $syms "," $rmsds >> $sym_out
# collect in pretty large rmsd and in post update to be more reasonable
done

