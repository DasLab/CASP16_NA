# write results to this file
out_f="casp16_rna_stochi.txt"
echo 'model num_chains' > $out_f


for f in R*/R0???o/*
do
echo $f
  # count number of res 1 there is
  countA=$(grep -o "P     . .   1" $f | wc -l)
  
  # write out this number of chains
  echo $(basename -- "$f" .pdb) $countA >> $out_f
  
  # if they self reported stoich, check it matches number res 1
  if grep -q "STOICH" $f
  then
    reported_count=$(sed -n -e 's/^.*STOICH[ ]*.//p' $f)
    if [ $countA -ne $reported_count ]
    then
      echo "ERROR in REPORTED COUNT" $f $countA $reported_count
    fi
  else
    echo $f "no reported count"
  fi

  # check the chain number matches count
  # turns out this is wrong a lot...
  for i in $(seq 0 $((${countA}-1)))
  do
    present=$(grep -o "P     . ${i}   1" $f | wc -l)
    # echo "PRESESNT" $present
    if [ $present -ne 1 ]
    then
      echo "ERROR in" $f "chain" $i "not present" $present
    fi
  done
done
