for i in `cat gene_list.txt`
do
  for j in `cat gene_list.txt`
  do
    echo $i,$j
    if $i!=$j;
    then
      perl data_clear.pl $i $j

  done
done
