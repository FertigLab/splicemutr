# combining data splicemutr data into a single file

sed -n 1p data_splicemutr1.txt > data_splicemutr.txt
awk 'FNR>1' data_splicemutr*.txt >> data_splicemutr.txt
