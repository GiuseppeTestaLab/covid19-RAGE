
LANG=C

## go to csv files directory----
cd 'data/DEGS.csv'


## setup DataMatrix file----
tail -n +2 $(ls -1 *.csv | head -1) | cut -d"," -f 1 | sort -k1,1 > ../DEGS.DataMatrix.csv


## join all columns of all files into a DataMatrix----
for i in $(ls *.csv) 
do 
join -t',' -a1 -a2 -o auto -e 'NA' ../DEGS.DataMatrix.csv <(tail -n +2 $i | sort -k1,1) >  ../DEGS.DataMatrix.tmp
mv ../DEGS.DataMatrix.tmp ../DEGS.DataMatrix.csv 
done


# setup headers----
echo " " > ../DEGS.DataMatrix.headers.txt

for i in $(ls *.csv) 
do 
i2=$(echo $i | sed -e 's/.csv//g') 
i3=$(head -1 $i | sed 's/\"//g' | sed -e 's/\,logFC/logFC/g' | sed 's/ /_/g') 
echo $i3 | awk -v i2b=$i2 'NR==1{print i2b"__"$1","i2b"__"$2","i2b"__"$3","i2b"__"$4","i2b"__"$5}' FS="," OFS="," >> ../DEGS.DataMatrix.headers.txt 
done


# add headers to the DataMatrix----
paste -sd, ../DEGS.DataMatrix.headers.txt > ../DEGS.DataMatrix.headers.b.txt

cat ../DEGS.DataMatrix.headers.b.txt ../DEGS.DataMatrix.csv | gzip -c - > ../DEGS.DataMatrix.csv.gz


# clean up temporary files----
rm ../DEGS.DataMatrix.headers.b.txt ../DEGS.DataMatrix.headers.txt  ../DEGS.DataMatrix.csv 


