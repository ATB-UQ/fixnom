v="test_run"
rm log_$v.txt
echo "$v"
rm -r $v
mkdir $v
cd input
for i in *.lib; do
echo "$i" >> ../log_$v.txt
echo "$i"
../../fixnom.pl -i $i -o out_$i >> ../log_$v.txt
mv out_$i ../$v/
rm translate*
diff ../$v/out_$i ../correct/out_$i >> ../log_$v.txt
diff -q ../$v/out_$i ../correct/out_$i 
done
cd ..
