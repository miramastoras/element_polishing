HAPPY_VCF_GZ=$1

mkdir ./counting/
cd counting

zcat $HAPPY_VCF_GZ | grep -v "^#"  | grep ":F" | sort | uniq | cut -f10 | cut -d":" -f6-7  > truth.FPFN.txt

zcat $HAPPY_VCF_GZ | grep -v "^#"  |  grep ":F" | sort | uniq | cut -f11 | cut -d":" -f6-7  > query.FPFN.txt

paste -d":" truth.FPFN.txt query.FPFN.txt | sort | uniq -c  | sed 's/ \+ /\t/g' | cut -f1 -d " "| sed -r 's/\s+//g' > counts.txt
paste -d":" truth.FPFN.txt query.FPFN.txt | sort | uniq | cut -d ":" -f1-2 > truth.sum.txt
paste -d":" truth.FPFN.txt query.FPFN.txt | sort | uniq | cut -d ":" -f3-4 > query.sum.txt

paste -d"\t" counts.txt truth.sum.txt query.sum.txt > sums.corrected.txt

# truth should be het, query is wrong het
awk '$2 ~ /het/ { print }' sums.corrected.txt | awk '$3 ~ /het/ { print }' | awk -F'\t' '{sum+=$1;} END{print sum;}' > truth_het_query_Fhet.txt

# truth should be het, query is hom or ref call
awk '$2 ~ /het/ { print }' sums.corrected.txt | awk '$3 ~ /hom/ { print }' | awk -F'\t' '{sum+=$1;} END{print sum;}' > truth_het_query_hom.txt
awk '$2 ~ /het/ { print }' sums.corrected.txt | awk '$3 ~ /no/ { print }' | awk -F'\t' '{sum+=$1;} END{print sum;}' >> truth_het_query_hom.txt
awk -F'\t' '{sum+=$1;} END{print sum;}' truth_het_query_hom.txt > tmp; mv tmp truth_het_query_hom.txt

# truth should be hom query is het
awk '$2 ~ /hom/ { print }' sums.corrected.txt | awk '$3 ~ /het/ { print }' | awk -F'\t' '{sum+=$1;} END{print sum;}' > truth_hom_query_het.txt
awk '$2 ~ /no/ { print }' sums.corrected.txt | awk '$3 ~ /het/ { print }' | awk -F'\t' '{sum+=$1;} END{print sum;}' >> truth_hom_query_het.txt
awk -F'\t' '{sum+=$1;} END{print sum;}' truth_hom_query_het.txt > tmp; mv tmp truth_hom_query_het.txt

# truth should be hom, query is wrong hom or no call
awk '$2 ~ /hom/ { print }' sums.corrected.txt | awk '$3 ~ /hom/ { print }' | awk -F'\t' '{sum+=$1;} END{print sum;}' > truth_hom_query_Fhom.txt
awk '$2 ~ /hom/ { print }' sums.corrected.txt | awk '$3 ~ /no/ { print }' | awk -F'\t' '{sum+=$1;} END{print sum;}' >> truth_hom_query_Fhom.txt
awk '$2 ~ /no/ { print }' sums.corrected.txt | awk '$3 ~ /hom/ { print }' | awk -F'\t' '{sum+=$1;} END{print sum;}' >> truth_hom_query_Fhom.txt
awk -F'\t' '{sum+=$1;} END{print sum;}' truth_hom_query_Fhom.txt > tmp ; mv tmp truth_hom_query_Fhom.txt

echo "TruthHet_Query_WrongHet,TruthHet_QueryHom,TruthHom_QueryHet,TruthHom_QueryWrongHom" > ../counts_from_happy.csv
paste -d"," truth_het_query_Fhet.txt truth_het_query_hom.txt truth_hom_query_het.txt truth_hom_query_Fhom.txt >> ../counts_from_happy.csv

rm -rf ./counting/
