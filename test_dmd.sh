#!/bin/bash

for x in "testtemp" "testtemp1" "testtemp2"; do
    if [ -f $x ]; then
	echo "Running tests would overwrite a temporary file:" $x
	exit 1
    fi
done

./bin/CaVEMaN --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 100000,4 > testtemp

if [[ $( sha1sum testtemp | awk {'print toupper($1)'}) == "FA9DDAD12C9DC98C447AB0660D574C2E061DF869" ]]; then
    echo "Passed: CaVEMaN."
else
    echo "Failed: CaVEMaN."
    exit 1
fi

if [[ $( ./bin/CaVEMaN --best testtemp | sha1sum | awk {'print toupper($1)'}) == "C530F1AC9E5D57C43F68FA64C7861780C2F742EE" ]]; then
    echo "Passed: extract best."
else
    echo "Failed: extract best."
    exit 1
fi

rm -f testtemp

if [[ $( ./bin/CaVEMaN --single-signal --eqtl data/eQTL --bed data/phenotype.bed --vcf data/genotype.vcf.gz | sha1sum | awk {'print toupper($1)'}) == "FDED43C25211773C54FB6F854FFED8D0A0CFEA9C" ]]; then
    echo "Passed: correct phenotypes."
else
    echo "Failed: correct phenotypes."
    exit 1
fi

if [[ $( ./bin/CaVEMaN --normal --single-signal --eqtl data/eQTL --bed data/phenotype.bed --vcf data/genotype.vcf.gz | sha1sum | awk {'print toupper($1)'}) == "BA58F5A4E604A5185270E074CC9BC754DD582C7E" ]]; then
    echo "Passed: correct with normalisation."
else
    echo "Failed: correct with normalisation."
    exit 1
fi

if [[ $( ./bin/CaVEMaN --single-signal --eqtl data/eQTL --bed data/phenotype.bed --vcf data/genotype.vcf.gz --cov data/covariates | sha1sum | awk {'print toupper($1)'}) == "798E1AD6FF67BCEE6C96B19E8297F107439E6609" ]]; then
    echo "Passed: correct with covariates."
else
    echo "Failed: correct with covariates."
    exit 1
fi

./bin/CaVEMaN --simulate --vcf data/genotype.vcf.gz --bed data/phenotype.bed --eqtl data/eQTL --perm 4 --out testtemp
./bin/CaVEMaN ./bin/CaVEMaN --bed testtemp --vcf data/genotype.vcf.gz --perm 10000,4 --job-number 1 --genes 10 --out testtemp2
./bin/CaVEMaN --get-weights --results testtemp2 --rank testtemp --weights testtemp1

if [[ ($( sha1sum testtemp | awk {'print toupper($1)'}) == "6FF7D9DFD6DD9BD4880BB8642B2EFB4A4C0878D9") &&
	  ($( sha1sum testtemp1 | awk {'print toupper($1)'}) == "26B8704D251F7B6B789399A5D8DB05174DB1E49E")]]; then
    echo "Passed: estimating ranks and weights."
else
    echo "Failed: estimating ranks and weights."
    exit 1
fi

rm -f testtemp*

echo "All tests completed successfully."
