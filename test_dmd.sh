#!/bin/bash

if [ -f "testtemp" ]; then
    echo "Running tests would overwrite testtemp file"
    exit 1
fi

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

if [[ $( ./bin/CaVEMaN --interval testtemp,0.6 | sha1sum | awk {'print toupper($1)'}) == "19E2745363F67F5DDFF93C70295F4530DB735E6F" ]]; then
    echo "Passed: interval."
else
    echo "Failed: interval."
    exit 1
fi

rm -f testtemp

if [[ $( ./bin/CaVEMaN --correct data/eQTL --bed data/phenotype.bed --vcf data/genotype.vcf.gz | sha1sum | awk {'print toupper($1)'}) == "FDED43C25211773C54FB6F854FFED8D0A0CFEA9C" ]]; then
    echo "Passed: correct phenotypes."
else
    echo "Failed: correct phenotypes."
    exit 1
fi

if [[ $( ./bin/CaVEMaN --normal --correct data/eQTL --bed data/phenotype.bed --vcf data/genotype.vcf.gz | sha1sum | awk {'print toupper($1)'}) == "BA58F5A4E604A5185270E074CC9BC754DD582C7E" ]]; then
    echo "Passed: correct with normalisation."
else
    echo "Failed: correct with normalisation."
    exit 1
fi

if [[ $( ./bin/CaVEMaN --correct data/eQTL --bed data/phenotype.bed --vcf data/genotype.vcf.gz --cov data/covariates | sha1sum | awk {'print toupper($1)'}) == "798E1AD6FF67BCEE6C96B19E8297F107439E6609" ]]; then
    echo "Passed: correct with covariates."
else
    echo "Failed: correct with covariates."
    exit 1
fi

echo "All tests completed successfully."
