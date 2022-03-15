find -maxdepth 2 -name '*nr*' | cut -d'/' -f2 | awk '{printf "igvam_dbval -d %s/nr-snp-kmer.tsv -n %s -t 2 -L <(cut -f1 ../snp_calling_hqsubsets/%s.path.list) 1> %s/kmer_profiles.tsv 2> %s/kmer_profiles.log\n", $1, $1, $1, $1, $1}' | head -n 1

cat missing_species.list | xargs -I[] -n1 -P7 bash -c "python /home/ubuntu/proj/snpMLST/snp_mlst/build_db_new.py tt extract --ref-genome ./[]/reference.fna --vcf ./[]/core_snps.vcf --msa ./[]/temp/mummer4/[]/msa.fa --out ./[]/nr --kmer-type all --coords ./[]/coords.tsv --no-reduction &> ./[]/kmer_xtract.log" &
