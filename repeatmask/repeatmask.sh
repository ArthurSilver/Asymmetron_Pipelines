###############EDTA#################
EDTA.pl --genome al_canu_ref_patch_HIC.tidy.masked.rescaffolding.fasta \
        --cds Asymmetron_lucayanum.cds-transcripts.fa \
        --species others --step all --sensitive 1 --anno 1 --threads 8 2>&1 | tee Al.edta.log
###################################

###########RepeatMasker############
RepeatMasker -pa 8 -xsmall -gff -lib EDTA.TElib.fa -dir . al_canu_ref_patch_HIC.tidy.masked.rescaffolding.fasta
