###############FUNANNOTATE#############

GENOME_ASSEMBLE=al_canu_ref_patch_HIC.tidy.masked.fasta
OUTPUT=/project/amphioxus/Annotation/funanno/al_ref
LEFT=al_EJ01.deconseq_clean.PE_filtered.PE.R1.fq.gz
RIGHT=al_EJ01.deconseq_clean.PE_filtered.PE.R2.fq.gz
LEFT2=al_EJ03.deconseq_clean.PE_filtered.PE.R1.fq.gz
RIGHT2=al_EJ03.deconseq_clean.PE_filtered.PE.R2.fq.gz
LEFT3=al_EJ05.deconseq_clean.PE_filtered.PE.R1.fq.gz
RIGHT3=al_EJ05.deconseq_clean.PE_filtered.PE.R2.fq.gz

ISO_SEQ=hq_transcripts.fastq

SPECIES="Asymmetron_lucayanum"       

#########TRAIN#############

funannotate train -i $GENOME_ASSEMBLE -o $OUTPUT \
                  -l $LEFT $LEFT2 $LEFT3 \
                  -r $RIGHT $RIGHT2 $RIGHT3 \
                  --no_trimmomatic \
                  --no_normalize_reads \
                  --memory 50G --cpus 8 \
                  --pasa_db sqlite \
                  --max_intronlen 200000 \
                  --species $SPECIES
                  --pacbio_isoseq $ISO_SEQ \
                  
#########PREDICT#############
                    
funannotate predict -i $GENOME_ASSEMBLE \
            -o $OUTPUT -s $SPECIES --cpus 8 --repeats2evm --busco_db metazoa \
            --keep_no_stops --optimize_augustus --organism other --max_intronlen 200000 \
            --protein_evidence /public/renyifan/project/amphioxus/Sequence/Protein_seq/Bl.protein.faa
#           -w augustus:10 hiq:10 pasa:10 snap:0 genemark:0 glimmerhmm:0 proteins:10 transcripts:1

##########UPDATE##############

funannotate update -i $OUTPUT --cpus 8 \
                   -l $LEFT -r $RIGHT \ 
                   --pacbio_isoseq $ISO_SEQ

############Annotate##########

emapper.py -i $WORK_DICT/update_results/${SPECIES}.proteins.fa --output ./annotate_misc/${SPECIES} -d euk -m diamond --cpu 8

interproscan.sh -i $WORK_DICT/update_results/${SPECIES}.proteins.fa -f xml -dp 2>&1 | tee interproscan.log

funannotate annotate -i $WORK_DICT --busco_db metazoa --cpus 8 \
                     --iprscan $WORK_DICT/annotate_misc/${SPECIES}.proteins.fa.xml \
                     --eggnog $WORK_DICT/annotate_misc/${SPECIES}.emapper.annotations 2>&1 | tee funannotate-annotate.log
                   
############################################