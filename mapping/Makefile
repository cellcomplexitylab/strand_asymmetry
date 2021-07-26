DOCKER_RUN= docker run --rm -v $$(pwd):/tmp -u$$(id -u):$$(id -g) ff

INSERTIONS = \
	CT1.ins CT2.ins \
	CA1.ins CA2.ins GT1.ins GT2.ins GA1.ins GA2.ins \
	TC1.ins TC2.ins

TARGETS= $(INSERTIONS) black.lst

CT1= CT_ACG
CT2= CT_TCT
CA1= CA_AAC CA_TTC
CA2= CA_AGG CA_ACA 
GT1= GT_ATC GT_TGT
GT2= GT_TGA GT_ACT
GA1= GA_ATT GA_TAA
GA2= GA_TGC GA_CCG
TC1= TC_ACT TC_AAC
TC2= TC_TTC TC_CCG

# Note that we add a % sign here. This will generate
# a pattern rule with multiple targets, which is the
# only way to tell Make that all fasta files are generated
# by the same command.
FASTA_BATCH_1= $(addsuffix %fasta.gz, $(CT1) $(CT2))
FASTA_BATCH_2= $(addsuffix %fasta.gz, $(CA1) $(CA2) $(GT1) $(GT2) $(GA1) $(GA2))
FASTA_BATCH_3= $(addsuffix %fasta.gz, $(TC1) $(TC2))

all: $(TARGETS)

mm9_pT2_unmasked.fasta.gz:
	$(DOCKER_RUN) wget -q https://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz -O $@
	$(DOCKER_RUN) cat pT2.fa.gz >> $@

mm9_pT2_unmasked.gem: mm9_pT2_unmasked.fasta.gz
	$(DOCKER_RUN) gunzip mm9_pT2_unmasked -c > tmp-genome.fasta
	$(DOCKER_RUN) gem-indexer -i tmp-genome.fasta -o mm9_pT2_unmasked
	$(DOCKER_RUN) rm tmp-genome.fasta

# Download 86 GB of data from the European Nucleotide Archive.
iPCR1_read1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/049/SRR10564749/SRR10564749_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 47eb08d47327de2e46e0175a5fb39c23"
iPCR1_read2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/049/SRR10564749/SRR10564749_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep f99a29fe459a3e9649399998faedb5f4"
iPCR2_read1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/050/SRR10564750/SRR10564750_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 30dcc4f681f98189a6dd3138b4f85399"
iPCR2_read2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/050/SRR10564750/SRR10564750_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8dccba2efbb71294babaa05ab791c31c"
iPCR3_read1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/051/SRR10564751/SRR10564751_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 54a1d6b9535ea9b99dce15363b9c1855"
iPCR3_read2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/051/SRR10564751/SRR10564751_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d09a4512788c9c9b816ce8aca079417f"

# This creates all fasta files for extra 10 GB on disk.
$(FASTA_BATCH_1): iPCR1_1.fastq.gz iPCR1_2.fastq.gz
	$(DOCKER_RUN) python scripts/makefasta_nextseq.py $^
$(FASTA_BATCH_2): iPCR2_1.fastq.gz iPCR2_2.fastq.gz
	$(DOCKER_RUN) python scripts/makefasta_hiseq.py $^ CA_GT_GA
$(FASTA_BATCH_3): iPCR3_1.fastq.gz iPCR3_2.fastq.gz
	$(DOCKER_RUN) python scripts/makefasta_hiseq.py $^ TC

# Map individual files.
%.map.gz: %.fasta.gz mm9_pT2_unmasked.gem
	$(DOCKER_RUN) gunzip $< -c > tmp-$*.fasta
	$(DOCKER_RUN) /bin/bash -c "gem-mapper -I mm9_pT2_unmasked.gem -q ignore --unique-mapping -T4 -i tmp-$*.fasta | gzip > $@"
	$(DOCKER_RUN) rm tmp-$*.fasta

.SECONDARY:
%.stc.gz: %.map.gz
	$(DOCKER_RUN) /bin/bash -c "zcat $< | cut -f1 | starcode -t4 -d2 --print-clusters | gzip > $@"

# The secondary expansions will allow to replace the matching
# stem of the rule by the variables declared above.
.SECONDEXPANSION:
%.ins: $$(addsuffix .map.gz, $$($$*)) $$(addsuffix .stc.gz, $$($$*))
	$(DOCKER_RUN) python scripts/collect.py $^ > $@

black.lst: $(INSERTIONS)
	$(DOCKER_RUN) python scripts/make_blacklist.py $^ > $@

cleantmp:
	rm -f *.fa *.map.gz *.stc.gz

clean:
	rm -f *fasta.gz *.map.gz *.stc *.ins *.bcd
