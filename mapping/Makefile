DOCKER_RUN= docker run --rm -v $$(pwd):/tmp -u$$(id -u):$$(id -g) ff

INSERTIONS = \
	TG1.ins.gz TG2.ins.gz \
	AG1.ins.gz AG2.ins.gz AC1.ins.gz AC2.ins.gz TC1.ins.gz TC2.ins.gz \
	GT1.ins.gz GT2.ins.gz

TARGETS= $(INSERTIONS)

TG1= TG_ACG
TG2= TG_TCT
AG1= AG_AAC AG_TTC
AG2= AG_AGG AG_ACA 
AC1= AC_ATT AC_TAA
AC2= AC_TGC AC_CCG
TC1= TC_ATC TC_TGT
TC2= TC_TGA TC_ACT
GT1= GT_ACT GT_AAC
GT2= GT_TTC GT_CCG

# Note that we add a % sign here. This will generate
# a pattern rule with multiple targets, which is the
# only way to tell Make that all fasta files are generated
# by the same command.
FASTA_BATCH_1= $(addsuffix %fasta.gz, $(TG1) $(TG2))
FASTA_BATCH_2= $(addsuffix %fasta.gz, $(AG1) $(AG2) $(AC1) $(AC2) $(TC1) $(TC2))
FASTA_BATCH_3= $(addsuffix %fasta.gz, $(GT1) $(GT2))

all: $(TARGETS)

mm9_pT2_unmasked.fasta.gz:
	$(DOCKER_RUN) wget -q https://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz -O $@
	$(DOCKER_RUN) cat pT2.fasta.gz >> $@

mm9_pT2_unmasked.gem: mm9_pT2_unmasked.fasta.gz
	$(DOCKER_RUN) gunzip mm9_pT2_unmasked -c > tmp-genome.fasta
	$(DOCKER_RUN) gem-indexer -i tmp-genome.fasta -o mm9_pT2_unmasked
	$(DOCKER_RUN) rm tmp-genome.fasta

# Download 86 GB of data from the European Nucleotide Archive.
iPCR1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/049/SRR10564749/SRR10564749_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 47eb08d47327de2e46e0175a5fb39c23"
iPCR1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/049/SRR10564749/SRR10564749_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep f99a29fe459a3e9649399998faedb5f4"
iPCR2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/050/SRR10564750/SRR10564750_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 30dcc4f681f98189a6dd3138b4f85399"
iPCR2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/050/SRR10564750/SRR10564750_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8dccba2efbb71294babaa05ab791c31c"
iPCR3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/051/SRR10564751/SRR10564751_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 54a1d6b9535ea9b99dce15363b9c1855"
iPCR3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/051/SRR10564751/SRR10564751_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d09a4512788c9c9b816ce8aca079417f"

# This creates all fasta files for extra 10 GB on disk.
$(FASTA_BATCH_1): iPCR1_1.fastq.gz iPCR1_2.fastq.gz
	$(DOCKER_RUN) python scripts/makefasta_nextseq.py $^
$(FASTA_BATCH_2): iPCR2_1.fastq.gz iPCR2_2.fastq.gz
	$(DOCKER_RUN) python scripts/makefasta_hiseq.py $^ AG_AC_TC
$(FASTA_BATCH_3): iPCR3_1.fastq.gz iPCR3_2.fastq.gz
	$(DOCKER_RUN) python scripts/makefasta_hiseq.py $^ GT

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
%.ins-raw.gz: $$(addsuffix .map.gz, $$($$*)) $$(addsuffix .stc.gz, $$($$*)) mm9_pT2_unmasked.fasta.gz
	$(DOCKER_RUN) python scripts/collect.py $^ | gzip > $@

# Remove contaminants (the order matters in the calls below).
%1.ins.gz: %1.ins-raw.gz %2.ins-raw.gz
	$(DOCKER_RUN) python scripts/remove_contaminants.py $^ | gzip > $@
%2.ins.gz: %2.ins-raw.gz %1.ins-raw.gz
	$(DOCKER_RUN) python scripts/remove_contaminants.py $^ | gzip > $@

cleantmp:
	rm -f *.fasta.gz *.map.gz *.stc.gz

clean:
	rm -f *fasta.gz *.map.gz *.stc.gz *.ins-raw.gz *.ins.gz
