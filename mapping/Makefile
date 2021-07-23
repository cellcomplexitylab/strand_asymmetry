DOCKER_RUN= docker run --rm -v $$(pwd):/tmp -u$$(id -u):$$(id -g) ff

INSERTIONS = \
	CT1.ins CA1.ins GT1.ins GA1.ins \
	CT2.ins CA2.ins GT2.ins GA2.ins \
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
FASTA= $(addsuffix %fasta, $(CA1) $(GT1) $(GA1) $(CT1) \
       $(CA2) $(GT2) $(GA2) $(CT2) $(TC1) $(TC2))

iPCR_FASTQ= iPCR1_read1.fastq.gz iPCR1_read2.fastq.gz \
	    iPCR2_read1.fastq.gz iPCR2_read2.fastq.gz \
	    iPCR3_read1.fastq.gz iPCR3_read2.fastq.gz

all: $(TARGETS)

mm9_pT2_unmasked.fa:
	wget https://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz
	zcat mm9.fa.gz pT2.fa.gz > mm9_pT2_unmasked.fa

mm9_pT2_unmasked.gem: mm9_pT2_unmasked.fa
	$(DOCKER_RUN) gem-indexer -i $< -o mm9_pT2_unmasked


# Download 86 GB of data from the European Nucleotide Archive.
iPCR1_read1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/049/SRR10564749/SRR10564749_1.fastq.gz -O iPCR1_read1.fastq.gz && \
	md5sum iPCR1_read1.fastq.gz | grep 47eb08d47327de2e46e0175a5fb39c23
iPCR1_read2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/049/SRR10564749/SRR10564749_2.fastq.gz -O iPCR1_read2.fastq.gz && \
	md5sum iPCR1_read2.fastq.gz | grep f99a29fe459a3e9649399998faedb5f4
iPCR2_read1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/050/SRR10564750/SRR10564750_1.fastq.gz -O iPCR2_read1.fastq.gz && \
	md5sum iPCR2_read1.fastq.gz | grep 30dcc4f681f98189a6dd3138b4f85399
iPCR2_read2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/050/SRR10564750/SRR10564750_2.fastq.gz -O iPCR2_read2.fastq.gz && \
	md5sum iPCR2_read2.fastq.gz | grep 8dccba2efbb71294babaa05ab791c31c
iPCR3_read1.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/051/SRR10564751/SRR10564751_1.fastq.gz -O iPCR3_read1.fastq.gz && \
	md5sum iPCR3_read1.fastq.gz | grep 54a1d6b9535ea9b99dce15363b9c1855
iPCR3_read2.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/051/SRR10564751/SRR10564751_2.fastq.gz -O iPCR3_read2.fastq.gz && \
	md5sum iPCR3_read2.fastq.gz | grep d09a4512788c9c9b816ce8aca079417f

# This creates all fasta files (see definition of 'FASTA').
$(FASTA): $(iPCR_FASTQ)
	$(DOCKER_RUN) python makefasta_nextseq.py iPCR1_read1.fastq.gz iPCR1_read2.fastq.gz
	$(DOCKER_RUN) python makefasta_hiseq.py iPCR2_read1.fastq.gz iPCR2_read2.fastq.gz CA_GT_GA
	$(DOCKER_RUN) python makefasta_hiseq.py iPCR3_read1.fastq.gz iPCR3_read2.fastq.gz TC

#%.map: %.fasta
#	gem-mapper -I /data/mm10_pT2_unmasked.gem \
#		-q ignore --unique-mapping -T4 -i $< > $@
%.map.gz: %.fasta mm9_pT2_unmasked.gem
	$(DOCKER_RUN) gem-mapper -I mm9_pT2_unmasked.gem \
		-q ignore --unique-mapping -T4 -i $< | gzip > $@


.SECONDARY:
%.stc: %.map.gz
	$(DOCKER_RUN) zcat $< | cut -f1 | starcode -t4 -d2 --print-clusters -o $@

# The secondary expansions will allow to replace the matching
# stem of the rule by the variables declared above.
.SECONDEXPANSION:
%.ins: $$(addsuffix .map.gz, $$($$*)) $$(addsuffix .stc, $$($$*))
	$(DOCKER_RUN) python collect.py $^ > $@

black.lst: $(INSERTIONS)
	$(DOCKER_RUN) python make_blacklist.py $^ > $@

cleantmp:
	rm -f *.fa *.map.gz *.stc

clean:
	rm -f *fasta *.map.gz *.stc *.ins *.bcd
