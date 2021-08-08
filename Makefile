DOCKER_RUN= docker run --rm -v $$(pwd):/tmp -u$$(id -u):$$(id -g) ff

INSERTIONS= $(addprefix mapping/, TG1.ins.gz TG2.ins.gz \
	AG1.ins.gz AG2.ins.gz AC1.ins.gz AC2.ins.gz \
	TC1.ins.gz TC2.ins.gz GT1.ins.gz GT2.ins.gz)

MISMATCHES= $(addprefix mismatches., insertmismatches)

all:
	echo "$(INSERTIONS)"
#	+$(MAKE) -C mapping
#	+$(MAKE) -C mismatch
#	+$(MAKE) -C figures



# Base analyses
$(INSERTIONS):
	cd mapping && $(MAKE)

black.lst.gz: $(INSERTIONS) # $(MISMATCHES)
	$(DOCKER_RUN) /bin/bash scripts/make_blacklist.sh | gzip > $@

barcode_view_all_events_with_mapping.txt: black.lst.gz
	$(DOCKER_RUN) /bin/bash scripts/extract_all_events_with_mapping.sh > $@

barcode_view_all_events_without_mapping.txt: black.lst.gz
	$(DOCKER_RUN) /bin/bash scripts/extract_all_events_without_mapping.sh > $@

bias_by_experiments_with_mapping.txt: barcode_view_all_events_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_bias_by_experiments_with_mapping.R

bias_by_experiments_without_mapping.txt: barcode_view_all_events_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_bias_by_experiments_without_mapping.R

conflicts_by_experiments_with_mapping.txt: barcode_view_all_events_with_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_conflicts_by_experiments_with_mapping.R

conflicts_by_experiments_without_mapping.txt: barcode_view_all_events_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_conflicts_by_experiments_without_mapping.R


# Figure 2.
figures/integ_circos.pdf: $(INSERTIONS)
	$(DOCKER_RUN) R -f scripts/plot_circos_integrations.R

figures/insertion_rates.pdf: $(INSERTIONS)
	$(DOCKER_RUN) R -f scripts/plot_insertion_rates.R

figures/spie_genes.pdf: $(INSERTIONS) misc/GSE93238_gene.fpkm.txt.gz misc/Mus_musculus.NCBIM37.67.gtf.gz
	$(DOCKER_RUN) R -f scripts/plot_spie_genes.R

figures/insertions_in_repeats.pdf: $(INSERTIONS)
	$(DOCKER_RUN) R -f scripts/plot_insertions_in_repeats.R

# Figure 3.
figures/bias.pdf: bias_by_experiments_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/plot_bias.R

# Figure 4.
figures/conflicts.pdf: conflicts_by_experiments_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/plot_conflicts.R

# Figure 6.
figures/CRISPR_circos.pdf: misc/gRNA_counts_GT1.txt misc/gRNA_counts_GT2.txt
	$(DOCKER_RUN) R -f scripts/plot_circos_CRISPR.R

# Complementary data sets
misc/GSE93238_gene.fpkm.txt.gz:
	cd misc && $(MAKE) GSE93238_gene.fpkm.txt.gz

misc/Mus_musculus.NCBIM37.67.gtf.gz:
	cd misc && $(MAKE) Mus_musculus.NCBIM37.67.gtf.gz

# !!!! This part still has the original FASTQ identifers !!!!
misc/gRNA_counts_GT1.txt:
	$(DOCKER_RUN) python scripts/gather_gRNA_counts.py mismatches/15D*.co.gz mismatches/17F*.co.gz > $@

misc/gRNA_counts_GT2.txt:
	$(DOCKER_RUN) python scripts/gather_gRNA_counts.py mismatches/16E*.co.gz mismatches/18G*.co.gz > $@
