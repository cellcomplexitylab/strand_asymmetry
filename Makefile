DOCKER_RUN= docker run --rm -v $$(pwd):/tmp -u$$(id -u):$$(id -g) ff

INSERTIONS= $(addprefix mapping/, AG1.ins.gz AG2.ins.gz \
	TG1.ins.gz TG2.ins.gz AC1.ins.gz AC2.ins.gz \
	TC1.ins.gz TC2.ins.gz CA1.ins.gz CA2.ins.gz)

MISMATCHES= $(addprefix mismatches., insertmismatches)

all:
	echo "$(INSERTIONS)"
#	+$(MAKE) -C mapping
#	+$(MAKE) -C mismatch
#	+$(MAKE) -C figures

# Base analyses
$(INSERTIONS):
	cd mapping && $(MAKE)

barcode_view_all_events_with_mapping.txt: $(INSERTIONS) # $(MISMATCHES)
	$(DOCKER_RUN) /bin/bash scripts/extract_all_events_with_mapping.sh > $@

barcode_view_all_events_without_mapping.txt: $(INSERTIONS) # $(MISMATCHES)
	$(DOCKER_RUN) /bin/bash scripts/extract_all_events_without_mapping.sh > $@

stats_by_experiments_with_mapping.txt: barcode_view_all_events_with_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_stats_by_experiments_with_mapping.R

stats_by_experiments_without_mapping.txt: barcode_view_all_events_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_stats_by_experiments_without_mapping.R


# Figure 2.
figures/integ_circos.pdf: $(INSERTIONS)
	$(DOCKER_RUN) R -f scripts/plot_circos_integrations.R

figures/spie_genes.pdf: $(INSERTIONS) misc/GSE93238_gene.fpkm.txt.gz misc/Mus_musculus.NCBIM37.67.gtf.gz
	$(DOCKER_RUN) R -f scripts/plot_spie_genes.R

# Figure 6.
figures/CRISPR_circos.pdf: misc/gRNA_counts_CA1.txt misc/gRNA_counts_CA2.txt
	$(DOCKER_RUN) R -f scripts/plot_circos_CRISPR.R

# Complementary data sets
misc/GSE93238_gene.fpkm.txt.gz:
	cd misc && $(MAKE) GSE93238_gene.fpkm.txt.gz

misc/Mus_musculus.NCBIM37.67.gtf.gz:
	cd misc && $(MAKE) Mus_musculus.NCBIM37.67.gtf.gz

# !!!! This part still has the original FASTQ identifers !!!!
misc/gRNA_counts_CA1.txt:
	$(DOCKER_RUN) python misc/gather_gRNA_counts.py mismatches/15D*co.gz mismatches/17F*.co.gz > $@
misc/gRNA_counts_CA2.txt:
	$(DOCKER_RUN) python misc/gather_gRNA_counts.py mismatches/16E*co.gz mismatches/18G*.co.gz > $@
