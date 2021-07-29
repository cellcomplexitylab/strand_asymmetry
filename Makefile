DOCKER_RUN= docker run --rm -v $$(pwd):/tmp -u$$(id -u):$$(id -g) ff

INSERTIONS:= $(addprefix mapping/, AG1.ins AG2.ins TG1.ins TG2.ins \
	AC1.ins AC2.ins TC1.ins TC2.ins CA1.ins CA2.ins)

all:
	echo "$(INSERTIONS)"
#	+$(MAKE) -C mapping
#	+$(MAKE) -C mismatch
#	+$(MAKE) -C figures

# Figure 2.
figures/integ_circos.pdf: $(INSERTIONS)
	$(DOCKER_RUN) R -f plot_circos_integrations.R

figures/spie_genes.pdf: $(INSERTIONS) misc/GSE93238_gene.fpkm.txt.gz misc/Mus_musculus.NCBIM37.67.gtf.gz
	$(DOCKER_RUN) R -f plot_spie_genes.R

# Figure 6.
figures/CRISPR_circos.pdf: misc/gRNA_counts_CA1.txt misc/gRNA_counts_CA2.txt
	$(DOCKER_RUN) R -f plot_circos_CRISPR.R

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

# Base analyses
$(INSERTIONS):
	cd mapping && $(MAKE)
