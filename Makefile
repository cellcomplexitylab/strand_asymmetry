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

barcode_view_all_events_with_mapping.txt: $(INSERTIONS) black.lst.gz
	$(DOCKER_RUN) /bin/bash scripts/extract_all_events_with_mapping.sh > $@

barcode_view_CRISPR_with_mapping.txt: mapping/GT1.ins.gz mapping/GT2.ins.gz black.lst.gz
	$(DOCKER_RUN) /bin/bash scripts/extract_CRISPR_events_with_mapping.sh > $@

barcode_view_all_events_without_mapping.txt: black.lst.gz
	$(DOCKER_RUN) /bin/bash scripts/extract_all_events_without_mapping.sh > $@

barcode_view_CRISPR_without_mapping.txt: black.lst.gz
	$(DOCKER_RUN) /bin/bash scripts/extract_CRISPR_events_without_mapping.sh > $@

barcode_view_GT_with_genic.txt: barcode_view_all_events_with_mapping.txt misc/GSE93238_gene.fpkm.txt.gz misc/Mus_musculus.NCBIM37.67.gtf.gz
	$(DOCKER_RUN) R -f scripts/mapped_GT_events_to_genic.R

barcode_view_CRISPR_GT_with_genic.txt: barcode_view_CRISPR_with_mapping.txt misc/GSE93238_gene.fpkm.txt.gz misc/Mus_musculus.NCBIM37.67.gtf.gz
	$(DOCKER_RUN) R -f scripts/mapped_CRISPR_GT_events_to_genic.R

barcode_view_GT_with_histones.txt: barcode_view_all_events_with_mapping.txt misc/CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt.gz
	$(DOCKER_RUN) R -f scripts/mapped_GT_events_to_histones.R

barcode_view_CRISPR_GT_with_histones.txt: barcode_view_CRISPR_with_mapping.txt misc/CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt.gz
	$(DOCKER_RUN) R -f scripts/mapped_CRISPR_GT_events_to_histones.R

barcode_view_all_features_with_mapping.txt: barcode_view_all_events_with_mapping.txt misc/CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt.gz
	$(DOCKER_RUN) R -f scripts/make_feature_table.R

bias_by_experiments_with_mapping.txt: barcode_view_all_events_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_bias_by_experiments_with_mapping.R

bias_by_experiments_without_mapping.txt: barcode_view_all_events_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_bias_by_experiments_without_mapping.R

bias_by_CRISPR_without_mapping.txt: barcode_view_CRISPR_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_bias_by_CRISPR_without_mapping.R

bias_for_GT_by_genic.txt: barcode_view_GT_with_genic.txt
	$(DOCKER_RUN) R -f scripts/get_bias_for_GT_by_genic.R

bias_for_GT_by_histones.txt: barcode_view_GT_with_histones.txt
	$(DOCKER_RUN) R -f scripts/get_bias_for_GT_by_histones.R

conflicts_by_experiments_with_mapping.txt: barcode_view_all_events_with_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_conflicts_by_experiments_with_mapping.R

conflicts_by_experiments_without_mapping.txt: barcode_view_all_events_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/get_conflicts_by_experiments_without_mapping.R

mutual_information.txt: barcode_view_all_events_without_mapping.txt
	$(DOCKER_RUN) python scripts/compute_mutual_info.py $< > $@

ctrl_GT_barcodes_with_errors.txt: # $(MISMATCHES)
	$(DOCKER_RUN) /bin/bash scripts/extract_ctrl_GT_barcodes_with_errors.sh > $@

all_ctrl_GT_barcodes.txt: # $(MISMATCHES)
	$(DOCKER_RUN) /bin/bash scripts/extract_all_ctrl_GT_barcodes.sh > $@

features_ctrl_GT_barcodes.txt: ctrl_GT_barcodes_with_errors.txt all_ctrl_GT_barcodes.txt
	$(DOCKER_RUN) R -f scripts/map_ctrl_GT_barcodes.R

mutations.txt: # $(MISMATCHES)
	$(DOCKER_RUN) /bin/bash scripts/extract_mutations.sh > $@

#learning_curves_full.txt: barcode_view_all_features_with_mapping.txt
#	$(DOCKER_RUN) python scripts/compute_mutual_info.py $< > $@

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

figures/mutual_information.pdf: mutual_information.txt
	$(DOCKER_RUN) R -f scripts/plot_mutual_information.R

# Figure 5.
figures/bias_GT_genic.pdf: bias_for_GT_by_genic.txt
	$(DOCKER_RUN) R -f scripts/plot_bias_GT_genic.R

figures/bias_GT_histones.pdf: bias_for_GT_by_histones.txt
	$(DOCKER_RUN) R -f scripts/plot_bias_GT_histones.R

figures/learning_curves.pdf: learning_curves_full.txt learning_curves_null.txt
	$(DOCKER_RUN) R -f scripts/plot_learning_curves.R

figures/learning_curves_GC.pdf: learning_curves_GC_full.txt learning_curves_GC_null.txt
	$(DOCKER_RUN) R -f scripts/plot_learning_curves_GC.R

figures/learning_curves_conflict.pdf: learning_curves_24h-PCR_full.txt learning_curves_24h-PCR_null.txt
	$(DOCKER_RUN) R -f scripts/plot_learning_curves_conflict.R

# Figure 6.
figures/oligo_mutations.pdf:
	$(DOCKER_RUN) R -f scripts/plot_oligo_mutations.R

figures/C_to_T_transitions.pdf:
	$(DOCKER_RUN) R -f scripts/plot_C_to_T_transitions.R

# Figure 7.
figures/CRISPR_circos.pdf: misc/gRNA_counts_GT1.txt misc/gRNA_counts_GT2.txt
	$(DOCKER_RUN) R -f scripts/plot_circos_CRISPR.R

figures/bias_CRISPR.pdf: bias_by_CRISPR_without_mapping.txt
	$(DOCKER_RUN) R -f scripts/plot_bias_CRISPR.R

# Complementary data sets
misc/GSE93238_gene.fpkm.txt.gz:
	cd misc && $(MAKE) GSE93238_gene.fpkm.txt.gz

misc/Mus_musculus.NCBIM37.67.gtf.gz:
	cd misc && $(MAKE) Mus_musculus.NCBIM37.67.gtf.gz

misc/CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt.gz:
	wget http://epistemnet.bioinfo.cnio.es/download/others/CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt.gz -O $@

# !!!! This part still has the original FASTQ identifers !!!!
misc/gRNA_counts_GT1.txt:
	$(DOCKER_RUN) python scripts/gather_gRNA_counts.py mismatches/15D*.co.gz mismatches/17F*.co.gz > $@

misc/gRNA_counts_GT2.txt:
	$(DOCKER_RUN) python scripts/gather_gRNA_counts.py mismatches/16E*.co.gz mismatches/18G*.co.gz > $@
