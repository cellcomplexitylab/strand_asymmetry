DOCKER_RUN= docker run --rm -v $$(pwd):/tmp -u$$(id -u):$$(id -g) ff

all:
	echo "Hi!"
#	+$(MAKE) -C mapping
#	+$(MAKE) -C mismatch
#	+$(MAKE) -C figures

figures/integ_circos.pdf:
	$(DOCKER_RUN) R -f plot_circos_figure.R

figures/spie_genes.pdf: misc/GSE93238_gene.fpkm.txt.gz misc/Mus_musculus.NCBIM37.67.gtf.gz
	$(DOCKER_RUN) R -f plot_spie_genes.R

misc/GSE93238_gene.fpkm.txt.gz:
	cd misc && $(MAKE) GSE93238_gene.fpkm.txt.gz

misc/Mus_musculus.NCBIM37.67.gtf.gz:
	cd misc && $(MAKE) Mus_musculus.NCBIM37.67.gtf.gz
