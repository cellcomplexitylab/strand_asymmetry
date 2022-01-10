DOCKER_RUN= docker run --rm -v $$(pwd):/tmp -u$$(id -u):$$(id -g) ff

GSE93238_gene.fpkm.txt.gz:
	$(DOCKER_RUN) wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93238/suppl/GSE93238_gene.fpkm.txt.gz -O $@

Mus_musculus.NCBIM37.67.gtf.gz:
	$(DOCKER_RUN) wget -q http://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz -O $@

%.stc: %.txt
	$(DOCKER_RUN) starcode -t4 -d2 $< > $@

%.lst: %.stc
	$(DOCKER_RUN) seeq -ed2 TGCAACGAATTCATTAGTGCG $< | awk '{if ($$2 > 1) print}' | grep -o "^[GATC]\{32\}" > $@
