DOCKER_RUN= docker run --rm -v $$(pwd):/tmp -u$$(id -u):$$(id -g) ff

all:
	echo "Hi!"
#	+$(MAKE) -C mapping
#	+$(MAKE) -C mismatch
#	+$(MAKE) -C figures

figures/integ_circos.pdf:
	$(DOCKER_RUN) R -f do_circos_figure.R
