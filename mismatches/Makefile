DOCKER_RUN= docker run --rm -v $$(pwd):/tmp -u$$(id -u):$$(id -g) ff

TARGETS=                                \
	AC_6xPCR_24h_rep1.co.gz         \
	AC_6xPCR_24h_rep2.co.gz         \
	AC_6xPCR_24h_rep3.co.gz         \
	AC_6xPCR_24h_rep4.co.gz         \
	AC_6xPCR_48h_rep1.co.gz         \
	AC_6xPCR_48h_rep2.co.gz         \
	AC_6xPCR_48h_rep3.co.gz         \
	AC_6xPCR_48h_rep4.co.gz         \
	AC_6xPCR_no_ISceI_rep1.co.gz    \
	AC_6xPCR_no_ISceI_rep2.co.gz    \
	AC_6xPCR_no_ISceI_rep3.co.gz    \
	AC_6xPCR_no_ISceI_rep4.co.gz    \
	AC_6xPCR_no_ISceI_rep5.co.gz    \
	AC_6xPCR_no_ISceI_rep6.co.gz    \
	AC_6xPCR_no_ISceI_rep7.co.gz    \
	AC_6xPCR_no_ISceI_rep8.co.gz    \
	AC_6xPCR_no_ISceI_rep9.co.gz    \
	AC_LA_24h_rep1.co.gz            \
	AC_LA_24h_rep2.co.gz            \
	AC_LA_24h_rep3.co.gz            \
	AC_LA_24h_rep4.co.gz            \
	AC_LA_48h_rep1.co.gz            \
	AC_LA_48h_rep2.co.gz            \
	AC_LA_48h_rep3.co.gz            \
	AC_LA_48h_rep4.co.gz            \
	AC_LA_no_ISceI_rep1.co.gz       \
	AC_LA_no_ISceI_rep2.co.gz       \
	AC_LA_no_ISceI_rep3.co.gz       \
	AC_LA_no_ISceI_rep4.co.gz       \
	AG_6xPCR_24h_rep1.co.gz         \
	AG_6xPCR_24h_rep2.co.gz         \
	AG_6xPCR_24h_rep3.co.gz         \
	AG_6xPCR_24h_rep4.co.gz         \
	AG_6xPCR_48h_rep1.co.gz         \
	AG_6xPCR_48h_rep2.co.gz         \
	AG_6xPCR_no_ISceI_rep1.co.gz    \
	AG_6xPCR_no_ISceI_rep2.co.gz    \
	AG_6xPCR_no_ISceI_rep3.co.gz    \
	AG_6xPCR_no_ISceI_rep4.co.gz    \
	AG_6xPCR_no_ISceI_rep5.co.gz    \
	AG_6xPCR_no_ISceI_rep6.co.gz    \
	AG_6xPCR_no_ISceI_rep7.co.gz    \
	AG_LA_24h_rep1.co.gz            \
	AG_LA_24h_rep2.co.gz            \
	AG_LA_24h_rep3.co.gz            \
	AG_LA_24h_rep4.co.gz            \
	AG_LA_48h_rep1.co.gz            \
	AG_LA_48h_rep2.co.gz            \
	AG_LA_48h_rep3.co.gz            \
	AG_LA_48h_rep4.co.gz            \
	AG_LA_no_ISceI_rep1.co.gz       \
	AG_LA_no_ISceI_rep2.co.gz       \
	AG_LA_no_ISceI_rep3.co.gz       \
	AG_LA_no_ISceI_rep4.co.gz       \
	AG_LA_no_ISceI_rep5.co.gz       \
	AG_LA_no_ISceI_rep6.co.gz       \
	AG_LA_no_ISceI_rep7.co.gz       \
	AG_LA_no_ISceI_rep8.co.gz       \
	GT_6xPCR_24h_rep1.co.gz         \
	GT_6xPCR_24h_rep2.co.gz         \
	GT_6xPCR_24h_rep3.co.gz         \
	GT_6xPCR_24h_rep4.co.gz         \
	GT_6xPCR_48h_rep1.co.gz         \
	GT_6xPCR_48h_rep2.co.gz         \
	GT_6xPCR_48h_rep3.co.gz         \
	GT_6xPCR_48h_rep4.co.gz         \
	GT_6xPCR_no_ISceI_rep1.co.gz    \
	GT_6xPCR_no_ISceI_rep2.co.gz    \
	GT_6xPCR_no_ISceI_rep3.co.gz    \
	GT_6xPCR_no_ISceI_rep4.co.gz    \
	GT_LA_24h_rep1.co.gz            \
	GT_LA_24h_rep2.co.gz            \
	GT_LA_24h_rep3.co.gz            \
	GT_LA_24h_rep4.co.gz            \
	GT_LA_48h_rep1.co.gz            \
	GT_LA_48h_rep2.co.gz            \
	GT_LA_48h_rep3.co.gz            \
	GT_LA_48h_rep4.co.gz            \
	GT_LA_no_ISceI_rep1.co.gz       \
	GT_LA_no_ISceI_rep2.co.gz       \
	GT_LA_no_ISceI_rep3.co.gz       \
	GT_LA_no_ISceI_rep4.co.gz       \
	TC_6xPCR_24h_rep1.co.gz         \
	TC_6xPCR_24h_rep2.co.gz         \
	TC_6xPCR_24h_rep3.co.gz         \
	TC_6xPCR_24h_rep4.co.gz         \
	TC_6xPCR_48h_rep1.co.gz         \
	TC_6xPCR_48h_rep2.co.gz         \
	TC_6xPCR_48h_rep3.co.gz         \
	TC_6xPCR_48h_rep4.co.gz         \
	TC_6xPCR_no_ISceI_rep1.co.gz    \
	TC_6xPCR_no_ISceI_rep2.co.gz    \
	TC_6xPCR_no_ISceI_rep3.co.gz    \
	TC_6xPCR_no_ISceI_rep4.co.gz    \
	TC_6xPCR_no_ISceI_rep5.co.gz    \
	TC_6xPCR_no_ISceI_rep6.co.gz    \
	TC_6xPCR_no_ISceI_rep7.co.gz    \
	TC_6xPCR_no_ISceI_rep8.co.gz    \
	TC_LA_24h_rep1.co.gz            \
	TC_LA_24h_rep3.co.gz            \
	TC_LA_24h_rep4.co.gz            \
	TC_LA_48h_rep1.co.gz            \
	TC_LA_48h_rep2.co.gz            \
	TC_LA_48h_rep3.co.gz            \
	TC_LA_48h_rep4.co.gz            \
	TC_LA_no_ISceI_rep1.co.gz       \
	TC_LA_no_ISceI_rep2.co.gz       \
	TC_LA_no_ISceI_rep3.co.gz       \
	TC_LA_no_ISceI_rep4.co.gz       \
	TC_LA_no_ISceI_rep5.co.gz       \
	TC_LA_no_ISceI_rep6.co.gz       \
	TC_LA_no_ISceI_rep7.co.gz       \
	TC_LA_no_ISceI_rep8.co.gz       \
	TG_6xPCR_24h_rep1.co.gz         \
	TG_6xPCR_24h_rep2.co.gz         \
	TG_6xPCR_24h_rep3.co.gz         \
	TG_6xPCR_24h_rep4.co.gz         \
	TG_6xPCR_48h_rep1.co.gz         \
	TG_6xPCR_48h_rep2.co.gz         \
	TG_6xPCR_48h_rep3.co.gz         \
	TG_6xPCR_48h_rep4.co.gz         \
	TG_6xPCR_no_ISceI_rep1.co.gz    \
	TG_6xPCR_no_ISceI_rep2.co.gz    \
	TG_6xPCR_no_ISceI_rep3.co.gz    \
	TG_6xPCR_no_ISceI_rep4.co.gz    \
	TG_LA_24h_rep1.co.gz            \
	TG_LA_24h_rep2.co.gz            \
	TG_LA_24h_rep3.co.gz            \
	TG_LA_24h_rep4.co.gz            \
	TG_LA_48h_rep1.co.gz            \
	TG_LA_48h_rep2.co.gz            \
	TG_LA_48h_rep3.co.gz            \
	TG_LA_48h_rep4.co.gz            \
	TG_LA_no_ISceI_rep1.co.gz       \
	TG_LA_no_ISceI_rep2.co.gz       \
	TG_LA_no_ISceI_rep3.co.gz       \
	TG_LA_no_ISceI_rep4.co.gz


all: $(TARGETS)

# Download 193 GB of sequencing data to disk.
AC_6xPCR_24h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/023/SRR10564723/SRR10564723_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2e766183d21ed3e7dd1362a1183c223b"
AC_6xPCR_24h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/023/SRR10564723/SRR10564723_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 614689f061b496182e7899e85f2e0ae5"
AC_6xPCR_24h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/024/SRR10564724/SRR10564724_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 086831b880becd8424a4bec9c93b0909"
AC_6xPCR_24h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/024/SRR10564724/SRR10564724_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep cf1814cf18cb60c268b7bb87c4af5fb7"
AC_6xPCR_24h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/025/SRR10564725/SRR10564725_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d8b30b8e7d0a7113563b952312f0c12e"
AC_6xPCR_24h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/025/SRR10564725/SRR10564725_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 07636149bf87af5ab2c210f9742aa93c"
AC_6xPCR_24h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/026/SRR10564726/SRR10564726_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep a2c91ba7cc8cfd019562b5c04265521c"
AC_6xPCR_24h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/026/SRR10564726/SRR10564726_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep eab41a49d552e9944951d65383bbcd4a"
AC_6xPCR_48h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/027/SRR10564727/SRR10564727_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 1898091b7549253b8355d87d21f36d2b"
AC_6xPCR_48h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/027/SRR10564727/SRR10564727_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3be2169881e061e12f2a9d7c96129397"
AC_6xPCR_48h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/028/SRR10564728/SRR10564728_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep fdf73c9e56f33ea715cd9a279fb449ba"
AC_6xPCR_48h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/028/SRR10564728/SRR10564728_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3b6d1f96f1ba8ba2b08cc9d7e7f11537"
AC_6xPCR_48h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/029/SRR10564729/SRR10564729_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep f2239a44d6aa3d5414c5874f43ccad70"
AC_6xPCR_48h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/029/SRR10564729/SRR10564729_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8f32e1897cc40793d84faa8564d05f45"
AC_6xPCR_48h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/008/SRR10564808/SRR10564808_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep f2a755df2f9dfb18d9751801f4cce78b"
AC_6xPCR_48h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/008/SRR10564808/SRR10564808_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c1dc4071a6e1f84f3f64f18992bd775c"
AC_6xPCR_no_ISceI_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/009/SRR10564809/SRR10564809_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7bf03e52ba534b22e4798a310745ee2a"
AC_6xPCR_no_ISceI_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/009/SRR10564809/SRR10564809_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 948e67b783af4db21f27cf74a3344c1c"
AC_6xPCR_no_ISceI_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/010/SRR10564810/SRR10564810_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 08a60c61779b9253f2f7dfbe1df78740"
AC_6xPCR_no_ISceI_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/010/SRR10564810/SRR10564810_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 46fb9bb1347ca1e424b08b661d345a1d"
AC_6xPCR_no_ISceI_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/011/SRR10564811/SRR10564811_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 668d4f1b2a724a0a50fdecc40d106e60"
AC_6xPCR_no_ISceI_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/011/SRR10564811/SRR10564811_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 19b35e559c4dcffa675902dd66fb9560"
AC_6xPCR_no_ISceI_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/012/SRR10564812/SRR10564812_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep a2bb8fa392a72af9f6a3b10b3c288b44"
AC_6xPCR_no_ISceI_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/012/SRR10564812/SRR10564812_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep b9e86dc7f0cb2b88f870e1e937d267ad"
AC_6xPCR_no_ISceI_rep5_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/013/SRR10564813/SRR10564813_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 9dda5a896d6d3e5a18de84f6b132d9fe"
AC_6xPCR_no_ISceI_rep5_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/013/SRR10564813/SRR10564813_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep e086f9a9961e65f8e3e81b1e79028fef"
AC_6xPCR_no_ISceI_rep6_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/014/SRR10564814/SRR10564814_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 29a364b1afdcccb28223ddb3cae297c6"
AC_6xPCR_no_ISceI_rep6_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/014/SRR10564814/SRR10564814_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3d5e1cb303dc1448e908bfb161a5c549"
AC_6xPCR_no_ISceI_rep7_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/015/SRR10564815/SRR10564815_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c8e2dcc436f35da5e4f77d75f3b36366"
AC_6xPCR_no_ISceI_rep7_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/015/SRR10564815/SRR10564815_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 839727c056a3726d9de918336c49e722"
AC_6xPCR_no_ISceI_rep8_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/016/SRR10564816/SRR10564816_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 25d7c5f4c85e9983de936ab4e1ddda28"
AC_6xPCR_no_ISceI_rep8_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/016/SRR10564816/SRR10564816_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 50d9869016e3d36a8362b2594011d160"
AC_6xPCR_no_ISceI_rep9_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/017/SRR10564817/SRR10564817_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2cdcc3c6b134d079058de3868fd594f9"
AC_6xPCR_no_ISceI_rep9_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/017/SRR10564817/SRR10564817_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 938e87bfe8666a7b2b17104ca9ecbfd2"
AC_LA_24h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/018/SRR10564818/SRR10564818_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3a00fcffa0cada1cdaa320bbd111f088"
AC_LA_24h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/018/SRR10564818/SRR10564818_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 04d7ea458226588701314f2b36073628"
AC_LA_24h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/019/SRR10564819/SRR10564819_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8546e12c20fc5ff14d4ffcc0c009dac7"
AC_LA_24h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/019/SRR10564819/SRR10564819_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 88ca4725f968f3fac155e717c3834f8a"
AC_LA_24h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/020/SRR10564820/SRR10564820_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 81e9faf6cdd588ffd6808f24407f46dd"
AC_LA_24h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/020/SRR10564820/SRR10564820_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep ec444be0dd478216ee061e87e8485d27"
AC_LA_24h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/021/SRR10564821/SRR10564821_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2473f3f5e8349a1f73971c571bb4f291"
AC_LA_24h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/021/SRR10564821/SRR10564821_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 083053d98ccf487704636440f6ba2230"
AC_LA_48h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/022/SRR10564822/SRR10564822_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep f35ed0c6201e6e6653e87b834f116bd1"
AC_LA_48h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/022/SRR10564822/SRR10564822_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep e724dd9e98e6a1d4cd4584121b507895"
AC_LA_48h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/023/SRR10564823/SRR10564823_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 6d28a95f6401a81a38f4fa1715290c96"
AC_LA_48h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/023/SRR10564823/SRR10564823_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 230eaf3fde2bc5c4c64c62f72d857455"
AC_LA_48h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/024/SRR10564824/SRR10564824_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 996d5a2401eb6a152268a6cbac83b2d3"
AC_LA_48h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/024/SRR10564824/SRR10564824_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 12c954cd9c49fdcb8acc416136d83e4c"
AC_LA_48h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/025/SRR10564825/SRR10564825_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 16d278f7986194501dfcf48b65d99601"
AC_LA_48h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/025/SRR10564825/SRR10564825_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2b220014ff33ce28bfd44ea947fc1a78"
AC_LA_no_ISceI_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/026/SRR10564826/SRR10564826_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c6d920996f093b3748878dd32524b053"
AC_LA_no_ISceI_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/026/SRR10564826/SRR10564826_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 09ee98f6409490b2ed02cb9b8e1b70da"
AC_LA_no_ISceI_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/027/SRR10564827/SRR10564827_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 528a77d35f88417145b1ea6b48136f90"
AC_LA_no_ISceI_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/027/SRR10564827/SRR10564827_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 83bff6f069ce9297603ac132b46048d0"
AC_LA_no_ISceI_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/028/SRR10564828/SRR10564828_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 9f81d4ee26e5a5d1b2999cd2194c1314"
AC_LA_no_ISceI_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/028/SRR10564828/SRR10564828_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2b63dee642decae3226b93a0940c102f"
AC_LA_no_ISceI_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/029/SRR10564829/SRR10564829_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 75d0c9f2503fc2149f4d2806c4accf77"
AC_LA_no_ISceI_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/029/SRR10564829/SRR10564829_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8584d0d95a62c7bc1d5e9566a35cc694"
AG_6xPCR_24h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/030/SRR10564830/SRR10564830_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2dc500dde2b0cd5ea3b13155a637e626"
AG_6xPCR_24h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/030/SRR10564830/SRR10564830_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 46b1b25c6229e04123c65abbfe5e5ebb"
AG_6xPCR_24h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/031/SRR10564831/SRR10564831_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep f321c1135e4f82fc5f306545b6b9992d"
AG_6xPCR_24h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/031/SRR10564831/SRR10564831_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 9bc0bc941a2de7fc5830e4ed6b12e03e"
AG_6xPCR_24h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/032/SRR10564832/SRR10564832_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d08b395df3e55957994b2a037f520d50"
AG_6xPCR_24h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/032/SRR10564832/SRR10564832_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 1192b448b163a95ee8681cd8ed21fcf1"
AG_6xPCR_24h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/033/SRR10564833/SRR10564833_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 57f6a0bcedd61df791aa17504ae190bb"
AG_6xPCR_24h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/033/SRR10564833/SRR10564833_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 719a4b99d6af4459609cba0d17fa0f86"
AG_6xPCR_48h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/034/SRR10564834/SRR10564834_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep dff9ec6cce74104290caf639adf6ff75"
AG_6xPCR_48h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/034/SRR10564834/SRR10564834_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2ec66777847111d0d53df4c4bdeeff51"
AG_6xPCR_48h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/035/SRR10564835/SRR10564835_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c3aefc28597648fc906d588c93d4a837"
AG_6xPCR_48h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/035/SRR10564835/SRR10564835_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 82685b20bedfd2b42f8694d1ce9d153b"
AG_6xPCR_no_ISceI_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/036/SRR10564836/SRR10564836_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep e4a2d19b62a2811f03d8c75a504f199e"
AG_6xPCR_no_ISceI_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/036/SRR10564836/SRR10564836_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep fd12ee2678284e789fb340084331b453"
AG_6xPCR_no_ISceI_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/037/SRR10564837/SRR10564837_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 71ef9b9b653f6561fd6b0d5d90060512"
AG_6xPCR_no_ISceI_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/037/SRR10564837/SRR10564837_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep a634e7bcab79b97590522d7690e58eef"
AG_6xPCR_no_ISceI_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/038/SRR10564838/SRR10564838_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 96c7407da1872d5add643db7213258bb"
AG_6xPCR_no_ISceI_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/038/SRR10564838/SRR10564838_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7f1b757f9b184992e925fbf68370426a"
AG_6xPCR_no_ISceI_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/039/SRR10564839/SRR10564839_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 0940b3b083bc66cdefecadfb576b5680"
AG_6xPCR_no_ISceI_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/039/SRR10564839/SRR10564839_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep af54742ade358218e38f3c6c8a3bb026"
AG_6xPCR_no_ISceI_rep5_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/040/SRR10564840/SRR10564840_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 9df3421798eb35b8443a21dca490ec92"
AG_6xPCR_no_ISceI_rep5_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/040/SRR10564840/SRR10564840_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 999a3f3cb33282669db7add0b84dba76"
AG_6xPCR_no_ISceI_rep6_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/041/SRR10564841/SRR10564841_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep feac17283b07ed2a265a498d33396fe8"
AG_6xPCR_no_ISceI_rep6_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/041/SRR10564841/SRR10564841_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep ff708c00fa08fd3477a2b6b1f0009358"
AG_6xPCR_no_ISceI_rep7_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/042/SRR10564842/SRR10564842_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 1ecfa467efaf0da1e0681326dafe7dcd"
AG_6xPCR_no_ISceI_rep7_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/042/SRR10564842/SRR10564842_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 70fb36aa0df82687ac8d95157f2e70d3"
AG_LA_24h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/043/SRR10564843/SRR10564843_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 919e80128acfb29cbf89c27d060b29f3"
AG_LA_24h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/043/SRR10564843/SRR10564843_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 00f35d09a3478781c108c0fd6b088edf"
AG_LA_24h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/044/SRR10564844/SRR10564844_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 9d12ff784842433beaae3bf47f0e7d11"
AG_LA_24h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/044/SRR10564844/SRR10564844_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 28df67049353ca81b4d89d587903283f"
AG_LA_24h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/045/SRR10564845/SRR10564845_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep a785a5166cb626dc364798504572d2e1"
AG_LA_24h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/045/SRR10564845/SRR10564845_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 155917251e2d6bcdd34655cd091d2230"
AG_LA_24h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/046/SRR10564846/SRR10564846_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 69375eff57cd9cb9212449caa9d7e66c"
AG_LA_24h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/046/SRR10564846/SRR10564846_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8fad40d5f8a2c057daf0671a6f469b9c"
AG_LA_48h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/047/SRR10564847/SRR10564847_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 9ee2c501594f6c7bd266747344e76f6b"
AG_LA_48h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/047/SRR10564847/SRR10564847_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 635430ec34526bd55f8a03255902aa4b"
AG_LA_48h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/048/SRR10564848/SRR10564848_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 5ceea12736acb76d3cb4a7dacf02bd84"
AG_LA_48h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/048/SRR10564848/SRR10564848_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep bd2ab76a63230b1c5ff0e2154bfe31c5"
AG_LA_48h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/049/SRR10564849/SRR10564849_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c0b7080db4a8b6841447fcc9e3263d75"
AG_LA_48h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/049/SRR10564849/SRR10564849_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d39da88cbd280ca24fec50bdb7618eee"
AG_LA_48h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/050/SRR10564850/SRR10564850_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep fda4b4b9e9177ddb3e866853fbf0b6c7"
AG_LA_48h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/050/SRR10564850/SRR10564850_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 4b055551e0fc193363852cfc7cc05809"
AG_LA_no_ISceI_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/051/SRR10564851/SRR10564851_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7c6427a3f137cf6f6aba7ac114d58637"
AG_LA_no_ISceI_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/051/SRR10564851/SRR10564851_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c4c23c863db5f6b89a89f5bc6f0b8ff7"
AG_LA_no_ISceI_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/052/SRR10564852/SRR10564852_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c11aa70f59c0c453cc4c71d7e39ba3db"
AG_LA_no_ISceI_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/052/SRR10564852/SRR10564852_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 37bdeb48d28025fec6eb99a7e025c004"
AG_LA_no_ISceI_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/053/SRR10564853/SRR10564853_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 40502a5030ec90d1bf2b2dbb44c26925"
AG_LA_no_ISceI_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/053/SRR10564853/SRR10564853_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c386df611d4b2008ed46f2ac35acf5b5"
AG_LA_no_ISceI_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/054/SRR10564854/SRR10564854_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep eeb465553f27d61148c464e743dec247"
AG_LA_no_ISceI_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/054/SRR10564854/SRR10564854_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 44facd2cb843772f3fdace7377e4f77b"
AG_LA_no_ISceI_rep5_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/055/SRR10564855/SRR10564855_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 6cee2f56897bd30fe72fbe44b66d80e1"
AG_LA_no_ISceI_rep5_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/055/SRR10564855/SRR10564855_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 11573f6ca0500195ba552b837512cc41"
AG_LA_no_ISceI_rep6_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/056/SRR10564856/SRR10564856_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2ef7ce6017bc9155daa8a13d4f154f6a"
AG_LA_no_ISceI_rep6_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/056/SRR10564856/SRR10564856_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7e21c6a5538bd72ce4e9affb5505921a"
AG_LA_no_ISceI_rep7_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/057/SRR10564857/SRR10564857_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep a9c3a7acbf54d45632c00da228c2f681"
AG_LA_no_ISceI_rep7_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/057/SRR10564857/SRR10564857_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c4969d9cb217692301f58a380732de28"
AG_LA_no_ISceI_rep8_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/058/SRR10564858/SRR10564858_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep fde20ba7d59142a16d4adb47eaabbc89"
AG_LA_no_ISceI_rep8_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/058/SRR10564858/SRR10564858_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 5cb677cb7a9f6a703b40a3a8aeaa9d06"
#CA_6xPCR_24h_rep1_1.fastq.gz:
GT_6xPCR_24h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/059/SRR10564859/SRR10564859_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c1bda4a351705a54c9968faf5114319f"
#CA_6xPCR_24h_rep1_2.fastq.gz:
GT_6xPCR_24h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/059/SRR10564859/SRR10564859_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 1a5af16e1ebb41c78fc6cffc13b45cb7"
#CA_6xPCR_24h_rep2_1.fastq.gz:
GT_6xPCR_24h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/060/SRR10564860/SRR10564860_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8f6f07ef9142765550ff00faac4d03a9"
#CA_6xPCR_24h_rep2_2.fastq.gz:
GT_6xPCR_24h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/060/SRR10564860/SRR10564860_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 4ccf6ad8b18f387aab60d32a240c776a"
#CA_6xPCR_24h_rep3_1.fastq.gz:
GT_6xPCR_24h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/061/SRR10564861/SRR10564861_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3247098818815c6f3df27fe843f18e13"
#CA_6xPCR_24h_rep3_2.fastq.gz:
GT_6xPCR_24h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/061/SRR10564861/SRR10564861_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep caf3c327dbe3f9a942e0b1b47912987b"
#CA_6xPCR_24h_rep4_1.fastq.gz:
GT_6xPCR_24h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/062/SRR10564862/SRR10564862_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7029c89254d23c9f0dc842f13234ebc8"
#CA_6xPCR_24h_rep4_2.fastq.gz:
GT_6xPCR_24h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/062/SRR10564862/SRR10564862_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 27a368937712bbb3d51c1f553bc42983"
#CA_6xPCR_48h_rep1_1.fastq.gz:
GT_6xPCR_48h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/063/SRR10564863/SRR10564863_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7903085fdde2b24315519fbb8415f268"
#CA_6xPCR_48h_rep1_2.fastq.gz:
GT_6xPCR_48h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/063/SRR10564863/SRR10564863_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c9e6798548fe646919f6ddf90826ce22"
#CA_6xPCR_48h_rep2_1.fastq.gz:
GT_6xPCR_48h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/052/SRR10564752/SRR10564752_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 9877129566aa12720cd9afb922bf3cc6"
#CA_6xPCR_48h_rep2_2.fastq.gz:
GT_6xPCR_48h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/052/SRR10564752/SRR10564752_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 0ffd9ac4e34552a90662d43700da7e64"
#CA_6xPCR_48h_rep3_1.fastq.gz:
GT_6xPCR_48h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/053/SRR10564753/SRR10564753_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 366bc322a2dc083cede40dcb0c171627"
#CA_6xPCR_48h_rep3_2.fastq.gz:
GT_6xPCR_48h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/053/SRR10564753/SRR10564753_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d6abf73d8c93d8c77d71fe466a731878"
#CA_6xPCR_48h_rep4_1.fastq.gz:
GT_6xPCR_48h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/054/SRR10564754/SRR10564754_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep a925dd40c121a84826f953abf8fdaff5"
#CA_6xPCR_48h_rep4_2.fastq.gz:
GT_6xPCR_48h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/054/SRR10564754/SRR10564754_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep b95ba2bdb8dd3ff00d3821f1108c0852"
#CA_6xPCR_no_ISceI_rep1_1.fastq.gz:
GT_6xPCR_no_ISceI_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/055/SRR10564755/SRR10564755_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 9c904a0ba2bf752ccb63e6475943ae45"
#CA_6xPCR_no_ISceI_rep1_2.fastq.gz:
GT_6xPCR_no_ISceI_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/055/SRR10564755/SRR10564755_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep e7a8698d74ab6c2496d3a22b0b164eab"
#CA_6xPCR_no_ISceI_rep2_1.fastq.gz:
GT_6xPCR_no_ISceI_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/056/SRR10564756/SRR10564756_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c332a68a88c8b9a0805cf4ff3d3cc1a9"
#CA_6xPCR_no_ISceI_rep2_2.fastq.gz:
GT_6xPCR_no_ISceI_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/056/SRR10564756/SRR10564756_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 0dfd254a0a45b1dfcf494b53f8a42d45"
#CA_6xPCR_no_ISceI_rep3_1.fastq.gz:
GT_6xPCR_no_ISceI_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/057/SRR10564757/SRR10564757_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3760d077bdf6324dcfbc2dadfe2f7795"
#CA_6xPCR_no_ISceI_rep3_2.fastq.gz:
GT_6xPCR_no_ISceI_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/057/SRR10564757/SRR10564757_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 86f8a4573d2b3fb3a13cd9d2d6ca8da1"
#CA_6xPCR_no_ISceI_rep4_1.fastq.gz:
GT_6xPCR_no_ISceI_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/058/SRR10564758/SRR10564758_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 0e7a03607121bb76d7967dd3ed9ee22d"
#CA_6xPCR_no_ISceI_rep4_2.fastq.gz:
GT_6xPCR_no_ISceI_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/058/SRR10564758/SRR10564758_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8f521ad2af9789a56417dad68504d7af"
#CA_LA_24h_rep1_1.fastq.gz:
GT_LA_24h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/059/SRR10564759/SRR10564759_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep cfff042b077d03ecabcc2e1961b9b1ab"
#CA_LA_24h_rep1_2.fastq.gz:
GT_LA_24h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/059/SRR10564759/SRR10564759_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep a448a34a16174215a490a43ce28d0845"
#CA_LA_24h_rep2_1.fastq.gz:
GT_LA_24h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/060/SRR10564760/SRR10564760_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep cb78b3576661eb38ff7f0d91ded76465"
#CA_LA_24h_rep2_2.fastq.gz:
GT_LA_24h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/060/SRR10564760/SRR10564760_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep bef5da94a32a19d3b459e2a92722e547"
#CA_LA_24h_rep3_1.fastq.gz:
GT_LA_24h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/061/SRR10564761/SRR10564761_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d18a0149e1d150e9dec404162fd50abc"
#CA_LA_24h_rep3_2.fastq.gz:
GT_LA_24h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/061/SRR10564761/SRR10564761_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d35fdc18fb5f076b21990515a6feb56f"
#CA_LA_24h_rep4_1.fastq.gz:
GT_LA_24h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/062/SRR10564762/SRR10564762_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep b903a5bf7cf98baae750603c3153c9ac"
#CA_LA_24h_rep4_2.fastq.gz:
GT_LA_24h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/062/SRR10564762/SRR10564762_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d236ca9b669a86f7175e38df59a0a3d2"
#CA_LA_48h_rep1_1.fastq.gz:
GT_LA_48h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/063/SRR10564763/SRR10564763_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 0b5153a5f73be0b51734e12021599a19"
#CA_LA_48h_rep1_2.fastq.gz:
GT_LA_48h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/063/SRR10564763/SRR10564763_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 0de7b02385b07209e7735dfa4d6d09bd"
#CA_LA_48h_rep2_1.fastq.gz:
GT_LA_48h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/064/SRR10564764/SRR10564764_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 03ca516f8400224f304752602f462838"
#CA_LA_48h_rep2_2.fastq.gz:
GT_LA_48h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/064/SRR10564764/SRR10564764_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep e378cf2f238d3148f59e1e990d95ff64"
#CA_LA_48h_rep3_1.fastq.gz:
GT_LA_48h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/065/SRR10564765/SRR10564765_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 829dd08f7829c6447d386a174e7653ff"
#CA_LA_48h_rep3_2.fastq.gz:
GT_LA_48h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/065/SRR10564765/SRR10564765_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep ce49a094adfde4e37a855e42673f6c64"
#CA_LA_48h_rep4_1.fastq.gz:
GT_LA_48h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/066/SRR10564766/SRR10564766_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep ee9d509b50515a82ca1c5ffeefa84c33"
#CA_LA_48h_rep4_2.fastq.gz:
GT_LA_48h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/066/SRR10564766/SRR10564766_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 49b40553ce5eaa444ebf039f3d8aba78"
#CA_LA_no_ISceI_rep1_1.fastq.gz:
GT_LA_no_ISceI_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/067/SRR10564767/SRR10564767_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 69a7ba84b5727185c0dafc82b1d62bec"
#CA_LA_no_ISceI_rep1_2.fastq.gz:
GT_LA_no_ISceI_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/067/SRR10564767/SRR10564767_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8d339c6d9d26e7a9c3803d0ae253186d"
#CA_LA_no_ISceI_rep2_1.fastq.gz:
GT_LA_no_ISceI_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/068/SRR10564768/SRR10564768_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 51b7749e94c48fa2379f7f60ebf5d057"
#CA_LA_no_ISceI_rep2_2.fastq.gz:
GT_LA_no_ISceI_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/068/SRR10564768/SRR10564768_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 6d2ee4ee8e06a20dfe74e5892432d127"
#CA_LA_no_ISceI_rep3_1.fastq.gz:
GT_LA_no_ISceI_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/069/SRR10564769/SRR10564769_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep ed0da856dbe4172cf4f50ebb14a06824"
#CA_LA_no_ISceI_rep3_2.fastq.gz:
GT_LA_no_ISceI_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/069/SRR10564769/SRR10564769_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7ecdb1b1ca6b02bed235263c1a1772fd"
#CA_LA_no_ISceI_rep4_1.fastq.gz:
GT_LA_no_ISceI_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/070/SRR10564770/SRR10564770_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c4edccb8edec8aa20f9d69f5f4360d77"
#CA_LA_no_ISceI_rep4_2.fastq.gz:
GT_LA_no_ISceI_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/070/SRR10564770/SRR10564770_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 73f939d2cda7f8df98dda433afe4e2df"
TC_6xPCR_24h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/071/SRR10564771/SRR10564771_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 6be5e322a8a3564cb44455999c16e6b8"
TC_6xPCR_24h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/071/SRR10564771/SRR10564771_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 48a69d77110ed96e363056923cac9b47"
TC_6xPCR_24h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/072/SRR10564772/SRR10564772_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 71e15250917bdc7aeb965d6225ee1108"
TC_6xPCR_24h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/072/SRR10564772/SRR10564772_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 955c36d30baa9dbd8c2856b9a1be0e13"
TC_6xPCR_24h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/073/SRR10564773/SRR10564773_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep b8558827e9b6f0e8fae532891a35daf5"
TC_6xPCR_24h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/073/SRR10564773/SRR10564773_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 842aec4106c3e9003805add0f3c80a85"
TC_6xPCR_24h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/074/SRR10564774/SRR10564774_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8ecae01924215c825aa8a72abe94ffd8"
TC_6xPCR_24h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/074/SRR10564774/SRR10564774_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 863ae6e7f3e323ccc289b571d607e479"
TC_6xPCR_48h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/075/SRR10564775/SRR10564775_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 32678907282ce58613ed5872cd87c062"
TC_6xPCR_48h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/075/SRR10564775/SRR10564775_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 28885a3feb2cd1c2192efdf105b354a9"
TC_6xPCR_48h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/076/SRR10564776/SRR10564776_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8b87a36f525a95e1dc3de7791eabd89d"
TC_6xPCR_48h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/076/SRR10564776/SRR10564776_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep b67037ad568af0cec9c461b6088e8755"
TC_6xPCR_48h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/077/SRR10564777/SRR10564777_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep e614a00a45c03b38db4a96910c87639c"
TC_6xPCR_48h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/077/SRR10564777/SRR10564777_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8ac69b94f3c00c9ae4e9edd144606135"
TC_6xPCR_48h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/078/SRR10564778/SRR10564778_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 1e1ae3ac153dead84e1da5a667b77b64"
TC_6xPCR_48h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/078/SRR10564778/SRR10564778_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8817b7ea003e49a27f2491e26fdae2bf"
TC_6xPCR_no_ISceI_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/079/SRR10564779/SRR10564779_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2ff8218de02e3fb74d06f08909a88ce6"
TC_6xPCR_no_ISceI_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/079/SRR10564779/SRR10564779_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep b07111c8a64c9c3b7c207770260b6432"
TC_6xPCR_no_ISceI_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/080/SRR10564780/SRR10564780_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 50672d191d5dabf520d7fac93a9960de"
TC_6xPCR_no_ISceI_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/080/SRR10564780/SRR10564780_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 85fc86ce3fdcc0529a33d7f59595c60f"
TC_6xPCR_no_ISceI_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/081/SRR10564781/SRR10564781_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep fb1634bbcb065e860fb2964935f0ff65"
TC_6xPCR_no_ISceI_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/081/SRR10564781/SRR10564781_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep ea9fdba911e5d1d392c1a3ddba2237f0"
TC_6xPCR_no_ISceI_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/082/SRR10564782/SRR10564782_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 6d1ecbf43dd3951677d3198f0426d51b"
TC_6xPCR_no_ISceI_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/082/SRR10564782/SRR10564782_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 1c292ca261d19551f2512edfe8b068c0"
TC_6xPCR_no_ISceI_rep5_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/083/SRR10564783/SRR10564783_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 5a610908de59b787bad231e4c9b5e4aa"
TC_6xPCR_no_ISceI_rep5_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/083/SRR10564783/SRR10564783_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep e8479eeb3af3ddfa203b70361a92a4c5"
TC_6xPCR_no_ISceI_rep6_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/084/SRR10564784/SRR10564784_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7314f41a20a4fbe72de864ab7d1e9f3f"
TC_6xPCR_no_ISceI_rep6_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/084/SRR10564784/SRR10564784_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3d6874e38d81ad1f1120fb67e698ff8a"
TC_6xPCR_no_ISceI_rep7_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/085/SRR10564785/SRR10564785_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep af212e0bd7cffa440b9364982cffb7b0"
TC_6xPCR_no_ISceI_rep7_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/085/SRR10564785/SRR10564785_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2654399d7dc30943a2622a21cf7d4b4f"
TC_6xPCR_no_ISceI_rep8_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/086/SRR10564786/SRR10564786_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 0a3d75fd2cbbeafb62b6fe6b091108df"
TC_6xPCR_no_ISceI_rep8_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/086/SRR10564786/SRR10564786_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 5bd928baaf582c14022212417ce88d49"
TC_LA_24h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/087/SRR10564787/SRR10564787_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 4dc882be293bf02b5d93f8cd32d3bd21"
TC_LA_24h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/087/SRR10564787/SRR10564787_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep b2e4b11b6364ff832316eb6663bb56d2"
TC_LA_24h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/088/SRR10564788/SRR10564788_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 47e8d9ba44c2c0401771e86103caa7d0"
TC_LA_24h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/088/SRR10564788/SRR10564788_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 4b5c2a5858361f5ce3b58ecae6f63f20"
TC_LA_24h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/089/SRR10564789/SRR10564789_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3e2874b2df00798f47ffe707b73e40b4"
TC_LA_24h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/089/SRR10564789/SRR10564789_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d7eaf150434e467e52efdd03918c0764"
TC_LA_24h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/090/SRR10564790/SRR10564790_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep ceffb2ef36f36aed7a748b2b922e11df"
TC_LA_24h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/090/SRR10564790/SRR10564790_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep f8b59d00f935bcf2d685bd9bbff4d91f"
TC_LA_48h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/091/SRR10564791/SRR10564791_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep a32c83796e8924b54d4c768c4356002c"
TC_LA_48h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/091/SRR10564791/SRR10564791_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep ffd9140de08ccfd918c5af587be9a815"
TC_LA_48h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/092/SRR10564792/SRR10564792_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep ae1d7ab1e0ede0055c9f4f020a6b2310"
TC_LA_48h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/092/SRR10564792/SRR10564792_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 816ddfceafa0032db9ddb72267e5bf3e"
TC_LA_48h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/093/SRR10564793/SRR10564793_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep e950f939e362f49332b27c10c3e38f23"
TC_LA_48h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/093/SRR10564793/SRR10564793_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d5f9f01a8269e2e6bdd1c397d8fdd559"
TC_LA_48h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/094/SRR10564794/SRR10564794_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8c67bc2d1976c49ce590f52a2c2d5b69"
TC_LA_48h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/094/SRR10564794/SRR10564794_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 42e993b9c6136408ccbde3dab478a959"
TC_LA_no_ISceI_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/095/SRR10564795/SRR10564795_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c31722e7a766e3c46c57d7094faa1732"
TC_LA_no_ISceI_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/095/SRR10564795/SRR10564795_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep f2a9ac820fa38c93222c1cad086b3f09"
TC_LA_no_ISceI_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/096/SRR10564796/SRR10564796_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 438429c604e49f6425b17db17e5fe959"
TC_LA_no_ISceI_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/096/SRR10564796/SRR10564796_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 4dfc53b22deedcdac639cc038f7b0163"
TC_LA_no_ISceI_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/097/SRR10564797/SRR10564797_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 75c6972eb626f3e8683934e91dbc22cf"
TC_LA_no_ISceI_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/097/SRR10564797/SRR10564797_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 5bd9661ac9636f0d0b7105c912b235c3"
TC_LA_no_ISceI_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/098/SRR10564798/SRR10564798_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8ddba4d7527a8550b14013a526973b66"
TC_LA_no_ISceI_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/098/SRR10564798/SRR10564798_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep dc712971c31e812967fcdc1aab039da4"
TC_LA_no_ISceI_rep5_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/099/SRR10564799/SRR10564799_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep bca6d36d9bb8bc92965fca4094e97979"
TC_LA_no_ISceI_rep5_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/099/SRR10564799/SRR10564799_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 881af27b18712675ea49c770ab1b0bda"
TC_LA_no_ISceI_rep6_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/000/SRR10564800/SRR10564800_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep bb4b857d76ee3e90167a3a3053fd19da"
TC_LA_no_ISceI_rep6_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/000/SRR10564800/SRR10564800_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 5e61155d6a3a902c50b8dcd778f51a1d"
TC_LA_no_ISceI_rep7_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/001/SRR10564801/SRR10564801_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep eef5bd5cedd030400c43ad12c7b7377c"
TC_LA_no_ISceI_rep7_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/001/SRR10564801/SRR10564801_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d21be56123346fdac0d6d173a7b58a65"
TC_LA_no_ISceI_rep8_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/002/SRR10564802/SRR10564802_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 4836e30635c582ec3bf5393e0f5c115f"
TC_LA_no_ISceI_rep8_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/002/SRR10564802/SRR10564802_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 6a07765c1127edc36dbbfda1dfb2c4eb"
TG_6xPCR_24h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/003/SRR10564803/SRR10564803_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7a4d9b253100c95a6cd4f8a19636adab"
TG_6xPCR_24h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/003/SRR10564803/SRR10564803_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d73c961a0044707e545c79ed106c5d8c"
TG_6xPCR_24h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/004/SRR10564804/SRR10564804_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7abd3a064de878d64fd2ee9ede8cc9a4"
TG_6xPCR_24h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/004/SRR10564804/SRR10564804_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3e502fc863ca60ce5e0836b344b36c1b"
TG_6xPCR_24h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/005/SRR10564805/SRR10564805_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep da4ec37d0da9e0c9646ed13c2ea3d30c"
TG_6xPCR_24h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/005/SRR10564805/SRR10564805_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 43cf314760837f7b95260f621e8b5402"
TG_6xPCR_24h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/006/SRR10564806/SRR10564806_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3d2a99ffea6d326ee4e14e3e3384f32d"
TG_6xPCR_24h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/006/SRR10564806/SRR10564806_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 490c15ed527259b2cfc7e438acf8f244"
TG_6xPCR_48h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/007/SRR10564807/SRR10564807_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c8d3bc8fa5da3684cafa805a99bf85f7"
TG_6xPCR_48h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/007/SRR10564807/SRR10564807_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 7331276f2eb979216a6909c06b4ac843"
TG_6xPCR_48h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/030/SRR10564730/SRR10564730_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 63a8cd45b6db572157ce2068a80f646b"
TG_6xPCR_48h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/030/SRR10564730/SRR10564730_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 1bb6a802de57ea962797918e9b088f4d"
TG_6xPCR_48h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/031/SRR10564731/SRR10564731_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 31ccc2231e5745ffd30a04719359fb6f"
TG_6xPCR_48h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/031/SRR10564731/SRR10564731_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 265d35b79388ccb24d4c78c3db601be4"
TG_6xPCR_48h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/032/SRR10564732/SRR10564732_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep ecd6552b169b4e6d0498da919ed738cd"
TG_6xPCR_48h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/032/SRR10564732/SRR10564732_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep e04493703f15e4c2a176771860a7c889"
TG_6xPCR_no_ISceI_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/033/SRR10564733/SRR10564733_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 86280d3fdfc5905d002155ea85f5f8be"
TG_6xPCR_no_ISceI_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/033/SRR10564733/SRR10564733_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep d7f4e4314d87068081ba02c0a9ae0e48"
TG_6xPCR_no_ISceI_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/034/SRR10564734/SRR10564734_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep f2789d839ce270e9f27ed28ed45373f2"
TG_6xPCR_no_ISceI_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/034/SRR10564734/SRR10564734_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 0d5f256e38c9da4e070154a7d791144c"
TG_6xPCR_no_ISceI_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/035/SRR10564735/SRR10564735_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 9008636ac04ef77a59c7328a9698a3ac"
TG_6xPCR_no_ISceI_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/035/SRR10564735/SRR10564735_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 1d4c275f27e429554eafe17ca2636c60"
TG_6xPCR_no_ISceI_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/036/SRR10564736/SRR10564736_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 2e7f102c0c8fa4afdb9596ce0390c96d"
TG_6xPCR_no_ISceI_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/036/SRR10564736/SRR10564736_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep b955b495a51895a7a070336576f652a2"
TG_LA_24h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/037/SRR10564737/SRR10564737_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep aa5320499ba58df6e894c5003402e903"
TG_LA_24h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/037/SRR10564737/SRR10564737_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 89c2527bb8407d5d87ab17ecca953e03"
TG_LA_24h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/038/SRR10564738/SRR10564738_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 992eeb0a4c70fd30dc0f615ff26e7f36"
TG_LA_24h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/038/SRR10564738/SRR10564738_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep cd9288cb1699285aad431d02a332f695"
TG_LA_24h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/039/SRR10564739/SRR10564739_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep de14dab805cc4a0449029a3a5e71994b"
TG_LA_24h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/039/SRR10564739/SRR10564739_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep aec4326380228b0d5280b8ab27795611"
TG_LA_24h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/040/SRR10564740/SRR10564740_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 896a43cc323571913d9c28483bcaf19b"
TG_LA_24h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/040/SRR10564740/SRR10564740_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 6ff02146e81820b6cf96ef9a64ce81e6"
TG_LA_48h_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/041/SRR10564741/SRR10564741_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 10c0684ec1683688bba30e7a75dd81c2"
TG_LA_48h_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/041/SRR10564741/SRR10564741_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 71026d5e4050d8645592d5a628562f16"
TG_LA_48h_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/042/SRR10564742/SRR10564742_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c076fa8ca8a58caaf8ec24afd4714055"
TG_LA_48h_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/042/SRR10564742/SRR10564742_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c4fa8b0bf6071d7e7ebab5e4e95aa103"
TG_LA_48h_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/043/SRR10564743/SRR10564743_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep cebeb9021aa77f25db9e39b4c5b8693b"
TG_LA_48h_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/043/SRR10564743/SRR10564743_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep c56e3d3acb7703af33ebe184263d5854"
TG_LA_48h_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/044/SRR10564744/SRR10564744_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 8a0d5838d2a0b873c540b95c3bb22b85"
TG_LA_48h_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/044/SRR10564744/SRR10564744_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 94621f46703a803951187fdb52b27c55"
TG_LA_no_ISceI_rep1_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/045/SRR10564745/SRR10564745_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep b3cc2d4ed0cbb9e6d248dbf6fb55ee2f"
TG_LA_no_ISceI_rep1_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/045/SRR10564745/SRR10564745_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 5927a5c07abd3b2aa1fcefca25dd7ecb"
TG_LA_no_ISceI_rep2_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/046/SRR10564746/SRR10564746_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep a979fe0892be87cd4e9501935fac38c8"
TG_LA_no_ISceI_rep2_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/046/SRR10564746/SRR10564746_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 0233d11bd2f5dfc50887c541b399f806"
TG_LA_no_ISceI_rep3_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/047/SRR10564747/SRR10564747_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep e066c7fb452764c393854632db615df6"
TG_LA_no_ISceI_rep3_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/047/SRR10564747/SRR10564747_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep a66e26933f6579d4e269fc2ab1f4906e"
TG_LA_no_ISceI_rep4_1.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/048/SRR10564748/SRR10564748_1.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 39a59ed4a14157347340be42add91dc4"
TG_LA_no_ISceI_rep4_2.fastq.gz:
	$(DOCKER_RUN) wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/048/SRR10564748/SRR10564748_2.fastq.gz -O $@
	$(DOCKER_RUN) /bin/bash -c "md5sum $@ | grep 3100685bb564c19a3f6df14e6d642993"


%.co.gz: %_1.fastq.gz %_2.fastq.gz
	$(DOCKER_RUN) /bin/bash -c "python scripts/preprocess.py $^ | gzip > $*.pps.gz"
	$(DOCKER_RUN) /bin/bash -c "zcat $*.pps.gz | cut -f1 | starcode -d2 --print-clusters | gzip > $*.stc.gz"
	$(DOCKER_RUN) /bin/bash -c "python scripts/count.py $*.pps.gz $*.stc.gz $*_seq.txt.gz | gzip > $@ && rm $*.pps.gz $*.stc.gz"

check:
	md5sum -c md5sums.txt

clean:
	rm -f *.co.gz *_seq.txt.gz *fastq.gz
