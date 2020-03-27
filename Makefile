all:
	##################################################
	###               MapleCOMSPS                  ###
	##################################################
	if [ -d mapleCOMSPS/m4ri-20140914 ]; then : ; \
	else cd mapleCOMSPS && tar zxvf m4ri-20140914.tar.gz && \
	cd m4ri-20140914 && ./configure; fi
	+ $(MAKE) -C mapleCOMSPS/m4ri-20140914
	+ $(MAKE) -C mapleCOMSPS

	##################################################
	###                 PaInleSS                   ###
	##################################################
	+ $(MAKE) -C painless-src
	mv painless-src/painless painless-mcomsps

clean:
	##################################################
	###               MapleCOMSPS                  ###
	##################################################
	rm -rf mapleCOMSPS/m4ri-20140914
	+ $(MAKE) -C mapleCOMSPS clean

	##################################################
	###                 PaInleSS                   ###
	##################################################
	+ $(MAKE) clean -C painless-src
	rm -f painless-mcomsps
