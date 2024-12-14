#optional: use a custom build directory 
ifdef PREFIX
    VASP_BUILD_DIR=$(PREFIX)
else
    VASP_BUILD_DIR=build
endif

VERSIONS = std gam ncl
.PHONY: all veryclean test test_all versions $(VERSIONS)
all: std gam ncl
versions: $(VERSIONS)
$(VERSIONS):
	if [ ! -d $(VASP_BUILD_DIR)/$@ ] ; then mkdir -p $(VASP_BUILD_DIR)/$@  ; fi
	cp src/makefile src/.objects src/makedeps.awk makefile.include $(VASP_BUILD_DIR)/$@ 

	$(MAKE) -C $(VASP_BUILD_DIR)/$@ VERSION=$@ check

ifdef DEPS
	$(MAKE) -C $(VASP_BUILD_DIR)/$@ VERSION=$@ dependencies -j1
else
	$(MAKE) -C $(VASP_BUILD_DIR)/$@ VERSION=$@ cleandependencies -j1
endif


ifdef MODS
	$(MAKE) -C $(VASP_BUILD_DIR)/$@ VERSION=$@ modfiles -j1
endif
	$(MAKE) -C $(VASP_BUILD_DIR)/$@ VERSION=$@ all

veryclean: 
	rm -rf $(VASP_BUILD_DIR)/std
	rm -rf $(VASP_BUILD_DIR)/gam
	rm -rf $(VASP_BUILD_DIR)/ncl

test:
	$(MAKE) -C testsuite test

test_all:
	$(MAKE) -C testsuite test_all
