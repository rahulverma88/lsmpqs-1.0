SRC_DIR = .

BUILD_DIR = ..
include $(BUILD_DIR)/config/Makefile.config

# Set MAKE

APPL_DIR = applications

SUBDIRS = geometry                  \
          level_set_evolution       \
          reinitialization          \
          spatial_derivatives       \
          time_integration          \
          connectivity              \
          utilities 
          
all:  includes library progs

includes: 
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

library: 
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

progs:
	cd $(APPL_DIR); make $@ || exit 1

clean:
	rm -f *.o
	for DIR in $(SUBDIRS); do (cd $$DIR; make $@) || exit 1; done

spotless: clean
	make clean
