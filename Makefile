SRC_DIR = .
VPATH = .
BUILD_DIR = .
include $(BUILD_DIR)/config/Makefile.config

# Set MAKE


all:  library progs 

directories:
	if [ ! -d "include" ]; then mkdir include; fi
	if [ ! -d "lib" ]; then mkdir lib; fi

includes:  directories
	cp -f config/QSSLIB_config.h include
	cd src; make $@ || exit 1


library:  directories includes 
	cd src; make $@ || exit 1
	make lsmpqs_lib || exit 1

lsmpqs_lib:
	find src -name "*.o" > objs_file.tmp
	cat objs_file.tmp | xargs ar -ru $(BUILD_DIR)/lib/liblsmpqs.a 
	ranlib lib/liblsmpqs.a
	rm -f objs_file.tmp 

progs:
	cd src; make $@ || exit 1

clean:
	rm -f lib/*
	cd src; make $@ || exit 1
	rm -f *.o
	rm -f *.tmp

spotless: clean
	cd src; make $@ || exit 1
	cd examples; make $@ || exit 1
	if [ -d "include" ]; then rm -f -rf include/*; rmdir include; fi
	if [ -d "lib" ]; then rm -f -rf lib/*; rmdir lib; fi

