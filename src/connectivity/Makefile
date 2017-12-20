SRC_DIR = .

BUILD_DIR = ../..
include $(BUILD_DIR)/config/Makefile.config

all:  includes library

includes:
	cp -f $(SRC_DIR)/connectivity.h $(BUILD_DIR)/include/

library: connectivity.o trapComponents.o

clean:
	$(RM) *.o 

