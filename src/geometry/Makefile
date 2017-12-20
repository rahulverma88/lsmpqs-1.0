SRC_DIR = .

BUILD_DIR = ../..
include $(BUILD_DIR)/config/Makefile.config

all:  includes library

includes:  
	cp -f $(SRC_DIR)/qss_grid.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_initialization2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_initialization3d.h $(BUILD_DIR)/include/

library:   qss_grid.o         \
	    qss_initialization2d.o         \
	    qss_initialization3d.o

clean:
	rm -f *.o 
