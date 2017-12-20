SRC_DIR = .

BUILD_DIR = ../..
include $(BUILD_DIR)/config/Makefile.config

all:  includes library

includes:  
	cp -f $(SRC_DIR)/qss_reinitialization2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_reinitialization3d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/reinitialize_top.h $(BUILD_DIR)/include/

library:  qss_reinitialization2d.o         \
	   qss_reinitialization3d.o         \
	   reinitialize_top.o       

clean:
	rm -f *.o 
