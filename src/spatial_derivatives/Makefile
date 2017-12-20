SRC_DIR = .

BUILD_DIR = ../..
include $(BUILD_DIR)/config/Makefile.config

all:  includes library

includes:  
	cp -f $(SRC_DIR)/qss_spatial_derivatives2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_spatial_derivatives3d.h $(BUILD_DIR)/include/

library:  qss_spatial_derivatives2d.o       \
          qss_spatial_derivatives3d.o       

clean:
	rm -f *.o 


