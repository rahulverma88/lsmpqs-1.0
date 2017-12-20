SRC_DIR = .

BUILD_DIR = ../..
include $(BUILD_DIR)/config/Makefile.config

all:  includes library

includes:  
	cp -f $(SRC_DIR)/qss_tvd_runge_kutta2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_tvd_runge_kutta3d.h $(BUILD_DIR)/include/

library:  qss_tvd_runge_kutta2d.o       \
          qss_tvd_runge_kutta3d.o       

clean:
	rm -f *.o 

