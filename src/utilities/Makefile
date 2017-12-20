SRC_DIR = .

BUILD_DIR = ../..
include $(BUILD_DIR)/config/Makefile.config

all:  includes library

includes:  
	cp -f $(SRC_DIR)/qss_data_arrays.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_file.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_general_util.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_options.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_macros.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_util2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/qss_util3d.h $(BUILD_DIR)/include/


library:   qss_data_arrays.o         \
	    qss_file.o         \
	    qss_general_util.o         \
	    qss_options.o       \
	    qss_util2d.o        \
	    qss_util3d.o

clean:
	rm -f *.o 
