SRC_DIR = .

BUILD_DIR = ../..
include $(BUILD_DIR)/config/Makefile.config

# Set MAKE
    

all:  includes library

includes: 
	cp -f $(SRC_DIR)/drain_model2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/drain_model3d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/drain_top.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/imbibe_model2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/imbibe_model3d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/imbibe_top.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/constCurvModel2d.h $(BUILD_DIR)/include/
	cp -f $(SRC_DIR)/constCurvModel3d.h $(BUILD_DIR)/include/

library: constCurvModel2d.o         \
	constCurvModel3d.o      \
	drain_model2d.o         \
	drain_model3d.o         \
	drain_top.o             \
	imbibe_model2d.o         \
	imbibe_model3d.o         \
	imbibe_top.o

clean:
	$(RM) *.o
