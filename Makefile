.PHONY: all clean buildall

SRC_DIR = src

all:
	cd $(SRC_DIR) && $(MAKE)

clean:
	cd $(SRC_DIR) && $(MAKE) clean

buildall: clean all
