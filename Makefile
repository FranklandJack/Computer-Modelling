# Makefile for the Solar-System-Simulation project

SRC_DIR=source
OUT_DIR=output

SRC_FILES=$(wildcard $(SRC_DIR)/*.java)
CLS_FILES=$(patsubst $(SRC_DIR)/%.java, %.class, $(SRC_FILES))
OUT_FILES=$(wildcard $(OUT_DIR)/*.txt)


## all       : Compile code and create executible 
.PHONY : all
all : $(CLS_FILES)

$(CLS_FILES) : $(SRC_FILES)
	javac -d . $(SRC_FILES)


## clean     : remove auto generated files
.PHONY : clean
clean :
	rm -f $(CLS_FILES)
	rm -f $(OUT_FILES)

## variables : Print variables
.PHONY :variables
variables:
	@echo SRC_DIR:   $(SRC_DIR)
	@echo SRC_FILES: $(SRC_FILES)
	@echo OBJ_FILES: $(CLS_FILES)
	


## help      : Print help
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<