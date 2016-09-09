NAME    = atomorph
CC      = gcc
PROF    = -O2 
C_FLAGS = -std=c++11 -Wall $(PROF)
L_FLAGS = -lm -lstdc++ -pthread $(PROF)
OBJ_DIR = obj
OPENCV  = 0

#Change the LIB_DIR to the path where the MODULES are located:
LIB_DIR = ../
MODULES = atomorph.a

SRC_FILES := $(wildcard *.cpp)
LIB_FILES := $(patsubst %,$(LIB_DIR)/lib%,$(MODULES))
O_FILES   := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

OUT = ./$(NAME)

all:
	@printf "\033[0mHINT: On errors, try \033[1;33m-std=gnu++11 -stdlib=libc++\033[0m compiler flags.\033[0m\n"
	@printf "\033[0mHINT: On linker errors, try \033[1;33mmake opencv\033[0m to link with OpenCV.\033[0m\n"
	@$(MAKE) independently_linked -s

opencv:
	@$(MAKE) opencv_linked -s

opencv_linked: L_FLAGS+= -lopencv_core -lopencv_flann
opencv_linked: $(O_FILES)
	@printf "\033[1;33mMaking \033[37m   ...."
	$(CC) -o ./$(NAME) $(O_FILES) $(LIB_FILES) $(L_FLAGS)
	@printf "\033[1;32m OpenCV dependent %s DONE!\033[0m\n" $(NAME)

independently_linked: $(O_FILES)
	@printf "\033[1;33mMaking \033[37m   ...."
	$(CC) -o ./$(NAME) $(O_FILES) $(LIB_FILES) $(L_FLAGS)
	@printf "\033[1;32m %s DONE!\033[0m\n" $(NAME)

$(OBJ_DIR)/%.o: %.cpp
		@printf "\033[1m\033[31mCompiling \033[37m....\033[34m %-20s\t\033[33m%6s\033[31m lines\033[0m \n" $*.cpp "`wc -l $*.cpp | cut -f1 -d' '`"
		@$(CC) $< $(DEFINES) $(C_FLAGS) -c -o $@

clean:
	@printf "\033[1;36mCleaning \033[37m ...."
	@rm -f $(O_FILES) $(OUT) *~ *.bak *.orig *.rej
	@printf "\033[1;37m $(NAME) cleaned!\033[0m\n"
