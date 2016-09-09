NAME    = atomorph-gui
CC      = gcc
PROF    = -O2 
C_FLAGS = -std=c++11 -Wall $(PROF)
L_FLAGS = ../../libatomorph.a -lm -lstdc++ $(PROF)
OBJ_DIR = obj
DEFINES = -D ATOMORPH_DEPRECATED

#Change the LIB_DIR to the path where the MODULES are located:
LIB_DIR = /usr/lib/
MODULES = allegro.so allegro_image.so allegro_font.so allegro_primitives.so

SRC_FILES := $(wildcard *.cpp)
LIB_FILES := $(patsubst %,$(LIB_DIR)/lib%,$(MODULES))
O_FILES   := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

OUT     = ./$(NAME)

all:
	@printf "\033[0mHINT: libatomorph.a has to be compiled using the \033[1;33mmake deprecated\033[0m command.\033[0m\n"
	@printf "\033[0mHINT: Use \033[1;33mmake static\033[0m to link statically.\033[0m\n"
	@$(MAKE) dynamically_linked -s	

static:
	@$(MAKE) statically_linked -s

statically_linked: $(O_FILES)
	@printf "\033[1;33mMaking \033[37m   ...."
	$(CC) -static -o ./$(NAME) $(O_FILES) $(LIB_FILES) $(L_FLAGS)
	@printf "\033[1;32m Statically linked %s DONE!\033[0m\n" $(NAME)

dynamically_linked: $(O_FILES)
	@printf "\033[1;33mMaking \033[37m   ...."
	$(CC) -o ./$(NAME) $(O_FILES) $(LIB_FILES) $(L_FLAGS)
	@printf "\033[1;32m Dynamically linked %s DONE!\033[0m\n" $(NAME)

$(OBJ_DIR)/%.o: %.cpp
		@printf "\033[1m\033[31mCompiling \033[37m....\033[34m %-20s\t\033[33m%6s\033[31m lines\033[0m \n" $*.cpp "`wc -l $*.cpp | cut -f1 -d' '`"
		@$(CC) $< $(DEFINES) $(C_FLAGS) -c -o $@

clean:
	@printf "\033[1;36mCleaning \033[37m ...."
	@rm -f $(O_FILES) $(OUT) *~ *.bak *.orig *.rej
	@printf "\033[1;37m $(NAME) cleaned!\033[0m\n"
