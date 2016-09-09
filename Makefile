NAME    = atomorph
CC      = gcc
PROF    = -O2
C_FLAGS = -std=c++11 -Wall -pedantic $(PROF)
OBJ_DIR = obj

SRC_FILES := $(wildcard *.cpp)
O_FILES   := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

OUT = ./lib$(NAME).a

all:
	@printf "\033[0mHINT: On errors, try \033[1;33m-std=gnu++11 -stdlib=libc++\033[0m compiler flags.\033[0m\n"
	@printf "\033[0mHINT: Use \033[1;33mmake opencv\033[0m for OpenCV optimizations (experimental).\033[0m\n"
	@printf "\033[0mHINT: Use \033[1;33mmake deprecated\033[0m to compile the old version.\033[0m\n"
	@$(MAKE) $(OUT) -s

opencv: DEFINES = -D ATOMORPH_OPENCV
opencv: $(O_FILES)
	@ar rcs $(OUT) $(O_FILES)
	@printf "\033[1;32mOpenCV dependent lib$(NAME).a DONE!\033[0m\n"

deprecated: DEFINES = -D ATOMORPH_DEPRECATED
deprecated: $(O_FILES)
	@ar rcs $(OUT) $(O_FILES)
	@printf "\033[1;32mDeprecated lib$(NAME).a DONE!\033[0m\n"

$(OUT): $(O_FILES)
	@ar rcs $(OUT) $(O_FILES)
	@printf "\033[1;32mlib$(NAME).a DONE!\033[0m\n"

$(OBJ_DIR)/%.o: %.cpp
		@printf "\033[1m\033[31mCompiling \033[37m....\033[34m %-20s\t\033[33m%6s\033[31m lines\033[0m \n" $*.cpp "`wc -l $*.cpp | cut -f1 -d' '`"
		@$(CC) $(INCLUDE) $< $(DEFINES) $(C_FLAGS) -c -o $@

clean:
	@printf "\033[1;36mCleaning \033[37m ...."
	@rm -f $(O_FILES) $(OUT) *~ *.bak *.orig *.rej
	@printf "\033[1;37m lib$(NAME).a cleaned!\033[0m\n"
