FUNC := gcc
copt := -c
OBJ_DIR := ./bin/
FLAGS := -O3 -lm -g -Werror -fopenmp -std=gnu11 -fopt-info-vec-all

C_FILES := $(wildcard src/*.c)
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(C_FILES:.c=.obj)))

all:
	$(FUNC) ./main.c -o ./main.exe $(FLAGS)

clean:
	rm -f ./*.exe