HEADER_DIR = headers
SRC_DIR = src
LIB_DIR = lib

CC ?= gcc
CFLAGS ?= -I$(HEADER_DIR)

TARGET = main.exe

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(LIB_DIR)/%.o, $(SRCS))

all: $(TARGET) clean_aux

libcml: $(LIB_DIR)/libcml.a

libcml-basic: $(LIB_DIR)/libcml-basic.a


$(TARGET): main.o libcml libcml-basic
	$(CC) -o $@ main.o -L$(LIB_DIR) -lcml

main.o: main.c | $(LIB_DIR)
	$(CC) $(CFLAGS) -c main.c -o $@

$(LIB_DIR)/libcml.a: $(LIB_DIR)/CMatLib.o
	ar rcs $@ $^

$(LIB_DIR)/libcml-basic.a: $(LIB_DIR)/basic.o
	ar rcs $@ $^

$(LIB_DIR)/CMatLib.o: $(SRC_DIR)/CMatLib.c $(LIB_DIR)/basic.o | $(LIB_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB_DIR)/basic.o: $(SRC_DIR)/basic.c | $(LIB_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB_DIR):
	test -d $(LIB_DIR) || mkdir $(LIB_DIR)

clean_aux:
	rm -f main.o $(LIB_DIR)/*.o

clean:
	rm -rf $(TARGET) main.o $(LIB_DIR)/*.o $(LIB_DIR)/*.a
