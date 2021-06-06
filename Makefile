CC=gcc
CFLAGS=-Wall-std=gnu99

TARGET=factotial
SRCS=alg2compress.c Huffman_decode.c Huffman_encode.c

OBJS=$(SRCS:.c=.o)

all:$(TARGET)

$(TARGET):$(OBJS)
$(CC)=0 $@ $^

.PHONY:clean
clean:
rm-rf $(TARFET) $(OBJS)

%.o:%.c
$(CC) $(CFLAGS)-o $@-c $<