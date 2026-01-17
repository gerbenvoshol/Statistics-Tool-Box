# Makefile for dim_reduce example program

CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c99 -D_GNU_SOURCE
LDFLAGS = -lm -lz

TARGET = dim_reduce
SOURCE = dim_reduce.c
HEADER = stb_stats.h

all: $(TARGET)

$(TARGET): $(SOURCE) $(HEADER)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCE) $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o

test: $(TARGET)
	@echo "Testing PCA..."
	./$(TARGET) test_data.txt -a pca -o test_pca.txt
	@echo "Testing t-SNE..."
	./$(TARGET) test_data.txt -a tsne -o test_tsne.txt
	@echo "Testing UMAP..."
	./$(TARGET) test_data.txt -a umap -o test_umap.txt

.PHONY: all clean test
