# Makefile for Statistics-Tool-Box example programs and test suite

CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c99 -D_GNU_SOURCE -I.
LDFLAGS = -lm -lz

# Example programs
TARGET = dim_reduce
TARGET2 = deseq2_example
TARGET5 = spearman
SOURCE = examples/dim_reduce.c
SOURCE2 = examples/deseq2_example.c
SOURCE5 = examples/spearman.c

# Test programs
TARGET3 = test_stb_stats
TARGET4 = test_isolated
SOURCE3 = test_stb_stats.c
SOURCE4 = test_isolated.c

# Header file
HEADER = stb_stats.h

all: $(TARGET) $(TARGET2) $(TARGET5) $(TARGET3) $(TARGET4)

$(TARGET): $(SOURCE) $(HEADER)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCE) $(LDFLAGS)

$(TARGET2): $(SOURCE2) $(HEADER)
	$(CC) $(CFLAGS) -o $(TARGET2) $(SOURCE2) $(LDFLAGS)

$(TARGET3): $(SOURCE3) $(HEADER)
	$(CC) $(CFLAGS) -o $(TARGET3) $(SOURCE3) $(LDFLAGS)

$(TARGET4): $(SOURCE4) $(HEADER)
	$(CC) $(CFLAGS) -o $(TARGET4) $(SOURCE4) $(LDFLAGS)

$(TARGET5): $(SOURCE5) $(HEADER)
	$(CC) $(CFLAGS) -o $(TARGET5) $(SOURCE5) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) *.o

test: $(TARGET3)
	@echo "Running comprehensive stb_stats test suite..."
	./$(TARGET3)

test-isolated: $(TARGET4)
	@echo "Running isolated tests for problematic functions..."
	./$(TARGET4)

test-dim: $(TARGET)
	@echo "Testing PCA..."
	./$(TARGET) test_data.txt -a pca -o test_pca.txt
	@echo "Testing t-SNE..."
	./$(TARGET) test_data.txt -a tsne -o test_tsne.txt
	@echo "Testing UMAP..."
	./$(TARGET) test_data.txt -a umap -o test_umap.txt

.PHONY: all clean test test-isolated test-dim
