CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -Wall
BIN=.

all: csvExtractRepeatsAndBorders reverseComplement indexmaker

reverseComplement: reverseComplement.c
	$(CC) $(CFLAGS) reverseComplement.c commonFunctions.c -o $(BIN)/reverseComplement

indexmaker: indexmaker.c
	$(CC) $(CFLAGS) indexmaker.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/indexmaker

csvExtractRepeatsAndBorders: csvExtractRepeatsAndBorders.c
	$(CC) $(CFLAGS) csvExtractRepeatsAndBorders.c fragmentv2.c commonFunctions.c comparisonFunctions.c -lm -o $(BIN)/csvExtractRepeatsAndBorders

clean:
	rm -rf $(BIN)/csvExtractRepeatsAndBorders $(BIN)/indexmaker $(BIN)/reverseComplement
