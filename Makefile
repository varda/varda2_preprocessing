.PHONY: all

all:
	gcc -Wall -Wextra -Wpedantic -std=c99 -O3 -I/exports/rse/varda2/samtools/include -L/exports/rse/varda2/samtools/lib main.c -lhts
