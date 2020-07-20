TARGET   = gvcf2coverage

CC       = gcc
CFLAGS   = -std=c99 -march=native -Wall -Wextra -Wpedantic \
           -Wformat=2 -Wshadow -Wwrite-strings -Wstrict-prototypes \
           -Wold-style-definition -Wredundant-decls -Wnested-externs \
           -O3 -DNDEBUG


$(TARGET): src/main.c
	$(CC) -o $@ $(CFLAGS) -I$(HTSLIB_INCDIR) -L$(HTSLIB_LIBDIR) $< -lhts

.PHONY: clean

clean:
	rm $(TARGET)
