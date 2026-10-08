TARGET   = gvcf2coverage

CC       = gcc
CFLAGS   = -std=c99 -march=native -Wall -Wextra -Wpedantic \
           -Wformat=2 -Wshadow -Wwrite-strings -Wstrict-prototypes \
           -Wold-style-definition -Wredundant-decls -Wnested-externs \
           -O3 -DNDEBUG

.PHONY: clean install

$(TARGET): src/main.c
	$(CC) -o $@ $(CFLAGS) -I$(HTSLIB_INCDIR) -L$(HTSLIB_LIBDIR) $< -lhts

clean:
	rm $(TARGET)

install: $(TARGET)
	mkdir -p $(PREFIX)/bin
	install -m 755 $< $(PREFIX)/bin
