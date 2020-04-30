.PHONY: install

gvcf2coverage: src/main.c
	$(CC) -o $@ -Wall -Wextra -Wpedantic -std=c99 -O3 -I$(HTSLIB_INCDIR) -L$(HTSLIB_LIBDIR) src/main.c -lhts

install: gvcf2coverage
	mkdir -p $(PREFIX)/bin
	install -m 755 gvcf2coverage $(PREFIX)/bin
