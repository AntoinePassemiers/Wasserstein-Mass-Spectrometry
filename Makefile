all: src

src:
	$(MAKE) -C $@

.PHONY: src release vectorized clean

release:
	cd src && $(MAKE) release

clean:
	cd src && $(MAKE) clean