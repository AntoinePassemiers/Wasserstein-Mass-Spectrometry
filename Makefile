all: src

src:
	$(MAKE) -C $@

.PHONY: src debug release clean wassms deconvms

debug:
	cd src && $(MAKE) debug

release:
	cd src && $(MAKE) release

wassms:
	cd src && $(MAKE) wassms

deconvms:
	cd src && $(MAKE) deconvms

clean:
	cd src && $(MAKE) clean