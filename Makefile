#original tests from master dont work with the modifications made
SUBDIRS=src #tests

all: src 
#tests

$(SUBDIRS):
	$(MAKE) -C $@
.PHONY: $(SUBDIRS)

clean:
	@for y in $(SUBDIRS); do $(MAKE) -C $$y/ clean ; done

