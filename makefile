PROGRAM = tedna
C_FILES := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp, %.o, $(C_FILES))
#CC = clang++
CC = g++
# CFLAGS = -Wall -fstack-protector-all -Wstack-protector -fno-omit-frame-pointer -pedantic -std=c++11 -pthread -fsanitize=thread
CFLAGS = -Wall -fstack-protector-all -Wstack-protector -fno-omit-frame-pointer -pedantic -std=c++11 -pthread
# CFLAGS = -Wall -pedantic -std=c++11 -pthread
LDFLAGS =
# LDFLAGS = -static -static-libgcc -static-libstdc++

k ?= 61
NB_BLOCKS := $(shell echo \(${k}+1\)/32+1 | bc)
CFLAGS += -DPRE_NB_BLOCKS=$(NB_BLOCKS)

ifdef DEBUG
	# CFLAGS += -O0 -p -pg -g -ggdb
	CFLAGS += -O0 -g -ggdb
else
	CFLAGS += -O3
endif

HASH ?= SLOW
ifeq ($(HASH), SLOW)
	CFLAGS += -DHASH_SLOW
endif
ifeq ($(HASH), MID)
	CFLAGS += -DHASH_MID
endif
ifeq ($(HASH), FAST)
	CFLAGS += -DHASH_FAST
endif

all: $(PROGRAM)

$(PROGRAM): depend $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(PROGRAM) $(CFLAGS) 

depend: .depend

release: .release

.depend: cmd = $(CC) -std=c++11 -MM -MF depend $(var); cat depend >> .depend;
.depend:
	@echo "Generating dependencies..."
	@$(foreach var, $(C_FILES), $(cmd))
	@rm -f depend

.release: rel = tedna_$(RELEASE)
.release:
	@mkdir $(rel)
	@cp -r makefile *.cpp *.hpp *.h Evaluator sparsehash NEWS README COPYING ChangeLog AUTHORS INSTALL $(rel)
	@tar cvf $(rel).tar $(rel)
	@gzip -9 $(rel).tar
	@rm -rf $(rel)

-include .depend

# These are the pattern matching rules. In addition to the automatic
# variables used here, the variable $* that matches whatever % stands for
# can be useful in special cases.
%.o: %.cpp
	$(CC) -c $< -o $@ $(CFLAGS) 

clean:
	rm -f .depend $(OBJS) $(PROGRAM)

.PHONY: clean depend
