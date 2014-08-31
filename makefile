#compiler
CXX=g++

#compiler flags
COPTI= -O3 -std=c++11
CDEBG= -g -std=c++11
CWARN= -Wall -Wextra -Wno-sign-compare
CBAMH= -D_LIBBAM
CSHV1=

#linker flags
LOPTI= -O3 -std=c++11
LDEBG= -g -std=c++11
LSTDD= -lm -lpthread -lboost_iostreams -lboost_program_options
LSTDS= -Wl,-Bstatic -lboost_iostreams -lz -lbz2 -Wl,-Bdynamic -lm -lpthread
LCL3S= -Wl,-Bstatic -lz -lbz2 -Wl,-Bdynamic -lm -lpthread
LBAMD= -lbam -lz

#versions
VMAJ = $(shell grep VERSION_MAJOR src/version.h | cut -f3 -d' ')
VMIN = $(shell grep VERSION_MINOR src/version.h | cut -f3 -d' ')

#executable file
EFBASE:=bin/measureAggregateRsquared
EFILE:= $(EFBASE).$(VMAJ).$(VMIN)
EDFILE:=$(EFBASE).$(VMAJ).$(VMIN).dbg


#header files
HFILE= $(shell find src -name *.h)

#source files
CFILE= $(shell find src -name *.cpp)

#source path
VPATH= $(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#include path
ISTDP= -Isrc
IBAMP= -Ilib
ICL3P= -I/users/delaneau/BOOST/include

#library path
LSTDP= -Llib

#object files
OFILE= $(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
OBOST=

#default
all: dynamic

#dynamic release
dynamic: CFLAG=$(COPTI) $(CWARN) $(CSHV1)
dynamic: LFLAG=$(LOPTI) $(LSTDD)
dynamic: IFLAG=$(ISTDP)
dynamic: $(EFILE)

# testing release
test: CFLAG=$(CDEBG) $(CWARN) $(CSHV1) $(CDEBG)
test: LFLAG=$(LDEBG) $(LSTDD)
test: IFLAG=$(ISTDP)
test: $(EDFILE)

#static release
static: CFLAG=$(COPTI) $(CWARN) $(CSHV1)
static: LFLAG=$(LOPTI) $(LSTDS)
static: IFLAG=$(ISTDP)
static: $(EFILE)

#cluster release
cluster: CFLAG=$(COPTI) $(CWARN) $(CSHV1)
cluster: LFLAG=$(LOPTI) $(LCL3S)
cluster: IFLAG=$(ISTDP) $(ICL3P)
cluster: OBOST=~/BOOST/lib/libboost_iostreams.a ~/BOOST/lib/libboost_program_options.a
cluster: $(EFILE)

$(EFILE): $(OFILE)
	$(CXX) $^ $(OBOST) -o $@ $(LFLAG)
	rm -f $(EFBASE) && ln -s $(shell basename $@) $(EFBASE)

$(EDFILE): $(OFILE)
	$(CXX) $^ $(OBOST) -o $@ $(LFLAG)
	rm -f $(EFBASE).dbg && ln -s $(shell basename $@) $(EFBASE).dbg

obj/%.o: %.cpp $(HFILE)
	$(CXX) -o $@ -c $< $(CFLAG) $(IFLAG)

clean: 
	rm -f obj/*.o $(EFILE) $(EDFILE)

test: $(EDFILE)
	perl t/runtests.pl

oxford:
	cp $(EFILE) ~/bin/
	rm -f ~/$(EFBASE) && ln -s $(shell basename $(EFILE)) ~/bin/$(shell basename $(EFBASE)) 

install:
	cp $(EFILE) /usr/local/bin/.
