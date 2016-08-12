VIS=\
		Visualizer

SRC=\
		Spatiocyte\
		Model\
		Stepper\
		Compartment\
		Lattice\
		Species\
		Diffuser\
		VisualLogger\
		Random\
		Optimization

IFLAGS = -I. -I$(HOME)/root/include -I/usr/include/x86_64-linux-gnu/c++/5.2.1
LDFLAGS = # -L$(HOME)/root/lib -lRandom
#CXXFLAGS = -O3 -march=native -mavx -Werror -Wfatal-errors -Wall -std=c++0x #-fprofile-use #-pg -fprofile-generate
CXXFLAGS = -O3 -march=core-avx2 -Wfatal-errors -Wall -std=c++0x #-fprofile-use #-pg -fprofile-generate
CXXEFLAGS = -Wextra 
#CXX = icc
CXX = g++
GUILIBS = $(shell pkg-config --libs gtkmm-2.4 gtkglextmm-x11-1.2 libpng)
GUIFLAGS = $(shell pkg-config --cflags gtkmm-2.4 gtkglextmm-x11-1.2) -I.
CPPFLAGS = -DG_DISABLE_DEPRECATED -DGDK_PIXBUF_DISABLE_DEPRECATED -DPNG_SKIP_SETJMP_CHECK
OBJECTS=${SRC:=.o}
SPATIOCYTE_CORE = spatiocyte-core
VISUALIZER = visualizer


all: $(SPATIOCYTE_CORE) $(VISUALIZER)

$(SPATIOCYTE_CORE): $(OBJECTS)
		$(CXX) $(CXXFLAGS) $(CXXEFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

$(VISUALIZER):
		$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(IFLAGS) $(GUIFLAGS) -o $@ $(VIS).cpp $(GUILIBS)

%.o: %.cpp
		$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o $@ $<

clean:
		rm -f $(SPATIOCYTE_CORE) $(VISUALIZER) *.o
