all: binary-binary binary-single v-dist bs-phase bb-phase

VPATH =${MY_REPO_PATH}/SpaceX:${MY_REPO_PATH}/SpaceHub/src/
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall  
OBJS = dynamicSystem.o InitFileCreator.o \
	 chain.o chainSystem.o tools.o\
	 IO.o symplectic.o regularizator.o BSIterator.o\
	particleSystem.o integrator.o seqIterator.o Wrapper.o
	
INCL = dynamicSystem.h macros.h chain.h\
	regularizator.h particleSystem.h interaction.h chainSystem.h \
	integrator.h BSIterator.h errhand.h InitFileCreator.h\
	seqIterator.h odeIterator.h vector3d.h tools.h

binary-binary : $(OBJS) binary-binary.o
	$(CXX) -std=c++17 -O3 -o spacehub.ex $(OBJS) binary-binary.o -pthread

binary-single : $(OBJS) binary-single.o
	$(CXX) -std=c++17 -O3 -o binary-single.ex $(OBJS) binary-single.o -pthread

bb-phase : $(OBJS) binary-binary-phase.o
	$(CXX) -std=c++17 -O3 -o bb-phase.ex $(OBJS) binary-binary-phase.o -pthread

bs-phase : $(OBJS) binary-single-phase.o
	$(CXX) -std=c++17 -O3 -o bs-phase.ex $(OBJS) binary-single-phase.o -pthread

v-dist: $(OBJS) v-dist.o
	$(CXX) -std=c++17 -O3 -o v-dist.ex $(OBJS) v-dist.o -pthread

$(OBJS) : $(INCL)
.PHONY : clean	
clean:
	-rm $(OBJS) *.o *.ex
