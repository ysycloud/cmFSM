#programs,flags,etc.
INCLUDE	=	-I include
TARGET = cmFSM
SOURCE = src/main.cpp src/Graph.cpp src/Traveler.cpp src/SubTraveler.cpp src/Methods.cpp src/Global.cpp src/Supervisor.cpp src/IO.cpp

#ALL Phony Targets
.PHONY : everything clean all

#Default starting position
everything:	$(TARGET)
clean:	
	rm -f $(TARGET)
all:	clean everything
	
cmFSM:	src/main.cpp	\
			src/Graph.cpp include/Graph.h	\
			src/Traveler.cpp include/Traveler.h	\
			src/SubTraveler.cpp include/SubTraveler.h	\
			src/Methods.cpp include/Methods.h	\
			src/Global.cpp include/Global.h	\
			src/Supervisor.cpp include/Supervisor.h	\
			src/IO.cpp include/IO.h	\
			include/EdgeFrequency.h
	mpicxx $(INCLUDE) -g -fopenmp -o cmFSM $(SOURCE) -O2