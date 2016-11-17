#programs,flags,etc.
INCLUDE	=	-I include
TARGET = paraGSpan
SOURCE = src/main.cpp src/Graph.cpp src/Traveler.cpp src/SubTraveler.cpp src/Methods.cpp src/Global.cpp

#ALL Phony Targets
.PHONY : everything clean all

#Default starting position
everything:	$(TARGET)
clean:	
	rm -f $(TARGET)
all:	clean everything
	
paraGSpan:	src/main.cpp	\
			src/Graph.cpp include/Graph.h	\
			src/Traveler.cpp include/Traveler.h	\
			src/SubTraveler.cpp include/SubTraveler.h	\
			src/Methods.cpp include/Methods.h	\
			src/Global.cpp include/Global.h	\
			include/EdgeFrequency.h
	mpicxx $(INCLUDE) -w -g -fopenmp -o paraGSpan $(SOURCE) -O2