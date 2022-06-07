CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=choose.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=choose

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
