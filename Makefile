TARGET = opg

OBJECTS = main.o Matrix.o algorithm.o

CC = g++

CFLAGS = -c




$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -lm -lpthread -o $(TARGET) -O3
	
	
main.o: main.cpp Matrix.cpp algorithm.cpp
	$(CC) $(CFLAGS) main.cpp -O3
	
	
Matrix.o: Matrix.cpp
	$(CC) $(CFLAGS) Matrix.cpp -O3
	
	
algorithm.o: algorithm.cpp
	$(CC) $(CFLAGS) algorithm.cpp -O3
	
	
.PHONY: clean
clean:
	-rm -f $(OBJECTS) $(TARGET)
