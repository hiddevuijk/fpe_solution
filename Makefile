TARGET = test.exe
OBJS = main.o  system.o system2d.o xy.o
CC = g++
#CFLAGS = -c -Wall -g -std=c++11
#LFLAGS = -Wall -g
CFLAGS = -c -Wall -O3 -DNDEBUG -std=c++11
LFLAGS = -Wall  -O3 -DNDEBUG

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS)  $(OBJS) -o $(TARGET)

xy.o: xy.cpp xy.h
	$(CC) $(CFLAGS) xy.cpp

system.o: system.cpp system.h 
	$(CC) $(CFLAGS) system.cpp

system2d.o: system2d.cpp system2d.h 
	$(CC) $(CFLAGS) system2d.cpp

main.o: main.cpp  system.h system2d.h
	$(CC) $(CFLAGS) main.cpp


.PHONY: clean
clean:
	rm -f  $(OBJS) $(TARGET) 

.PHONY: cleanObject
cleanObject:
	rm -f  $(OBJS)

