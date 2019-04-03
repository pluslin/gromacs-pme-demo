CC = g++

TARGET = example
SRCS = gromacs/ewald/pme-spread-new.cpp \
	   gromacs/ewald/pme-malloc-free.cpp \
	   example.cpp

OBJS = $(SRCS:.cpp=.o)

$(TARGET):$(OBJS)
	$(CC) $^ -o $@

%.o:%.cpp
	$(CC) -c $^ -o $@

clean:
	rm -rf $(TARGET) $(OBJS)

run:
	./example
