FLAGS=-pedantic -Wall -Werror -Wno-sign-compare -Wno-long-long -lm -std=c++11
COMPILLER=g++
LDFLAGS=-m64
LIBS=-lm
SOURCES=nm-lab5.cpp
EXECUTABLE=$(SOURCES:%.cpp=%)

TARGET=Debug
ifeq ($(TARGET), Debug)
	#FLAGS+=-Og -g
	FLAGS+=-g
else ifeq ($(TARGET), Release)
	FLAGS+=-fomit-frame-pointer -fexpensive-optimizations -flto -O3
	LDFLAGS+=-O3 -flto -s
endif

all: start

start: nm-lab5.cpp
	$(COMPILLER) $(FLAGS) $< -o nm-lab5

clean:
	rm -frd nm-lab5 output_5.csv error_5.csv
