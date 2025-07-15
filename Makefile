CC = gcc
CFLAGS = -Wall -O3 -flto -ffast-math -march=native `sdl2-config --cflags`
LDFLAGS = `sdl2-config --libs` -lSDL2_image -lm -lasound -lfftw3f
TARGET = VibeScope
SRC = VibeScope.c incbin.h

all: $(TARGET)

$(TARGET): VibeScope.c incbin.h
	$(CC) $(CFLAGS) -o $@ VibeScope.c $(LDFLAGS)

clean:
	rm -f $(TARGET)
