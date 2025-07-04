CC = gcc
CFLAGS = -Wall -O0 `sdl2-config --cflags`
LDFLAGS = `sdl2-config --libs` -lSDL2_image -lm
TARGET = VibeScope
SRC = VibeScope.c incbin.h

all: $(TARGET)

$(TARGET): VibeScope.c incbin.h
	$(CC) $(CFLAGS) -o $@ VibeScope.c $(LDFLAGS)

clean:
	rm -f $(TARGET)
