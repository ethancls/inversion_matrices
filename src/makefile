CC = gcc
CFLAGS = -O3

all: projet.exe courageux.exe mult.exe

projet.exe: projet.c
	$(CC) $(CFLAGS) -o projet.exe projet.c operations.c

courageux.exe: courageux.c
	$(CC) $(CFLAGS) -o courageux.exe courageux.c operations.c

mult.exe: mult.c
	$(CC) $(CFLAGS) -o mult.exe mult.c operations.c
clean:
	rm -f courageux.exe projet.exe mult.exe