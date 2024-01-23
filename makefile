EXE = FGM_PC_to_RL
SRC = lib/FGMlib/src/FGMlib.c lib/KdTreelib/src/KdTreelib.c lib/KNNlookuplib/src/KNNlookuplib.c lib/rectilinearMeshlib/src/rectilinearMeshlib.c src/main.c

OBJ = $(SRC:.c=.o)

# Update the include path to point to your new include directories
CFLAGS = -fPIC -fopenmp -Iinclude -Ilib/FGMlib/include -Ilib/KdTreelib/include -Ilib/KNNlookuplib/include -Ilib/rectilinearMeshlib/include
LINKFLAGS = -lm -fopenmp

# Compiler
CC = gcc

$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $(EXE) $(LINKFLAGS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJ) $(EXE)

