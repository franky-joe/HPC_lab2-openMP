# Nombre del compilador
CC = gcc

# Flags de compilación
CFLAGS = -Wall -O2 -fopenmp 

# Librerías necesarias (cfitsio, matemática)
LIBS = -lcfitsio -lm

# Archivo de salida
TARGET = lab2

# Nombre del archivo fuente
SRC = main.c

# Regla para compilar el programa
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LIBS)

# Regla para limpiar los archivos compilados
clean:
	rm -f $(TARGET)
