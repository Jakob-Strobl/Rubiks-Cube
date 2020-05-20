# Using clang for compiling
CC = clang

HEAD = matrix.h initShader.h shapes.h camera.h solve_rc.h
OBJ = matrix.o initShader.o main.o shapes.o camera.o solve_rc.o
EXC = main # The name of executable

# Libraries depending on OS
MACLIBS = -L/System/Library/Frameworks -framework GLUT -framework OpenGL
LINUXLIBS = -lglut -lGLEW -lGLU -lm -lGL -lpthread

default:
	@echo ""
	@echo "WARNING: Just calling 'make' is not supported."
	@echo ""
	@echo "Please type:"
	@echo "  'make mac' for MacOS devices,"
	@echo "  'make linux' for linux devices."
	@echo ""

# Create C object files from source code.
%.o: %.c $(HEAD)
	$(CC) -c -o $@ $<

# MacOS Make configs - Tested on MacOS10.14
# Test file for assignment req
mac: $(OBJ)
	$(CC) -o $(EXC) $^ $(MACLIBS)

# Unix Make Configs - I use OpenSUSE so mileage may vary.
linux: $(OBJ)
	$(CC) -o $(EXC) $^ $(LINUXLIBS)


# Clean config
.PHONY: clean
clean:
	rm -f $(OBJ) $(TESTOBJ) $(UNITOBJ) $(TEST) $(UNITTEST) $(EXC)