# Compiler and Flags
NVCC = nvcc

# Program and Output
PROGRAM = B210473CS_A2.cu
OUTPUT = program
CFLAGS = -O3

# Input and Output Files
INPUTS = indata50.txt indata80.txt indata100.txt indata200.txt indata500.txt
OUTPUTS = output50.txt output80.txt output100.txt output200.txt output500.txt
REPORTS= report50.txt report80.txt report100.txt report200.txt report500.txt

# Default target
.PHONY: all build run clean

all: build run

# Build the program
build:
	$(NVCC) $(CFLAGS) $(PROGRAM) -o $(OUTPUT)

# Run the program with each input and output file
run: build
	@echo "Running program with indata50.txt..."
	./$(OUTPUT) indata50.txt output50.txt report50.txt
	@echo "Running program with indata80.txt..."
	./$(OUTPUT) indata80.txt output80.txt report80.txt
	@echo "Running program with indata100.txt..."
	./$(OUTPUT) indata100.txt output100.txt report100.txt
	@echo "Running program with indata200.txt..."
	./$(OUTPUT) indata200.txt output200.txt report200.txt
	@echo "Running program with indata500.txt..."
	./$(OUTPUT) indata500.txt output500.txt report500.txt

# Clean the generated files
exec: build
		./$(OUTPUT) input.txt output.txt report.txt
clean:
	rm -f $(OUTPUT) $(OUTPUTS) $(REPORTS)
