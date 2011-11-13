CUDA = $(subst /bin/nvcc,,$(shell echo `which nvcc`))
NVCC = $(CUDA)/bin/nvcc
AP64 = $(shell if [ `uname` = Linux ]; then echo 64; fi)

LDFLAGS = $(addprefix -Xlinker , -rpath $(CUDA)/lib)$(AP64)
CFLAGS  = $(addprefix --compiler-options , -Wall) -O3

compile:
	mkdir -p bin
	$(NVCC) src/*.{cu,cpp} $(CFLAGS) $(LDFLAGS) -o bin/hydro

clean:
	-rm -f bin/hydro
	-rm -f */*~
	-if [ -z "`ls bin 2>&1`" ]; then rmdir bin; fi
