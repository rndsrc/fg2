CUDA = $(subst /bin/nvcc,,$(shell echo `which nvcc`))
NVCC = $(CUDA)/bin/nvcc
AP64 = $(shell if [ `uname` = Linux ]; then echo 64; fi)

LDFLAGS = $(addprefix -Xlinker , -rpath $(CUDA)/lib)$(AP64)
CFLAGS  = $(addprefix --compiler-options , -Wall) -O3

help:
	@echo 'The follow system of equations is avilable:'
	@echo
	@c=0; for F in src/*.cuh; do  \
	   f=$${F##src/};             \
	   c=`expr $$c + 1`;          \
	   echo "  $$c. $${f%%.cuh}"; \
	 done
	@echo
	@echo 'Use `make NAME` and `bin/NAME` to compile and run fg2.'

%:
	@mkdir -p bin
	@echo -n 'Compiling... '
	@$(NVCC) src/*.{cu,cpp} $(CFLAGS) $(LDFLAGS) -DEQNS=$@.cuh -o bin/$@
	@if [ -f bin/$@ ]; then            \
	   echo 'DONE';                    \
	   echo 'Use `bin/$@` to run fg2'; \
	 else                              \
	   echo 'FAIL!!!';                 \
	 fi

clean:
	-rm -f bin/hydro
	-rm -f */*~
	@if [ -z "`ls bin 2>&1`" ]; then rmdir bin; fi
