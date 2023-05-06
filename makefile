SHELL := /bin/bash

n ?= 5
file ?= matrix.txt

all: compile run clean

compile:
	/opt/nfs/config/station204_name_list.sh 1 16 >nodes && \
	/opt/nfs/mpich-4.0.1/bin/mpicxx QR.cpp -o QR.out

run:
	/opt/nfs/mpich-4.0.1/bin/mpiexec -f nodes -n $(n) ./QR.out $(file)

clean:
	rm -f QR; rm nodess