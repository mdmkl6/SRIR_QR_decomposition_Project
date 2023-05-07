SHELL := /bin/bash

file ?= matrix.txt
outfile ?= 

all: compile run clean

compile:
	/opt/nfs/config/station204_name_list.sh 1 16 >nodes; \
	/opt/nfs/mpich-4.0.1/bin/mpicxx QR.cpp -o QR.out

run:
	n=$$(wc -l nodes | awk '{print $$1}'); \
	/opt/nfs/mpich-4.0.1/bin/mpiexec -f nodes -n $$n ./QR.out $(file) $(outfile)

clean:
	rm -f QR.out; rm nodes