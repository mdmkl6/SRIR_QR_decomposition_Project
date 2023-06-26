SHELL := /bin/bash

file ?= ../matrix.txt
outfile ?= ../output.txt

all: compile run

compile:
	/opt/nfs/config/station204_name_list.sh 1 16 >nodes; \
  /opt/nfs/berkeley_upcxx-2023.3.0/bin/upcxx -O2 QR.cpp utils.h -o QR.out
 
run:
	/opt/nfs/berkeley_upcxx-2023.3.0/bin/upcxx-run -shared-heap 256M -n 4 $(shell upcxx-nodes nodes) QR.out $(file) $(outfile)

clean:
	rm -f QR.out; rm nodes

all-test: compile-test run-test clean-test

compile-test:
	g++ ../test.cpp -o ../test.out

run-test:
	./../test.out $(file) $(outfile)

clean-test:
	rm -f ./../test.out