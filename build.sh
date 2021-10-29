#!/bin/bash

gcc src/grogue.c -O2 -Wall -o bin/grogue -pthread -lm -lgsl -lmpfr -lgmp -lgslcblas -lexpat
