#!/bin/bash

make;
./simpar 1 3 10 1;
./simpar 1 3 1000000 20;
./simpar 1 10 2000000 10;
./simpar 1 30 20000000 10;