#!/bin/bash
## output in current dir, source in parent dir. 
sphinx-apidoc -o . ..

make clean html