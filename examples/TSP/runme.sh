#!/bin/sh
make
cd tspviewer
javac TSPViewer.java
java TSPViewer &
cd ..
bin/sciptsp.linux.x86.gnu.opt.spx -f tspdata/pr76.tsp
