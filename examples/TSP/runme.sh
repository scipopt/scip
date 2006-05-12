#!/bin/sh
make
cd tspviewer
javac TSPViewer.java
java TSPViewer &
cd ..
sleep 2 # wait for the TSPViewer
bin/sciptsp.linux.x86.gnu.opt.static.cpx -f tspdata/pr76.tsp
