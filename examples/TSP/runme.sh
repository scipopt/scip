#!/bin/sh
make
cd tspviewer
javac TSPViewer.java
java TSPViewer &
cd ..
bin/sciptsp -f tspdata/pr76.tsp
