#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

# This script prints all instancens from the given testset file
# that are missing in the given .res file.

import sys
import os
import math

if len(sys.argv) != 3:
    print '"Usage:'
    print sys.argv[0]+' testset resfile'
    quit()

testset = open(sys.argv[1], 'r')
resfile = open(sys.argv[2], 'r')
reslines = resfile.readlines()
resfile.close()

namelength = 18

for line in testset:
   linesplit = line.split('/')
   linesplit = linesplit[len(linesplit) - 1].rstrip(' \n').rstrip('.gz').rstrip('.GZ').rstrip('.z').rstrip('.Z')
   linesplit = linesplit.split('.')
   instancename = linesplit[0]
   for i in range(1, len(linesplit)-1):
      instancename = instancename + '.' + linesplit[i]
   length = len(instancename)
   if length > namelength:
      instancename = instancename[length-namelength-2:length-2]
   missing = True
   for result in reslines:
      if result.startswith(instancename):
         missing = False
   if missing:
      print line.rstrip(' \n')

quit()