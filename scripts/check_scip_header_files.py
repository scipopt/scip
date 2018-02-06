#! /usr/bin/python
# -*- coding: utf-8 -*-
import re
import sys

'''
This script checks for correct synchronization between a
header file and its corresponding source file.

Invokation via

./check_scip_header_files.py [HEADER_FILE.h ...]

where HEADER_FILE.h is the (relative or absolute)
path to a header file, e.g., ../src/scip/scip.h.
Optionally, more than one header file may be passed.

'''

docstartstring = "/\*\*"
docendstring = "\*/"
docstartexpr = re.compile(docstartstring)
docendexpr = re.compile(docendstring)
docexpr = re.compile("/\*\*(.*|(?<!\*/)\n)*\*/", re.MULTILINE)
funcexpr = re.compile("(extern|EXTERN)+\n+(([\S_]+)\s*([a-zA-Z0-9]+)\(\n(.*|(?<!;)\n)*\));",re.MULTILINE) #", re.MULTILINE)

# loop over every header file
for filename in sys.argv[1:]:
    if not filename.endswith(".h"):
        print "  WARNING:", filename, " is not a header file (.h!!!), skipping it"
        continue
    print
    print "  --checking ", filename

    if filename.find("misc_") >= 0:
        print "    --> skipped"
        continue

    print
    print

    # first round: check for comments
    with open(filename, 'r') as currentfile:

        for match in docexpr.finditer(currentfile.read()):
            cfilename = filename.replace("pub_","").replace(".h", ".c")
            with open(cfilename, 'r') as cfile:

                if not match.group(1) in cfile.read():
                    print "WARNING: The documentation string"
                    print match.group(1)
                    print "was not found in %s"%cfilename
                    print

    #second round: check for function signatures
    with open(filename, 'r') as currentfile:
        for match in funcexpr.finditer(currentfile.read()):
            cfilename = filename.replace("pub_","").replace(".h", ".c")
            with open(cfilename, 'r') as cfile:
                if not match.group(2) in cfile.read():
                    print "WARNING: The function signature "
                    print match.group(2)
                    print "was not found in %s"%cfilename
                    print

    print "  --finished", filename

  #indoc = False
  #lines = []
  #for line in currentfile:
    #startfound = docstartexpr.search(line)
    #if startfound:
      #indoc = True
    #if indoc:
      #lines.append(line)

    #endfound = docendexpr.search(line)
    #if endfound and indoc:
      #for theline in lines:
	#print theline
      #indoc = False
      #lines = []
