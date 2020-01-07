#! /usr/bin/env python2
# -*- coding: utf-8 -*-
import sys
import commands

'''
This script updates the SCIP version

Invokation via

./updateversion.py [-a] [-v <VERSION>]

where -a signals that an API change was done and the API version
number should be increased (this is the default when providing no arguments)
and -v <VERSION> set the SCIP version number to VERSION, e.g., 4.0.0, or 4.0.0.1.

'''

# search for parameters
vpos = -1
apos = -1
for i in range(1,len(sys.argv)):
    argument = sys.argv[i]
    if argument == "-a":
        apos = i
    if argument == "-v":
        vpos = i
if vpos == -1:
    apos = 1

# check correct usage
if len(sys.argv) == 1 or len(sys.argv) > 4 or vpos == len(sys.argv)-1 or vpos == apos -1:
    print "usage: scripts/updateversion.py [-a] [-v <MAJOR>.<MINOR>.<PATCH>[.<SUB>]]"
else:
    # old version numbers
    oldversion = commands.getoutput("grep SCIP_VERSION src/scip/def.h").split()[2]
    oldmajor = newmajor = oldversion[0]
    oldminor = newminor = oldversion[1]
    oldpatch = newpatch = oldversion[2]
    oldsubversion = newsubversion = commands.getoutput("grep SCIP_SUBVERSION src/scip/def.h").split()[2]
    if commands.getoutput("grep SCIP_APIVERSION src/scip/def.h") != "":
        oldapiversion = newapiversion = commands.getoutput("grep SCIP_APIVERSION src/scip/def.h").split()[2]
    else:
        oldapiversion = newapiversion = None
        if apos != -1:
            print "this SCIP version has no API version (yet), API version cannot be updated!"
            apos = -1

    if apos != -1:
        newapiversion = str(int(oldapiversion)+1)
        print "update SCIP API version from %s to %s" %(oldapiversion, newapiversion)

    if vpos != -1:
        tokens = sys.argv[vpos+1].split(".")
        assert len(tokens) == 3 or len(tokens) == 4, "The version number must consist of 3 or 4 numbers, divided by '.'s (is: "+sys.argv[vpos+1]+")"

        # new version numbers
        newmajor = tokens[0]
        newminor = tokens[1]
        newpatch = tokens[2]
        if len(tokens) == 4:
            newsubversion = tokens[3]
        else:
            newsubversion = "0"

    newversion = newmajor+newminor+newpatch

    oldversionstring = oldmajor + "." + oldminor + "." + oldpatch
    if oldsubversion != "0":
        oldversionstring += "." + oldsubversion

    newversionstring = newmajor + "." + newminor + "." + newpatch
    if newsubversion != "0":
        newversionstring += "." + newsubversion

    if newversionstring != oldversionstring:
        print "update SCIP version from %s to %s" %(oldversionstring, newversionstring)

    # update version numbers
    commands.getoutput('sed -i "s/\#define SCIP_VERSION.*%4s/\#define SCIP_VERSION               %4s/" src/scip/def.h' \
                       %(oldversion, newversion))
    commands.getoutput('sed -i "s/\#define SCIP_SUBVERSION.*%3s/\#define SCIP_SUBVERSION             %3s/" src/scip/def.h' \
                       %(oldsubversion, newsubversion))
    commands.getoutput('sed -i "s/\@version.*/\@version  %-s/" doc/xternal.c' %(newversionstring))
    commands.getoutput('sed -i "s/^SCIP_VERSION.*/SCIP_VERSION	=	%-s/" make/make.project' %(newversionstring))
    commands.getoutput('''sed -i 's/^VERSION=.*/VERSION=\"%s\"/' makedist.sh''' %(newversionstring))
    commands.getoutput('sed -i "s/set(SCIP_VERSION_MAJOR %s)/set(SCIP_VERSION_MAJOR %s)/" CMakeLists.txt' %(oldmajor, newmajor))
    commands.getoutput('sed -i "s/set(SCIP_VERSION_MINOR %s)/set(SCIP_VERSION_MINOR %s)/" CMakeLists.txt' %(oldminor, newminor))
    commands.getoutput('sed -i "s/set(SCIP_VERSION_PATCH %s)/set(SCIP_VERSION_PATCH %s)/" CMakeLists.txt' %(oldpatch, newpatch))
    commands.getoutput('sed -i "s/set(SCIP_VERSION_SUB %s)/set(SCIP_VERSION_SUB %s)/" CMakeLists.txt' %(oldsubversion, newsubversion))

    if newapiversion != oldapiversion:
        print("\nAPI versions before the change:")
        print(commands.getoutput('grep -e "APIVERSION" -e "VERSION_API" src/scip/def.h  CMakeLists.txt'))
        if newsubversion == "0":
            print "\nWarning: API version increased for what seems to be a bugfix version (%s)" %(newversionstring)
        commands.getoutput('sed -i "s/\#define SCIP_APIVERSION.*%3s/\#define SCIP_APIVERSION             %3s/" src/scip/def.h' \
                           %(oldapiversion, newapiversion))
        commands.getoutput('sed -i "s/set(SCIP_VERSION_API [0-9]*)/set(SCIP_VERSION_API %s)/" CMakeLists.txt' \
                           %(newapiversion))
        print("\nAPI versions after the change:")
        print(commands.getoutput('grep -e "APIVERSION" -e "VERSION_API" src/scip/def.h  CMakeLists.txt'))
