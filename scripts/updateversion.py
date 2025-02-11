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
and -v <VERSION> set the SCIP version number to VERSION, e.g., 4.0.0.

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
    print("usage: scripts/updateversion.py [-a] [-v <MAJOR>.<MINOR>.<PATCH>]")
else:
    # old version numbers
    oldmajor = newmajor = commands.getoutput("grep 'SCIP_VERSION_MAJOR =' make/make.project").split()[2]
    oldminor = newminor = commands.getoutput("grep 'SCIP_VERSION_MINOR =' make/make.project").split()[2]
    oldpatch = newpatch = commands.getoutput("grep 'SCIP_VERSION_PATCH =' make/make.project").split()[2]
    oldapiversion = newapiversion = commands.getoutput("grep 'SCIP_VERSION_API =' make/make.project").split()[2]

    if apos != -1:
        newapiversion = str(int(oldapiversion)+1)
        print("update SCIP API version from %s to %s" %(oldapiversion, newapiversion))

    if vpos != -1:
        tokens = sys.argv[vpos+1].split(".")
        assert len(tokens) == 3 or len(tokens) == 4, "The version number must consist of 3 or 4 numbers, divided by '.'s (is: "+sys.argv[vpos+1]+")"

        # new version numbers
        newmajor = tokens[0]
        newminor = tokens[1]
        newpatch = tokens[2]

    newversion = newmajor+newminor+newpatch

    oldversionstring = oldmajor + "." + oldminor + "." + oldpatch
    newversionstring = newmajor + "." + newminor + "." + newpatch

    if newversionstring != oldversionstring:
        print("update SCIP version from %s to %s" %(oldversionstring, newversionstring))

    # update version numbers
    commands.getoutput('sed -i "s/\@version.*/\@version  %-s/" doc/xternal.c' %(newversionstring))
    commands.getoutput('''sed -i 's/^VERSION=.*/VERSION=\"%s\"/' scripts/makedist.sh''' %(newversionstring))
    commands.getoutput('sed -i "s/^SCIP_VERSION_MAJOR.*/SCIP_VERSION_MAJOR = %-s/" make/make.project' % newmajor)
    commands.getoutput('sed -i "s/^SCIP_VERSION_MINOR.*/SCIP_VERSION_MINOR = %-s/" make/make.project' % newminor)
    commands.getoutput('sed -i "s/^SCIP_VERSION_PATCH.*/SCIP_VERSION_PATCH = %-s/" make/make.project' % newpatch)
    commands.getoutput('sed -i "s/set(SCIP_VERSION_MAJOR %s)/set(SCIP_VERSION_MAJOR %s)/" CMakeLists.txt' %(oldmajor, newmajor))
    commands.getoutput('sed -i "s/set(SCIP_VERSION_MINOR %s)/set(SCIP_VERSION_MINOR %s)/" CMakeLists.txt' %(oldminor, newminor))
    commands.getoutput('sed -i "s/set(SCIP_VERSION_PATCH %s)/set(SCIP_VERSION_PATCH %s)/" CMakeLists.txt' %(oldpatch, newpatch))

    if newapiversion != oldapiversion:
        print("\nAPI versions before the change:")
        print(commands.getoutput('grep -e "APIVERSION" -e "SCIP_VERSION_API" make/make.project CMakeLists.txt'))
        if newminor != "0" or newpatch != "0":
            print("\nWarning: API version increased for what does not seem to be the master branch (version %s)" %(newversionstring))
        commands.getoutput('sed -i "s/^SCIP_VERSION_API.*/SCIP_VERSION_API = %-s/" make/make.project' % newapiversion)
        commands.getoutput('sed -i "s/set(SCIP_VERSION_API [0-9]*)/set(SCIP_VERSION_API %s)/" CMakeLists.txt' \
                           %(newapiversion))
        print("\nAPI versions after the change:")
        print(commands.getoutput('grep -e "APIVERSION" -e "SCIP_VERSION_API" make/make.project CMakeLists.txt'))
