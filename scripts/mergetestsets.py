#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import commands

instances = []
paths = {}

def main():
    if len(sys.argv) < 4:
        print 'usage: mergetestsets.py <.test file 1> <.test file 2> <location of merged .test file>'
        return 0

    infile1 = file(sys.argv[1], "r")
    infile2 = file(sys.argv[2], "r")
    outfile = file(sys.argv[3], "w")
    solufile = file(sys.argv[3].replace(".test","")+".solu", "w")

    for line in infile1.readlines() + infile2.readlines():
        if line.startswith("#"):
            continue
        else:
            instance = line.split("/")[-1].strip().replace(".mps","").replace(".lp","").replace(".gz","")

            if not instance in instances:
                instances.append(instance)
                paths[instance] = line.strip()
            else:
                if paths[instance] != line.strip():
                    print "instance %s with two different paths:\n\t%s\n\t%s" %(instance, paths[instance], line.strip())

    for instance in sorted(instances):
        outfile.write("%s\n" %(paths[instance]))

        grepline = commands.getoutput("grep '%s ' %s %s" %(instance, sys.argv[1].replace(".test","")+".solu", sys.argv[2].replace(".test","")+".solu"))

        if grepline.strip() == "":
            grepline = commands.getoutput("grep '%s$' %s %s" %(instance, sys.argv[1].replace(".test","")+".solu", sys.argv[2].replace(".test","")+".solu"))

        line1 = grepline.strip().split("\n")[0].split(":",1)[1]

        if len(grepline.strip().split("\n")) != 1:
            assert len(grepline.strip().split("\n")) == 2, grepline
            line2 = grepline.strip().split("\n")[1].split(":",1)[1]

            if line1.split()[0] != line2.split()[0]\
               or (line1.split()[0] == "=opt=" and \
                   abs(float(line1.split()[2]) - float(line2.split()[2])) > 1e-04):
                print "different solution file entries:\n" + grepline

        solufile.write(line1+"\n")

if __name__ == "__main__":
        sys.exit( main() )
