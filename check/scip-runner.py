#! /usr/bin/env python

# This script allows to run PySCIPOpt using the SCIP make test environment.
# usage: make test[cluster] EXECUTABLE=check/pyscipopt-runner.py [options]
#
# You may use this file as template to extend it and add custom code.
# The executable needs to have the .py file extension to be recognized by
# the check scripts.

from pyscipopt import Model
import os
import sys

def parse_tmp_file(tmpfile):
   """helper method to parse TMPFILE that has been set up by configuration_tmp_setup_scip.sh"""
   instance = ""
   settings = ""
   outsettings = ""
   params = []
   parammap = {
      "misc referencevalue" : ("misc/referencevalue", "f"),
      "limits time" : ("limits/time", "f"),
      "limits nodes" : ("limits/nodes", "lld"),
      "limits memory" : ("limits/memory", "f"),
      "lp advanced" : ("lp/threads", "d"),
      "timing clocktype" : ("timing/clocktype", "d"),
      "display freq" : ("display/freq", "d"),
      "memory savefac" : ("memory/savefac", "f")
    }

   for l in tmpfile:
      splitstr = l.split()
      # instance
      if l.startswith("read"):
         instance = splitstr[1]
      # input settings
      elif l.startswith("set load"):
         settings = splitstr[2]
      # output settings
      elif l.startswith("set save"):
         outsettings = splitstr[2]
      # SCIP parameters
      elif l.startswith("set"):
         p = splitstr[1] + " " + splitstr[2]
         assert p in parammap
         img = parammap[p]
         val = splitstr[-1]
         # store paramets as triples of the form (name,type,value)
         params.append((img[0], img[1], val))

   return {"instance" : instance, "settings" : settings, "outsettings" : outsettings, "params" : params}


def prepare_model(model, data):
   """prepares a PySCIPOpt model for given input data that has been parsed by a TMPFILE"""

   # read settings file
   model.readParams(data["settings"])

   # set parameters
   for (name,t,val) in data["params"]:
      if t == "f":
         model.setRealParam(name, float(val))
      elif t in "d":
         model.setIntParam(name, int(val))
      elif t == "lld":
         model.setLongintParam(name, int(val))

   # read problem
   model.readProblem(data["instance"])

#
# example main
#
if __name__ == "__main__":

   # generate a model
   model = Model()

   # parse TMP file
   data = parse_tmp_file(sys.stdin)

   # prepare model
   prepare_model(model, data)

   # print version, optimize, and print statistics
   model.printVersion()

   # write parameters
   model.writeParams(data["outsettings"])

   # optimize and print statistics
   model.optimize()
   model.printStatistics()
