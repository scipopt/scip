#! /usr/bin/env python

# This script allows to run PySCIPOpt using the SCIP make test environment.
# usage: make test[cluster] EXECUTABLE=check/pyscipopt-runner.py PYTHON=python [options]
#
# You may use this file as template to extend it and add custom code.
# The executable needs to have the .py file extension to be recognized by
# the check scripts and _must_ be an executable file.
#
# It should _not_ be necessary to modify parse_tmp_file() nor build_model(),
# unless you need to add a plugin before the problem instance is read.

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
         # store parameters as triples of the form (name,type,value)
         params.append((img[0], img[1], val))

   return {"instance" : instance, "settings" : settings, "outsettings" : outsettings, "params" : params}


def build_model():
   """builds a model, reads input data provided by the scripts, and writes out setting file"""
   # generate a model
   model = Model()

   # parse TMP file
   data = parse_tmp_file(sys.stdin)

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

   # prepare model
   prepare_model(model, data)

   # print version, optimize, and print statistics
   model.printVersion()

   # write parameters
   model.writeParams(data["outsettings"])

   return model

#
# example main
#
if __name__ == "__main__":

   model = build_model()

   # optimize and print statistics
   model.optimize()
   model.printStatistics()
