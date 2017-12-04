#!/usr/bin/env python3
import sys
import re
import os
import tempfile
from gurobipy import *


if len(sys.argv) < 2:
   print('Usage: gurobi.sh checkcuts.py filename [cutnamefilter]')
   quit()

filename = sys.argv[1]
cutclassname = ""
if len(sys.argv) == 3:
   cutclassname = sys.argv[2]

base=os.path.basename(filename)
prefix=os.path.splitext(base)[0]

basemodelname=prefix + "_withoutcuts.lp"
lazyconss = None

with open (filename, "r") as myfile:
   data=myfile.readlines()
   lazyconsstart = data.index("lazy constraints\n")
   boundstart = data.index("Bounds\n")
   lines1 = data[:lazyconsstart]
   lines2 = data[boundstart:]
   lazyconss = data[lazyconsstart+1:boundstart]

   with open(basemodelname, "w") as outfile:
      for line in lines1:
         outfile.write(line)

      for line in lines2:
         outfile.write(line)

model = read(basemodelname)
os.remove(basemodelname)
model.setParam("OutputFlag", 0);
model.setParam("FeasibilityTol", 1e-9)

validcuts = 0
invalidcuts = 0
currcons = ""
for lazycons in lazyconss:
   currcons += lazycons
   m = re.match("([^:]*):([^<]*)(<=|>=)([^\n]*)", currcons)
   if m is None:
      currcons = currcons[:-1]
      continue
   currcons=""
   consname = m.group(1).strip()
   if not cutclassname in consname:
      continue
   cons = m.group(2).strip()
   rhs = float(m.group(4).strip())

   if m.group(3) == "<=":
      objsense = GRB.MAXIMIZE
   else:
      objsense = GRB.MINIMIZE

   objexpr = LinExpr()
   while len(cons) > 0:
      m = re.match("[\+-][0-9][^ ]*[^\+-]*", cons)
      length = len(m.group(0))
      term = m.group(0).strip()
      coef,varname = term.split()
      cons = cons[length:]
      var = model.getVarByName(varname)
      objexpr.add(var, float(coef))

   model.setObjective(objexpr, objsense)
   model.optimize()
   if model.status == GRB.Status.INF_OR_UNBD:
      print("\nError: model is infeasible or unbounded")
      quit()

   viol = model.objVal - rhs

   if viol > 1e-9:
      solname = prefix + "_" + consname + ".sol"
      modelfilename = prefix + "_" + consname + ".lp"
      model.write(solname)
      model.write(modelfilename)
      print("cut " + consname + " is invalid!")
      print("\tmodel to verify cut has been written to " + modelfilename)
      print("\tinteger solution that violates the cut by "+ str(viol) + " was written to " + solname)
      invalidcuts += 1
   else:
      validcuts += 1

print("found "+str(invalidcuts)+" invalid and "+str(validcuts) + " valid cuts")
