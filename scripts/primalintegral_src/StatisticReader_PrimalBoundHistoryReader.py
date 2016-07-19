from StatisticReader import StatisticReader
import re

class PrimalBoundHistoryReader(StatisticReader):
   name = 'PrimalBoundHistoryReader'
   regular_exp =  ' time | node  | left  |LP iter|LP it/n| mem |mdpt |frac |vars |cons |cols |rows |cuts |confs|strbr|  dualbound   | primalbound  |  gap   '
   columnwidth = 20
   datakey = 'PrimalBoundHistory'
   columnheaderstr = 'PrimalBoundHistory'.rjust(columnwidth)
   inTable = False
   lastPrimalBound = '--'
   listOfPoints = []
   easyCPLEX= True
   totalnumberofsols = 0
       
   testnumberofsols = 0
   gurobiextralist = []

   def extractStatistic(self, line):
      if StatisticReader.solvertype == StatisticReader.SOLVERTYPE_SCIP:
         if self.regular_exp in line:
            self.inTable = True
         elif self.inTable and re.search("s\|", line):
   #         print line
            mylist = line.split()
   #         print mylist
            if len(mylist) == 0:
               return None
            if re.search('s\|', mylist[0]):
               rawtime = mylist[0].split('s|')[0]
            else:
               assert(re.search('s|', mylist[1]))
               rawtime = mylist[1].split('s|')[0]
                     
            rawtime = re.split('[a-zA-Z*]', rawtime)
                  
            pointInTime = rawtime[len(rawtime) - 1]
            if re.search("\|",mylist[len(mylist) - 1]):
               PrimalBound = mylist[len(mylist) - 2].strip("|")
            else:
               PrimalBound = mylist[len(mylist) - 3].strip("|")
   #         print pointInTime, PrimalBound
            if PrimalBound != self.lastPrimalBound:
               self.lastPrimalBound = PrimalBound
               try:
                   self.listOfPoints.append((float(pointInTime),float(PrimalBound)))
               except ValueError:
                   pass

         elif "SCIP Status" in line and self.inTable:
            self.inTable = False

      elif StatisticReader.solvertype == StatisticReader.SOLVERTYPE_GUROBI:
         if "Found heuristic solution" in line:
            self.gurobiextralist.append(line.split()[-1])
         if "Expl Unexpl |  Obj  Depth" in line:
            self.inTable = True
         elif self.inTable and line.endswith("s\n") and self.gurobiextralist != []:
#            print "+++++++++++++++++++++"
            pointInTime = line.split()[-1].strip("s")
            self.listOfPoints.append((float(pointInTime), float(self.gurobiextralist[-1])))
            self.gurobiextralist = []
         elif self.inTable and line.startswith("H") or line.startswith("*"):
            self.readBoundAndTime(line, -5, -1, timestripchars="s")
             
         elif "Cutting planes:" in line and self.inTable:
            self.inTable = False
         elif self.gurobiextralist != [] and "Explored " in line:
#            print "-------------------------"
            pointInTime = line.split()[-2]
            self.listOfPoints.append((float(pointInTime), float(self.gurobiextralist[-1])))
            self.gurobiextralist = []
          
         return None
            
      elif StatisticReader.solvertype == StatisticReader.SOLVERTYPE_CBC:
         if "Integer solution of " in line:
            self.readBoundAndTime(line, 4, -2, timestripchars="(")
             
         return None
          
      elif StatisticReader.solvertype == StatisticReader.SOLVERTYPE_XPRESS:
         if "BestSoln" in line:
            self.xpresscutidx = line.index("BestSoln") + len("BestSoln")
         elif line.startswith("+") or line.startswith("*"):
            self.readBoundAndTime(line, -1, -1, cutidx=self.xpresscutidx)
         elif line.startswith(" *** Heuristic solution found: "):
            self.readBoundAndTime(line, -4, -2)
             
      elif StatisticReader.solvertype == StatisticReader.SOLVERTYPE_CPLEX:
         if "Solution pool: " in line:
            self.testnumberofsols = int(line.split()[2])
            # assert len(self.listOfPoints) >= self.testnumberofsols
         if self.easyCPLEX and "Found incumbent of value" in line:
            splitline = line.split()
            self.readBoundAndTime(line, splitline.index("Found") + 4, splitline.index("Found") + 6)
         elif not self.easyCPLEX:
            if "Welcome to IBM(R) ILOG(R) CPLEX(R)" in line:
               self.lastelapsedtime = 0.0
               self.nnodessincelastelapsedtime = 0
               self.lastnnodes = 0
               self.cpxprimals = []
                
            if "   Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap" in line:
               self.inTable = True

               self.cpxbestintegeridx = line.index("Best Integer") + 11
                
            elif self.inTable and ("cuts applied:" in line or "Root node processing" in line):
               self.inTable = False
            elif self.inTable and "Repeating presolve." in line:
               self.inTable = False
            elif self.inTable and "Covers:" in line:
               self.inTable = False
            elif self.inTable and len(line) > 0 and line.startswith(" ") or line.startswith("*"):
               if line=="\n":
                  return None
               nodeinlineidx = 7
               while line[nodeinlineidx] != " " and line[nodeinlineidx] != "+":
                  nodeinlineidx += 1
               nnodes = int(line[:nodeinlineidx].split()[-1].strip('*+')) + 1
               if line.startswith("*") or line.startswith("+"):
                  print line
                  primalbound = float(line.split()[-4])
                  print primalbound, nnodes
                  self.cpxprimals.append((nnodes, primalbound))
               self.lastnnodes = nnodes
            elif "Elapsed time = " in line:
               thetime = float(line.split()[3])
               self.processCpxprimals(thetime)
                
               self.nnodessincelastelapsedtime = self.lastnnodes
               self.lastelapsedtime = thetime

            elif "Solution time =" in line:
               thetime = float(line.split()[3])
               self.processCpxprimals(thetime)
       
      return None
       
   def readBoundAndTime(self, line, boundidx, timeidx, timestripchars="", cutidx=-1):
      splitline = line.split()
        
      primalbound = line[:cutidx].split()[boundidx]
        
      pointInTime = splitline[timeidx].strip(timestripchars)
             
#        print line 
#        print pointInTime, primalbound
        
        
      self.lastPrimalBound = primalbound
      self.listOfPoints.append((float(pointInTime), float(primalbound)))
       
   def processCpxprimals(self, currenttime):
      solvednodes = (self.lastnnodes - self.nnodessincelastelapsedtime)
      timespentonnode = (currenttime - self.lastelapsedtime)/max(1.0, float(solvednodes))
      assert currenttime >= self.lastelapsedtime
      for node, bound in self.cpxprimals:
         estimatedtime = self.lastelapsedtime + (node - self.nnodessincelastelapsedtime) * timespentonnode
                   
         if bound != self.lastPrimalBound:
            self.lastPrimalBound = bound
            if len(self.listOfPoints) > 0:
               othertime, otherbound = self.listOfPoints[-1]
               assert othertime <= estimatedtime
               self.listOfPoints.append((float(estimatedtime), float(bound)))
               
      self.cpxprimals = []

   def execEndOfProb(self):
      self.inTable = False
      self.testnumberofsols = 0
      theList = self.listOfPoints[:]
      self.listOfPoints = []
      self.lastPrimalBound = '--'
      PrimalBoundHistoryReader.totalnumberofsols += len(theList)
      #print " solutions added:", PrimalBoundHistoryReader.totalnumberofsols
      return theList