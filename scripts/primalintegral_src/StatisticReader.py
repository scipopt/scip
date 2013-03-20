# -*- coding: utf-8 -*-


class StatisticReader:
   '''
   base class for all statistic readers - readers should always inherit from this base class for minimal implementation
   effort
   
   readers only need to overwrite the methods extractStatistic() and perhaps execEndOfProb() 
   '''

   name = 'NO_NAME_DEFINED_YET'
   regular_exp = '???###!!!'
   INSTANCE_END_EXPRESSION = '=ready='
   datakey = 'NO KEY'
   #print data#
   columnwidth = -1
   columnheaderstr = 'NOCOLHEAD'
   problemname = ''
   problemnamelist = []

   noregexp_err = 'no reg exp defined yet'
   
   
   # the reader might behave differently depending on the solver type, due to the different output
   SOLVERTYPE_SCIP="SCIP"
   SOLVERTYPE_GUROBI="GUROBI"
   SOLVERTYPE_CPLEX="CPLEX"
   SOLVERTYPE_CBC="CBC"
   SOLVERTYPE_XPRESS="XPRESS"
   solvertype = SOLVERTYPE_SCIP

   def setProblemName(problemname):
      if problemname in StatisticReader.problemnamelist:
         counter = StatisticReader.problemnamelist.count(problemname)
         StatisticReader.problemname = problemname+'_'+repr(counter)
      else:
         StatisticReader.problemname = problemname
      StatisticReader.problemnamelist.append(problemname)
       
   setProblemName = staticmethod(setProblemName)

   def initializeForTestrun(self):
      StatisticReader.problemnamelist = []

   def getSplitLineWithRegexp(self, regular_exp, line, index = -1, startofline=False):
      length=len(regular_exp)
      tmpline = line
      if startofline:
         tmpline = line[0:length]
      list=[]
      if regular_exp in tmpline:
         list = line.split()
      else:
         return None

      if index == -1:
         return list
      else:
         assert 0 <= index
         assert index < len(list)
         return list[index]

      return None
       
   def getName(self):
      return self.name

   def endOfInstanceReached(self, line):
      if StatisticReader.INSTANCE_END_EXPRESSION in line:
         return True
      else:
         return False
   
   def extractStatistic(self, line):
      '''
      overwrite this method for own reader subclasses - make sure that a return statement different from 'None'
      is only accepted once on each instance of the log file
      '''
      return None
       
   def checkSolverType(self, line):
      changed = False
      if "Gurobi Optimizer version" in line and self.solvertype != StatisticReader.SOLVERTYPE_GUROBI:
         StatisticReader.solvertype = StatisticReader.SOLVERTYPE_GUROBI
         changed = True
      elif "SCIP version " in line and self.solvertype != StatisticReader.SOLVERTYPE_SCIP:
         StatisticReader.solvertype = StatisticReader.SOLVERTYPE_SCIP
         changed = True
      elif "Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer" in line and self.solvertype != StatisticReader.SOLVERTYPE_CPLEX:
         StatisticReader.solvertype = StatisticReader.SOLVERTYPE_CPLEX
         changed = True
      elif "Welcome to the CBC MILP Solver" in line and self.solvertype != StatisticReader.SOLVERTYPE_CBC:
         StatisticReader.solvertype = StatisticReader.SOLVERTYPE_CBC
         changed = True
      elif "FICO Xpress Optimizer" in line and self.solvertype != StatisticReader.SOLVERTYPE_XPRESS:
         StatisticReader.solvertype = StatisticReader.SOLVERTYPE_XPRESS
         changed = True
      if changed:
         print "changed solver type to", StatisticReader.solvertype
          

   def setTestRun(self, testrun):
      pass

   def execEndOfProb(self):
      '''
      overwrite this method to implement final behaviour at the end of each problem, such as setting flags
      '''
      return None

   def operateOnLine(self, line):
      
      if self.endOfInstanceReached(line):
         stat = self.execEndOfProb()
         if stat != None:
            assert StatisticReader.problemname != ''
            assert self.datakey != 'NO KEY'
#            print ' returning Statistic ', (StatisticReader.problemname, self.datakey, stat)
            return (StatisticReader.problemname, self.datakey, stat)
         else:
            return None
      else:
         self.checkSolverType(line)
         stat = self.extractStatistic(line)
         if stat != None:
            assert StatisticReader.problemname != ''
            assert self.datakey != 'NO KEY'
#            print ' returning Statistic ', (StatisticReader.problemname, self.datakey, stat)
            return (StatisticReader.problemname, self.datakey, stat)
         else:
            return None
