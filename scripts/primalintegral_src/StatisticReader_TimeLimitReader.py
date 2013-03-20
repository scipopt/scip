from StatisticReader import StatisticReader
import re

class TimeLimitReader(StatisticReader):
   '''
   extracts the time limit for an instance
   '''
   name = 'TimeLimitReader'
   solvingtimereadkeys = {StatisticReader.SOLVERTYPE_SCIP : "Solving Time (sec) :",
                  StatisticReader.SOLVERTYPE_CPLEX : "Solution time =",
                  StatisticReader.SOLVERTYPE_GUROBI : "Explored ",
                  StatisticReader.SOLVERTYPE_CBC : "Coin:Total time (CPU seconds):",
                  StatisticReader.SOLVERTYPE_XPRESS : " *** Search "}
                   
   timelimitreadkeys= {
                  StatisticReader.SOLVERTYPE_SCIP : '(SCIP> limits/time =|SCIP> set limits time)',
                  StatisticReader.SOLVERTYPE_CPLEX : 'CPLEX> New value for time limit in seconds',
                  StatisticReader.SOLVERTYPE_GUROBI : "Changed value of parameter TimeLimit to",
                  StatisticReader.SOLVERTYPE_CBC : "Coin:seconds has value",
                  StatisticReader.SOLVERTYPE_XPRESS : " @05"}
                   
   solvingtimelineindex = {
                  StatisticReader.SOLVERTYPE_SCIP : -1,
                  StatisticReader.SOLVERTYPE_CPLEX : 3,
                  StatisticReader.SOLVERTYPE_GUROBI : -2,
                  StatisticReader.SOLVERTYPE_CBC : 4,
                  StatisticReader.SOLVERTYPE_XPRESS : 5}
   columnwidth = 20
   datakey = 'TimeLimit'
   columnheaderstr = 'Timelimit'.rjust(columnwidth)
   #timelimit = 7200.00
   timelimit_reached = 'true'
   timelimit_not_reached = 'false'

   
   def __init__(self):
      self.timelimit=7200

   def extractStatistic(self, line):
        
      if re.search(self.timelimitreadkeys[StatisticReader.solvertype], line):
         self.timelimit = float(line.split()[-1])
      elif self.solvingtimereadkeys[StatisticReader.solvertype] in line:
         self.solvingtime = float(line.split()[self.solvingtimelineindex[StatisticReader.solvertype]])
      return None

   def execEndOfProb(self):
      return self.timelimit
   
   
   
