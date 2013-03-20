from StatisticReader import StatisticReader

class SolvingTimeReader(StatisticReader):
   '''
   extracts the solving time
   '''
   name = 'SolvingTimeReader'
   regular_exp = 'Solving Time (sec) :'
   columnwidth = 12
   datakey = 'SolvingTime'
   columnheaderstr = 'Time(s)'.rjust(columnwidth)
    
   solvingtimereadkeys = {
      StatisticReader.SOLVERTYPE_SCIP : "Solving Time (sec) :",
      StatisticReader.SOLVERTYPE_CPLEX : "Solution time =",
      StatisticReader.SOLVERTYPE_GUROBI : "Explored ",
      StatisticReader.SOLVERTYPE_CBC : "Coin:Total time (CPU seconds):",
      StatisticReader.SOLVERTYPE_XPRESS : " *** Search "
   }
                   
   solvingtimelineindex = {
      StatisticReader.SOLVERTYPE_SCIP : -1,
      StatisticReader.SOLVERTYPE_CPLEX : 3,
      StatisticReader.SOLVERTYPE_GUROBI : -2,
      StatisticReader.SOLVERTYPE_CBC : 4,
      StatisticReader.SOLVERTYPE_XPRESS : 5
   }
    
   returned = False
   DEFAULT_SOLVINGTIME=3600.0

   def extractStatistic(self, line):
      if self.solvingtimereadkeys[StatisticReader.solvertype] in line:
         solvingtime = line.split()[self.solvingtimelineindex[StatisticReader.solvertype]]
         self.returned = True
         return float(solvingtime)
      else:
         return None

   def execEndOfProb(self):
      if not self.returned:
         return SolvingTimeReader.DEFAULT_SOLVINGTIME
      else:
         self.returned = False
         return None