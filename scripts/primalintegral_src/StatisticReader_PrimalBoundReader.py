from StatisticReader import StatisticReader

class PrimalBoundReader(StatisticReader):
   name = 'PrimalBoundReader'
   regular_exp = 'Primal Bound       :'
   datakey = 'PrimalBound'
   columnwidth = 18
   columnheaderstr = 'Primal Bound'.rjust(columnwidth)
   maxndigits = 3

   def extractStatistic(self, line):
      if self.regular_exp in line:
         fullline=line.split()
         pb = float(fullline[3])
         return pb
      else:
         return None
