import os

from StatisticReader import StatisticReader

class ProblemNameReader(StatisticReader):
   '''
      reads a problem name
   '''
   name = 'ProblemNameReader'
   regular_exp = '@01'
   datakey = 'ProblemName'
   columnwidth = 40
   columnheaderstr = 'Problem Name'.ljust(columnwidth)
   extensions = ["mps", "cip", "fzn", "pip"]

   def extractStatistic(self, line):
      fullname = self.getSplitLineWithRegexp(self.regular_exp, line, index = 1, startofline = True)
      if fullname != None:
         namewithextension=os.path.basename(fullname);
         namewithextension = os.path.splitext(namewithextension)[0]
            
         if ".mps" in namewithextension[-4:]:
            namewithextension = namewithextension[:-4]
            
         if namewithextension.endswith('.lp'):
            namewithextension=namewithextension[:-3]
         #print namewithextension
         StatisticReader.setProblemName(namewithextension)
      return None

