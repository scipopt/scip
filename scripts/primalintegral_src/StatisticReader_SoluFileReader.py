from StatisticReader import StatisticReader

class SoluFileReader(StatisticReader):
   '''
   reader for solu file data
   '''
   name = 'SoluFileReader'
   regular_exp = ['=opt=', '=inf=', '=unkn=', '=best=', '=feas=']
   actions = {}
   datakey = 'SoluFileInformation'
   statistics = {}
   conflict_err = ' Error: instance already stored in solution statistics'
   columnwidth = 12
   columnheaderstr = 'SoluFile'.rjust(columnwidth)

   INFINITY = 1e20

   def __init__(self):
      self.__initActions()

   def setTestRun(self, testrun):
      self.testrun = testrun
      if testrun != None:
         self.statistics = self.testrun.solufiledata

   def extractStatistic(self, line):
      assert self.testrun != None
      for regexp in self.regular_exp:
         if regexp in line:
            newaction = self.actions[regexp]
            newaction(line)
      else:
         return None
            

   def __storeToStatistics(self, instance, objval, status):
      if not instance in self.statistics.keys():
         self.statistics[instance] = (objval, status)
      else:
         raise self.conflict_err


   def __newOptInstance(self, line):
      splittedline = line.split()
      assert splittedline[0] == '=opt='
      instance = splittedline[1]
      objval = splittedline[2]

      self.__storeToStatistics(instance, objval, status='opt')

   def __newInfInstance(self, line):
      splittedline = line.split()
      assert splittedline[0] == '=inf='
      instance = splittedline[1]
      objval = self.INFINITY

      self.__storeToStatistics(instance,  objval, status='inf')

   def __newUnknInstance(self, line):
      splittedline = line.split()
      assert splittedline[0] == '=unkn='
      instance = splittedline[1]
      objval = self.INFINITY

      self.__storeToStatistics(instance, objval, status='unkn')


   def __newBestInstance(self, line):
      splittedline = line.split()
      assert splittedline[0] == '=best='
      instance = splittedline[1]
      objval = splittedline[2]

      self.__storeToStatistics(instance, objval, status='best')

   def __newFeasInstance(self, line):
      splittedline = line.split()
      assert splittedline[0] == '=best='
      instance = splittedline[1]
      objval = self.INFINITY
        
      self.__storeToStatistics(instance, objval, status='feas')



   def __initActions(self):
      self.actions = {}
      self.actions['=opt='] = self.__newOptInstance
      self.actions['=inf='] = self.__newInfInstance
      self.actions['=unkn='] = self.__newUnknInstance
      self.actions['=best='] = self.__newBestInstance
      self.actions['=feas='] = self.__newFeasInstance
        
