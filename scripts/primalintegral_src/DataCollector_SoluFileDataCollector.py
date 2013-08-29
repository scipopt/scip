from DataCollector import DataCollector
from StatisticReader_SoluFileReader import SoluFileReader

class SoluFileDataCollector(DataCollector):
   '''
      a data collector to be run on .solu files.
   '''
      
   def __init__(self, testrun = None, index = 0):
      self.registerReader(SoluFileReader())
      if testrun != None:
         self.setTestrun(testrun)
         self.filestring = self.testrun.solufilename
      
   def setTestRun(self, testrun):
      if self.testrun != testrun:
         self.testrun = testrun
         self.filestring = self.testrun.solufilename
         for reader in self.listofreaders:
            reader.setTestRun(testrun)
         self.datacollected = False
