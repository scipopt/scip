from math import *
from TestRun import TestRun

class DataCollector:
   '''
      acquires test run data
   '''  
   listofreaders = []
   register_err = 'register error: Reader already in list of readers'
   readernotfound_err = 'Error: Reader is not registered'
   ndata_err = 'Number Error: Readers collected data for a different number of instances'
   filestring = ''
   index = -1
   datacollected = False
   testrun = None

   def __init__(self, testrun = None, index=0):
      self.index = index
      self.listofreaders = []
      self.testrun = testrun
      if self.testrun != None:
         self.filestring = self.testrun.filename

   def setTestRun(self, testrun):
      if self.testrun != testrun:
         self.testrun = testrun
         self.filestring = self.testrun.filename
         for reader in self.listofreaders:
            reader.setTestRun(testrun)
         self.datacollected = False

   def getNReaders(self):
      return len(self.listofreaders)
   
   
   def registerReader(self, reader):
      for handl in self.listofreaders:
         if reader.getName() == handl.getName():
            raise self.register_err
            
      self.listofreaders.append(reader)
      reader.setTestRun(self.testrun)

   def registerListOfReaders(self, readers):
      for reader in readers:
         self.registerReader(reader)


   def collectData(self):
      assert(self.testrun != None)
      f = None
      try:
         f = open(self.filestring, 'r')
      except IOError:
         print 'File', self.filestring, "doesn't exist!!!" 
      
      for reader in self.listofreaders:
         reader.initializeForTestrun()
      
      for line in f:
         for reader in self.listofreaders:
            data = reader.operateOnLine(line)
            if data != None:
               self.testrun.addData(data)
            
      f.close()
      #print("Collection of data finished")
      self.datacollected = True
      self.testrun.setDataCollected(True)
      return 1
      
   def printData(self):
      ninstances = -1;
      for reader in self.listofreaders:
         if ninstances == -1:
            ninstances = reader.getNInstances()
         else:
            if ninstances != reader.getNInstances():
               raise self.ndata_err
            
         
      for reader in self.listofreaders:
         print reader.columnheaderstr,
      else:
         print
      
      for reader in self.listofreaders:
         print reader.columnwidth * '-',
      else:
         print
      
      
      for i in range(0, ninstances):
         for reader in self.listofreaders:
            print reader.getStatisticForInstance(i),
         else:
            print
         
      return 1
      
      
