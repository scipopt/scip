from DataCollector import DataCollector
from StatisticReader import StatisticReader
from DataCollector_SoluFileDataCollector import SoluFileDataCollector
from StatisticReader_TimeLimitReader import TimeLimitReader
from StatisticReader_PrimalBoundReader import PrimalBoundReader
from TestRun import TestRun
import Misc

class Comparator:
   '''
   manages the collection of all log (.out) and .solu file data
   '''
   files = []
   testruns = []
   datacollector = None
   solufiledatacollector = None
   readers = []
   probnamelist = []
   datakeys = []

   datakey_err = 'ERROR: Datakey is not a suitable key'
   probname_err = 'ERROR: Probname is not in list of problem names'

   datakey_gap = 'SoluFileGap'

   EXCLUDE_REASON_TIMEOUT = 'timout'
   EXCLUDE_REASON_INFEASIBLE = 'infeasible'
   EXCLUDE_REASON_ZEROSOLUTION = 'zerosolution'
   EXCLUDE_REASON_NOSOLUTIONKNOWN = 'nosolution'
   EXCLUDE_REASON_NOOPTIMALSOLUTIONKNOWN = 'nooptsolution'
   INFINITY = 1e09
   COMP_SIGN_LE = 1
   COMP_SIGN_GE = -1

   def __listelementsdiffer(self, listname):
      for listitem in listname:
         if listname.count(listitem) > 1:
            return False
      else:
         return True

   def __init__(self, files, listofreaders = [], testsetname = '', solufilename = ''):
      self.files = files[0:len(files)]
      assert self.__listelementsdiffer(self.files)
      self.testruns = []
      self.readers = []
      self.datakeys = []
      for filename in files:
         testrun = TestRun(filename, solufilename, testsetname)
         self.testruns.append(testrun)
         testrun.settings = filename.split('.')[-2]
            
      self.datacollector = DataCollector()
      self.datacollector.registerListOfReaders(listofreaders)
      if solufilename != '':
#        print 'solufiledatacollector initialized for solufilename:', solufilename
         self.solufiledatacollector = SoluFileDataCollector()
        
      self.readers = listofreaders
      for reader in self.readers:
         self.addDataKey(reader.datakey)

   def addDataKey(self, datakey):
      self.datakeys.append(datakey)
    
   def setTestRun(self, testrun):
      self.datacollector.setTestRun(testrun)
      if self.solufiledatacollector != None:
         self.solufiledatacollector.setTestRun(testrun)

   def __makeProbNameList__(self):
      self.probnamelist = []
      for testrun in self.testruns:
         if self.probnamelist == []:
            self.probnamelist = testrun.problist
         else:
            if testrun.datacollected and testrun.problist == []:
               print testrun.getIdentification()
            assert not testrun.datacollected or testrun.problist != []
            for probname in testrun.problist:
               if not probname in self.probnamelist:
                  self.probnamelist.append(probname)
                    
                        
        

   def collectData(self):
      for testrun in self.testruns:
         self.setTestRun(testrun)
         self.datacollector.collectData()
         if self.solufiledatacollector != None:
            #print 'Collecting Solu File Data'
            assert self.solufiledatacollector.testrun == self.datacollector.testrun
            self.solufiledatacollector.collectData()
                
#      for testrun in self.testruns:
         #print testrun.problist
      self.__makeProbNameList__()
#        self.calculateData()

   def problemCompareData(self, probname, testrun1, testrun2, datakey, compsign):
      if not datakey in self.datakeys:
         raise self.datakey_err
            
#        if not testrun1.datacollected or not testrun2.datacollected:
#            print 'Collect all Data First'
        
      if not probname in self.probnamelist:
         raise self.probname_err

      data1 = float(testrun1.problemGetData(probname, datakey))
      data2 = float(testrun2.problemGetData(probname, datakey))

      return (data1 - data2) * compsign

   def getBestTestRunForProblem(self, probname, datakey, compsign):
      bestrun = None
      for testrun in self.testruns:
         if bestrun == None:
            bestrun = testrun
         elif self.problemCompareData(probname, bestrun, testrun, datakey, compsign) > 0:
            bestrun = testrun
      return bestrun

   def getBestRuns(self, datakey, compsign):
      bestruns = []
      for probname in self.probnamelist:
         bestruns.append(self.getBestTestRunForProblem(probname, datakey, compsign))
      return bestruns

   def excludeProb(self, probname, excludereasons=[]):
      for testrun in self.testruns:
         for reason in excludereasons:
            if Comparator.EXCLUDE_REASON_TIMEOUT == reason:
               if testrun.problemGetData(probname, TimeLimitReader.datakey) == TimeLimitReader.timelimit_reached:
                  return True
            elif Comparator.EXCLUDE_REASON_NOOPTIMALSOLUTIONKNOWN == reason:
               if testrun.problemGetSoluFileStatus(probname) != 'opt':
                  return True
            elif Comparator.EXCLUDE_REASON_NOSOLUTIONKNOWN == reason:
               if testrun.problemGetSoluFileStatus(probname) == 'unkn':
                  return True
            elif Comparator.EXCLUDE_REASON_ZEROSOLUTION == reason:
               if testrun.problemGetOptimalSolution(probname) == 0:
                  return True
            elif Comparator.EXCLUDE_REASON_INFEASIBLE == reason:
               if testrun.problemGetSoluFileStatus(probname) == 'inf':
                  return True
      else:
         return False

   def testrunGetKeyGeomMean(self, testrun, datakey, exclude=False):
      listofnumbers = []
      for probname in self.probnamelist:
         if not self.excludeProb(probname) or not exclude:
            listofnumbers.append(float(testrun.problemGetData(probname, datakey)))
      return Misc.listGetGeomMean(listofnumbers)

   def keyGetBestValue(self, probname, datakey, compsign):
      bestrun = self.getBestTestRunForProblem(probname, datakey, compsign)
      assert bestrun != None
      return bestrun.problemGetData(probname, datakey)

   def keyGetBestValues(self, datakey, compsign):
      values = []
      for probname in self.probnamelist:
         values.append(self.keyGetBestValue(probname, datakey, compsign))
      return values


   def testrunGetProbGapToOpt(self, testrun, probname):
      assert probname in testrun.problist
      optsol = testrun.problemGetOptimalSolution(probname)
      status = testrun.problemGetSoluFileStatus(probname)
      pb = testrun.problemGetData(probname, PrimalBoundReader.datakey)
      if status == 'opt' or status == 'best':
         return Misc.getGap(float(pb), float(optsol))
      else:
         return Misc.FLOAT_INFINITY

   def calculateData(self):
      self.calculateGaps()

   def calculateGaps(self):
      for probname in self.probnamelist:
         for testrun in self.testruns:
            gap = self.testrunGetProbGapToOpt(testrun, probname)
            data = (probname, Comparator.datakey_gap, '%8.2g'%(gap))
            testrun.addData(data)
      self.addDataKey(Comparator.datakey_gap)
