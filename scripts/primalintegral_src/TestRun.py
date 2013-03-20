import Misc
from StatisticReader_SolvingTimeReader import SolvingTimeReader
from StatisticReader_TimeLimitReader import TimeLimitReader

class TestRun:
   '''
   represents the collected data of a particular log (.out) file
   '''
   solufilename = ''
   testsetname = ''
   filename = ''
   data = {}
   solufiledata = {}
   problist = []
   nprobs = 0
   settings = 'default'
   version = ''
   lpsolver = ''
   solver = ''
   datacollected = False
   INFINITY = 1e20
   
   solufile_err = 'ERROR: No Solufile known'
   
   def __init__(self, filename = '', solufilename = '', testsetname = '' ):
      if solufilename != '':
         self.solufilename = solufilename 
   
      if testsetname != '':
         self.testsetname = testsetname
   
      if filename != '':
         self.filename = filename 
   
      self.data = {}
      self.solufiledata = {}
      self.problist= []
      self.nprobs = 0
      self.settings = ''
      self.version = ''
      self.lpsolver = ''
   
   def setDataCollected(self, datacollected):
      self.datacollected = datacollected
   
   def addProblem(self, probname):
      assert not (probname in self.problist)
      assert not (probname in self.data.keys())
      #print probname
      self.problist.append(probname)
      self.data[probname] = {}
      self.nprobs = self.nprobs + 1
   
   def __addData(self, probname, datakey, datum):
      assert probname in self.data.keys()
      if datakey in self.data[probname].keys():
         print self.data[probname].keys(), datakey
      
      assert not datakey in self.data[probname].keys()
   
      self.data[probname][datakey]= datum
   
   def addData(self, data):
      assert data != None
      probname, datakey, datum = data
      self.addInternalData(probname, datakey, datum)
      
   def addInternalData(self, probname, datakey, datum):
      if not (probname in self.problist):
         self.addProblem(probname)
      
      self.__addData(probname, datakey, datum)
      
      
   def problemGetName(self, index):
      return self.problist[index]
   
   def problemGetData(self, probname, datakey):
      return self.data[probname][datakey]
       
   def problemlistGetData(self, problemlist, datakey, numeric=True):
      array=[]
      for probname in problemlist:
         data = self.problemGetData(probname, datakey)
         if numeric:
            data = float(data)
         
         array.append(data)
      return array
   
   def getSettings(self):
      return self.settings
   
   def getVersion(self):
      return self.version
   
   def getLpSolver(self):
      return self.lpsolver
   
   def getIdentification(self):
      return self.solver + '('+self.getVersion() + ')' + self.getLpSolver() + ':' + self.getSettings()
   def getShortIdentification(self, char='_', maxlength=-1):
      return Misc.cutString(self.getSettings(), char, maxlength)
   
   def isProblemSolved(self, probname):
      return True
   
   def isProblemFeasible(self, probname):
      return True
   
   def problemGetIndex(self, probname):
      return self.problist.index(probname)
   
   def timeLimitHit(self, probname):
      return self.data[probname][SolvingTimeReader.datakey] - self.data[probname][TimeLimitReader.datakey] >= 0
   
   def getTestset(self):
      self.testsetname
       
   
   
   def problemGetOptimalSolution(self, solufileprobname):
      if self.solufiledata == {}:
         raise self.solufile_err
         
      objval = self.solufiledata[solufileprobname][0]
      
      return objval
       
   def problemGetSoluFileStatus(self, solufileprobname):
      if self.solufiledata == {}:
         raise self.solufile_err
      
      status = self.solufiledata[solufileprobname][1]
      
      return status
      
