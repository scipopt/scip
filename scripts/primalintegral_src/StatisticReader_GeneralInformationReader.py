# -*- coding: utf-8 -*-
from StatisticReader import StatisticReader

class GeneralInformationReader(StatisticReader):
   '''
      this reader extracts general information from a log file, that is, the solver version, the used LP solver
      in case of SCIP and the settings
   '''
   name = 'GeneralInformationReader'
   regular_exp = ['loaded parameter file', 'default.set']
   columnwidth = 20
   actions = {}
   testrun = None
   columnheaderstr = 'Version'.rjust(columnwidth)
   
   versionkeys = {StatisticReader.SOLVERTYPE_SCIP : 'SCIP version',
                  StatisticReader.SOLVERTYPE_CPLEX : "Welcome to IBM(R) ILOG(R) CPLEX(R)",
                  StatisticReader.SOLVERTYPE_GUROBI : "Gurobi Optimizer version",
                  StatisticReader.SOLVERTYPE_CBC : "Version:",
                  StatisticReader.SOLVERTYPE_XPRESS : "FICO Xpress Optimizer"}
                  
   versionlineindices = {StatisticReader.SOLVERTYPE_SCIP : 2,
                         StatisticReader.SOLVERTYPE_CPLEX : -1,
                         StatisticReader.SOLVERTYPE_GUROBI : -1,
                         StatisticReader.SOLVERTYPE_CBC : -1,
                         StatisticReader.SOLVERTYPE_XPRESS: 4}
   
   def __init__(self, testrun = None):
      self.testrun = testrun
   
   def setTestRun(self, testrun):
      self.testrun = testrun
   
   def extractStatistic(self, line):
      if self.testrun == None:
         return None
         
      if self.testrun.version == '' and self.versionkeys[StatisticReader.solvertype] in line:
         if self.testrun.solver == '':
            self.testrun.solver = StatisticReader.solvertype
            
         if StatisticReader.solvertype == StatisticReader.SOLVERTYPE_SCIP:
            self.__handleVersion(line)
         else:
            self.testrun.version = line.split()[self.versionlineindices[StatisticReader.solvertype]]
      elif self.testrun.settings == '' and StatisticReader.solvertype == StatisticReader.SOLVERTYPE_SCIP:
         if 'default.set' in line:
            self.testrun.settings = 'default'
         elif 'loaded parameter file' in line:
            splittedsettings = self.getSplitLineWithRegexp('loaded parameter file', line, index=3, startofline = True);
            if splittedsettings != None:
               splittedsettings = splittedsettings.split('/')
               settings = splittedsettings[- 1].split('.set')[0]
               self.testrun.settings = settings
      
      elif self.testrun.settings == '' and StatisticReader.solvertype == StatisticReader.SOLVERTYPE_CPLEX:
         if "CPLEX> Non-default parameters written to file" in line:
            self.testrun.settings=line.split('.')[-3]
            print line.split('.')
         
             
      return None
      
   def __handleVersion(self, line):
      
      splittedline = self.getSplitLineWithRegexp('SCIP version', line)
      if splittedline != None:
         version = splittedline[2]
         if self.testrun.version != '':
            assert self.testrun.version == version
         else:
            self.testrun.version = version
         
         lpsolver = splittedline[splittedline.index('[LP') + 2].rstrip(']')
         if self.testrun.lpsolver != '':
            assert self.testrun.lpsolver == lpsolver
         else:
            self.testrun.lpsolver = lpsolver

















    
