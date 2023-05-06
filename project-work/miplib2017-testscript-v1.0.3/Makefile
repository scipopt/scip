#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

VERSION    = 1.0.3
TIME       = 3600
TEST       = benchmark
SOLVER     = scip
MEM        = 8192
THREADS    = 0
PERMUTE    = 0
QUEUE      = M620
EXCLUSIVE  = false
CONTINUE   = false
FORCE      = false
JOBSIZE    = 1

DOXY       = doxygen

CHECKERDIR = checker

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

.PHONY: help
help:
		@echo "See README for details about the MIPLIB test environment"
		@echo
		@echo "VERSION:      $(VERSION)"
		@echo
		@echo "TARGETS:"
		@echo "** checker     -> compiles the solution checker"
		@echo "** clean       -> cleans the solution checker"
		@echo "** cmpres      -> generates solver comparison file"
		@echo "** doc         -> generates doxygen documentation"
		@echo "** eval        -> evaluate local test run"
		@echo "** evalcluster -> evaluate (slurm) cluster test run"
		@echo "** test        -> start automatic test runs local"
		@echo "** testcluster -> start automatic test runs on the (slurm) cluster"
		@echo
		@echo "PARAMETERS:"
		@echo "** MEM       -> maximum memory to use MB [8192]"
		@echo "** SOLVER    -> solver [scip]"
		@echo "** THREADS   -> number of threads (0: automatic) [0]"
		@echo "** TIME      -> time limit per instance in seconds [3600]"
		@echo "** TEST      -> test set [benchmark]"
		@echo "** PERMUTE   -> permutation to run (0 for original instance) [0]"
		@echo "** QUEUE     -> cluster queue to run on [M620]"
		@echo "** EXCLUSIVE -> exclusive cluster run [false]"
		@echo "** CONTINUE  -> continue interrupted cluster run [false]"
		@echo "** FORCE     -> run a instance on the cluster, even if files would be overwritten [false]"
		@echo "** JOBSIZE   -> number of instances per job on the cluster [1]"

.PHONY: checker
checker:
		@$(MAKE) -C $(CHECKERDIR) $^

.PHONY:		clean
clean:
		@$(MAKE) -C $(CHECKERDIR) clean

.PHONY: doc
doc:
		cd doc; $(DOXY) miplib.dxy;

.PHONY: test
test:
		@echo "run test: VERSION=$(VERSION) SOLVER=$(SOLVER) TEST=$(TEST) TIME=$(TIME) MEM=$(MEM) THREADS=$(THREADS) PERMUTE=$(PERMUTE)"
		@./scripts/run.sh $(VERSION) $(SOLVER) $(TEST) $(TIME) $(MEM) $(THREADS) $(PERMUTE);

.PHONY: testcluster
testcluster:
		@echo "run testcluster VERSION=$(VERSION) SOLVER=$(SOLVER) TEST=$(TEST) TIME=$(TIME) MEM=$(MEM) THREADS=$(THREADS) PERMUTE=$(PERMUTE) QUEUE=$(QUEUE) EXCLUSIVE=$(EXCLUSIVE) CONTINUE=$(CONTINUE) FORCE=$(FORCE) JOBSIZE=$(JOBSIZE)"
		@./scripts/runcluster.sh $(VERSION) $(SOLVER) $(TEST) $(TIME) $(MEM) $(THREADS) $(PERMUTE) $(QUEUE) $(EXCLUSIVE) $(CONTINUE) $(FORCE) $(JOBSIZE);

.PHONY: eval
eval:
		@echo "evaluate test: VERSION=$(VERSION) SOLVER=$(SOLVER) TEST=$(TEST)"
		@./scripts/evalrun.sh results/$(TEST).$(SOLVER).$(THREADS)threads.$(TIME)s.out;

.PHONY: evalcluster
evalcluster:
		@echo "evaluate testcluster: VERSION=$(VERSION) SOLVER=$(SOLVER) TEST=$(TEST) QUEUE=$(QUEUE)"
		@./scripts/evalruncluster.sh results/$(TEST).$(SOLVER).$(THREADS)threads.$(TIME)s.$(QUEUE).eval;

.PHONY: cmpres
cmpres:
		@echo "compare result tables: VERSION=$(VERSION) SOLVER=$(SOLVER) TEST=$(TEST)"
		@./scripts/allcmpres.sh results/$(TEST).$(SOLVER).res;


# --- EOF ---------------------------------------------------------------------
# DO NOT DELETE
