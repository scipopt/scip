#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      *
#*                                                                           *
#*  Licensed under the Apache License, Version 2.0 (the "License");          *
#*  you may not use this file except in compliance with the License.         *
#*  You may obtain a copy of the License at                                  *
#*                                                                           *
#*      http://www.apache.org/licenses/LICENSE-2.0                           *
#*                                                                           *
#*  Unless required by applicable law or agreed to in writing, software      *
#*  distributed under the License is distributed on an "AS IS" BASIS,        *
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
#*  See the License for the specific language governing permissions and      *
#*  limitations under the License.                                           *
#*                                                                           *
#*  You should have received a copy of the Apache-2.0 license                *
#*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# replace placeholders and probe perf events in EXECNAME
# expects: EXECNAME, ERRFILE, PERFFILE
# sets: EXECNAME (with placeholders resolved)

EXECNAME="${EXECNAME/ERRFILE_PLACEHOLDER/${ERRFILE}}"
EXECNAME="${EXECNAME/RRTRACEFOLDER_PLACEHOLDER/${ERRFILE}}"
EXECNAME="${EXECNAME/PERFFILE_PLACEHOLDER/${PERFFILE}}"

# determine available perf events on this node
if echo "${EXECNAME}" | grep -q 'PERFEVENTS_PLACEHOLDER'
then
    PERF_EVENTS=""
    for event in cpu-cycles:u cpu-cycles:k instructions:u instructions:k cache-references:u cache-misses:u branch-instructions:u branch-misses:u L1-dcache-loads:u L1-dcache-load-misses:u L1-icache-load-misses:u dTLB-loads:u dTLB-load-misses:u iTLB-load-misses:u minor-faults:u minor-faults:k major-faults:u major-faults:k context-switches cpu-migrations cycle_activity.cycles_mem_any cycle_activity.stalls_mem_any cycle_activity.stalls_total task-clock; do
        if perf stat -e "${event}" true 2>/dev/null; then
            PERF_EVENTS="${PERF_EVENTS:+${PERF_EVENTS},}${event}"
        fi
    done
    EXECNAME="${EXECNAME/PERFEVENTS_PLACEHOLDER/${PERF_EVENTS}}"
fi
