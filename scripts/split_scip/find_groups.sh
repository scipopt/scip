grep -A5 "@" ../../src/scip/scip.h | grep addtogroup | grep -oP "[^ ]+$" > groups.list
