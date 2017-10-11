# Analyze infeasible multi-commodity flow IP/LP via indicator constraints
# We relax the mutual capacity bounds.
#
#  Mutual capacity file (*.mut):
#
#  < mutual capacity pointer > , < mutual capacity >
#
#  Arc file (*.arc):
#
#  < arc name > , < from node > , < to node > , < commodity > , < cost > ,
#  < capacity > , < mutual capacity pointer >
#
#  Arc name is an integer between 1 and the number of arcs (differently from
#  the original mnetgen format), that is necessary to distinguish between
#  multiple instances of an arc (i, j) for the same commodity, that are
#  permitted
#
#  Node supply file (*.sup):
#
#  < node > , < commodity > , < supply >
#
#  Problem description file (*.nod, only if FOUR_F == 1)
#
#  < commodities > , < nodes > , < arcs > , < capacitated arcs >

# to get an infeasible instance multiply supplies by the following factor
param factor := 4;


# read full arcs (id, from, to, commodity)
set AF := { read "64-4-1.arc" as "<1n,2n,3n,4n>" };

# get nodes (projection onto second and third component)
set N := proj(AF, <2>) union proj(AF, <3>);

# get commodities
set C := proj(AF, <4>);

# read arcs (id, commodity)
set A := { read "64-4-1.arc" as "<1n,4n>" };


# read mutual capcities ID
set MID := { read "64-4-1.mut" as "<1n>" };

# read mutual capcities
param MC[MID] := read "64-4-1.mut" as "<1n> 2n";


# read supplies
param S[N*C] := read "64-4-1.sup" as "<1n,2n> 3n" default 0;

# read commodity capacities
param UK[A] := read "64-4-1.arc" as "<1n,4n> 6n";

# read mutual capacity ids
param MAID[A] := read "64-4-1.arc" as "<1n,4n> 7n";

# read cost
param COST[A] := read "64-4-1.arc" as "<1n,4n> 5n";


# flow variables
var x[<id,c> in A] real >= 0;

# indicator variables for each capacity
var y[MID] binary;

# objective
minimize obj: sum <mid> in MID: y[mid];

# flow conservation (unchanged)
subto fc:
forall <v> in N:
forall <c> in C:
sum <id,u,v,c> in AF: x[id,c] - sum <id,v,u,c> in AF: x[id,c] == -factor * S[v,c];

# mutual capacities (possibly relaxed)
subto mutCap:
forall <mid> in MID:
vif y[mid] == 0 then sum <id,c> in A with MAID[id,c] == mid: x[id,c] <= MC[mid] end, indicator;
