#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#@file    mps2zpl.awk
#@brief   converts mps file to zpl file
#@author  Thorsten Koch
#@author  Kati Wolter
#

function modify_str(x)
{
   if ( index(x,".") > 0 && index(x,".") == length(x))
      x = x "0";
   return x;
}
BEGIN {
   imode = 0;
   objsen = 1;
   objfound = 0;
}
{ gsub(/\r/,""); } # to remove ^M
/^\*/        { next; }
/^NAME/      { mode = 1; name = $2; }
/^OBJSENSE/  { mode = 2; }
/^ROWS/      {
   print "param objsen := " objsen ";";
   mode = 3;
}
/^COLUMNS/   {
   mode = 4; 
   first = 0;
   print "set NZO := {";  
}
/^RHS/       { mode = 5; }
/^RANGES/    { mode = 6; }
/^BOUNDS/    { mode = 7; }
/MARKER/ && /INTORG/ { imode = 1; next; }
/MARKER/ && /INTEND/ { imode = 0; next; }
NF == 1 {
   if (mode == 2) {
      if ($1 == "MAX") {
         objsen = -1;
      }
   } 
}
NF == 2 {
   if (mode == 3) { # Rows
      if ($1 == "N" && objfound)
         next;
      row_sen[$2] = $1; # N L E G
      row_rhs[$2] = 0.0;
      if ($1 == "N")
         objfound = 1;
   }
}
NF >= 3 {
   if (mode == 4 && NF == 3) {
      if (first == 1) {
         print ",";
      }
      first = 1;
      
      print " <\"" $2 "\",\"" $1 "\",", modify_str($3) ">";

      var_mode[$1]  = imode;
      var_lower[$1] = 0.0;
      var_upper[$1] = (imode == 1) ? 1.0 : 1e100;
      
   } else  if (mode == 4 && NF == 5) {
      if (first == 1) {
         print ",";
      }
      first = 1;
      print " <\"" $2 "\",\"" $1 "\",", modify_str($3) ">, <\"" $4 "\",\"" $1 "\",", modify_str($5) ">";
      var_mode[$1] = imode;
      var_lower[$1] = 0.0;
      var_upper[$1] = (imode == 1) ? 1.0 : 1e100;
   }
   # --- RHS ---------------------------------------------------
   if (mode == 5) { 
      if (NF == 2) {
         row_rhs[$1] = $2;
      } else if (NF == 3) {
         row_rhs[$2] = $3;
      } else  if (NF == 4) {
         row_rhs[$1] = $2;
         row_rhs[$3] = $4;
      } else  if (NF == 5) {
         row_rhs[$2] = $3;
         row_rhs[$4] = $5;
      }
   }
   # --- RANGES ---------------------------------------------------
   if (mode == 6) { 
      if (NF == 3) {
         row_range[$2] = $3;
      } else if (NF == 5) {
         row_range[$2] = $3;
         row_range[$4] = $5;
      }
   }

   # --- Bounds -----------------------------------------------------
   if (mode == 7) {
      if ($1 == "FX") {
         if (NF == 3) {
            var_lower[$2] = $3;
            var_upper[$2] = $3;
         } else if (NF == 4) {
            var_lower[$3] = $4;
            var_upper[$3] = $4;
         }
      } else if ($1 == "FR") {
         if (NF == 2) {
            var_lower[$2] = -1e100;
            var_upper[$2] =  1e100;
         } else if (NF == 3) {
            var_lower[$3] = -1e100;
            var_upper[$3] =  1e100;
         }
      } else if ($1 == "BV") {
         if (NF == 2) {
            var_mode [$2] = 1;
            var_lower[$2] = 0;
            var_upper[$2] = 1;
         } else if (NF == 3 && $3 != 1.0) {
            var_mode [$3] = 1;
            var_lower[$3] = 0;
            var_upper[$3] = 1;
         } else if (NF == 3 && $3 == 1.0) {
            var_mode [$2] = 1;
            var_lower[$2] = 0;
            var_upper[$2] = 1;
         } else if (NF == 4) {
            var_mode [$3] = 1;
            var_lower[$3] = 0;
            var_upper[$3] = 1;
         }
      } else if ($1 == "UP") {
         if (NF == 3) {
            var_upper[$2] = $3;
         } else if (NF == 4) {
            var_upper[$3] = $4;
         }
      } else if ($1 == "PL") {
         if (NF == 2) {
            var_upper[$2] = 1e100;
         } else if (NF == 3) {
            var_upper[$3] = 1e100;
         }
      } else if ($1 == "LO") {
         if (NF == 3) {
            var_lower[$2] = $3;
         } else if (NF == 4) {
            var_lower[$3] = $4;
         }
      } else if ($1 == "MI") {
         if (NF == 2) {
	    var_lower[$2] = -1e100;
	 } else if (NF == 3) {
            var_lower[$3] = -1e100;
         }       
      } else if ($1 == "UI") {
         if (NF == 3) {
            var_mode [$2] = 1;
            var_upper[$2] = $3;
         } else if (NF == 4) {
            var_mode [$3] = 1;
            var_upper[$3] = $4;
         }
      } else if ($1 == "LI") {
         if (NF == 3) {
            var_mode [$2] = 1;
            var_lower[$2] = $3;
         } else if (NF == 4) {
            var_mode [$3] = 1;
            var_lower[$3] = $4;
         }
      }  
   }
}

END {

   print "};";
   print "set IVAR := {";
   first = 0;
   for(v in var_mode) {      
      if (var_mode[v] == 1) {
         print (first == 1 ? "," : "") "<\"" v "\">";
         first = 1;
      }
   }
   has_ivar = first;
   print "};";
   print "set CVAR := {";
   first = 0;
   for(v in var_mode) {
      if (var_mode[v] == 0) {
         print (first == 1 ? "," : "") "<\"" v "\">";
         first = 1;
      }
   }
   has_cvar = first;
   print "};";

   print "set V := IVAR + CVAR;";
   
   print "param lower[V] := ";
   first = 0;
   for(v in var_lower) {
      if( v in var_mode )
      {
         print (first == 1 ? "," : "") "<\"" v "\">" modify_str(var_lower[v]);
         first = 1;
      }
   }
   print ";";

   print "param upper[V] := ";
   first = 0;
   for(v in var_upper) {
      if( v in var_mode )
      {     
         print (first == 1 ? "," : "") "<\"" v "\">" modify_str(var_upper[v]);
         first = 1;
      }
   }
   print ";";

   print "set RANGEROW := {";
   print "<\" dummy \">";
   for(r in row_range) {
      print ", <\"" r "\">";
      first = 1;
   }
   print "};";

   print "set ROW := { ";
   first = 0;
   for(r in row_sen) {
      print (first == 1 ? "," : "") "<\"" r "\",\"" row_sen[r] "\">";
      first = 1;
   }
   print "};";

   print "set R := proj(ROW, <1>);";
   print "param rhs[R] := ";

   first = 0;
   for(r in row_sen) {
      print (first == 1 ? "," : "") "<\"" r "\">" modify_str(row_rhs[r]);
      first = 1;
   }
   print ";";

   print "param range[RANGEROW] := ";
   print "<\" dummy \">0";
   for(r in row_range) {
      print ", <\"" r "\">" modify_str(row_range[r]);
   }
   print ";";

   if (has_ivar == 1) {
      print "var x[<v> in IVAR] integer >= if lower[v] > -1e20 then lower[v] else -infinity end <= if upper[v] < 1e20 then upper[v] else infinity end;";
   }
   if (has_cvar == 1) {
      print "var y[<v> in CVAR] real >= if lower[v] > -1e20 then lower[v] else -infinity end <= if upper[v] < 1e20 then upper[v] else infinity end;";
   }
   print "minimize obj: sum <r,\"N\"> in ROW : (";
   if (has_ivar == 1) {
      print "sum <r,c,a> in NZO with <c> in IVAR: objsen * a * x[c] ";
   }
   if (has_cvar == 1) {
      print "+ sum <r,c,a> in NZO with <c> in CVAR: objsen * a * y[c]";
   }
   print ");";

   # ==
   print "subto ce: forall <r,\"E\"> in ROW do";
   # ranged row, range sign +
   print "if <r> in RANGEROW and sgn(range[r]) > 0 then ";
   print "rhs[r] <= ";
   if (has_ivar == 1) {
      print "sum <r,c,a> in NZO with <c> in IVAR : a * x[c]";
   }
   if (has_cvar == 1) {
      print "+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c]";
   }
   print "<= rhs[r] + abs(range[r])";
   # ranged row, range sign -
   print "else if <r> in RANGEROW and sgn(range[r]) < 0 then ";
   print "rhs[r] - abs(range[r]) <= ";
   if (has_ivar == 1) {
      print "sum <r,c,a> in NZO with <c> in IVAR : a * x[c]";
   }
   if (has_cvar == 1) {
      print "+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c]";
   }
   print "<= rhs[r]";
   # normal row
   print " else ";
   if (has_ivar == 1) {
      print "sum <r,c,a> in NZO with <c> in IVAR : a * x[c]";
   }
   if (has_cvar == 1) {
      print "+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c]";
   }
   print "== rhs[r]";
   print " end end;";

   # <=
   print "subto cl: forall <r,\"L\"> in ROW do ";
   # ranged row
   print "if <r> in RANGEROW then ";
   print "rhs[r] - abs(range[r]) <= ";
   if (has_ivar == 1) {
      print "sum <r,c,a> in NZO with <c> in IVAR : a * x[c] ";
   }
   if (has_cvar == 1) {
      print "+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c] ";
   }
   print "<= rhs[r]";
   # normal row
   print " else ";
   if (has_ivar == 1) {
      print "sum <r,c,a> in NZO with <c> in IVAR : a * x[c] ";
   }
   if (has_cvar == 1) {
      print "+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c] ";
   }
   print "<= rhs[r]";
   print " end;";

   # >= 
   print "subto cg: forall <r,\"G\"> in ROW do ";
   # ranged row
   print "if <r> in RANGEROW then ";
   print "rhs[r] + abs(range[r]) >= ";
   if (has_ivar == 1) {
      print "sum <r,c,a> in NZO with <c> in IVAR : a * x[c] ";
   }
   if (has_cvar == 1) {
      print "+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c] ";
   }
   print ">= rhs[r]";
   # normal row
   print " else ";
   if (has_ivar == 1) {
      print "sum <r,c,a> in NZO with <c> in IVAR : a * x[c] ";
   }
   if (has_cvar == 1) {
      print "+ sum <r,c,a> in NZO with <c> in CVAR : a * y[c] ";
   }
   print ">= rhs[r]";
   print " end;";
}
