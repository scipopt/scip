/* allocated memory for profit array, solution array, non solution array, and items array */
   SCIP_CALL( SCIPallocBufferArray(scip, &profits, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, nitems) );

   /* copy the weights and ids array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &weights, pricerdata->weights, nitems) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &ids, pricerdata->ids, nitems) );
   
   nsolitems = 0;
   nnonsolitems = 0;

   /* fill profit array with dual solution value */
   for( i = 0; i < nitems; i++ ) 
   {
      profits[i] = SCIPgetDualsolSetppc(scip, conss[i]);
   }
   
#ifdef SCIP_DEBUG
   /* debug output; this one turned on using "#define SCIP_DEBUG" */

   SCIPdebugMessage("   nitems: %d capacity: %"SCIP_LONGINT_FORMAT"\n", nitems, capacity);
   SCIPdebugMessage(" %4s %4s %7s\n", "#", "dual", "weights");
   for( i = 0; i < nitems; ++i ) 
   {
      SCIPdebugMessage("%4d %4.2f %7"SCIP_LONGINT_FORMAT"\n", ids[i], profits[i], weights[i]);
   }
#endif
   
   /* solve the knapsack sack problem */
   SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, weights, profits, capacity, ids, 
         solitems, nonsolitems, &nsolitems, &nnonsolitems, &redcost, &success) );
   
#ifdef SCIP_DEBUG
   /* debug output; this one turned on using "#define SCIP_DEBUG" */

   SCIPdebugMessage("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
   SCIPdebugMessage("X Knapsack algorithm calculated reduced costs of %f X \n",redcost);
   SCIPdebugMessage("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");   

   if( success )
   {
      int item;
      SCIP_Longint usedcapacity;

      usedcapacity = 0;
      
      SCIPdebugMessage(" nitems: %d capacity: %"SCIP_LONGINT_FORMAT"\n", nitems, capacity);
      SCIPdebugMessage(" %4s %7s\n", "#", "weights");
      for( i = 0; i < nsolitems; ++i ) 
      {
         item  = solitems[i];
         usedcapacity += weights[item];
         SCIPdebugMessage(" %4d %7"SCIP_LONGINT_FORMAT"\n", item, pricerdata->weights[item]);
      }
      assert(usedcapacity <= capacity);
   }
   else
   {
      SCIPdebugMessage("the knapsack problem was not solved correctly\n");
   }
#endif
   
   if( success )
   {
      /* check if we found a solution with negative reduced cost */
      if( SCIPisFeasGT(scip, redcost, 1.0) ) 
      {
         SCIP_VAR* var;
         SCIP_VARDATA* vardata;
         char strtmp[SCIP_MAXSTRLEN];
      
         /* create a suitable variable name; this might help for debugging */
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "items"); 
      
         /* create variable name and collect constraint ids for the variable data */
         for( i = 0; i < nsolitems; ++i )
         {
            (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", solitems[i]); 
            strcat(varname, strtmp); 
         }
         
         /* create the variable data for the variable; the variable data contains the information in which constraints
          * the variable appears */
         SCIP_CALL( SCIPvardataCreateBinpacking(scip, &vardata, solitems, nsolitems) );
            
         /* create variable for a new column */
         SCIP_CALL( SCIPcreateVarBinpacking(scip, &var, varname, 1.0, FALSE, TRUE, vardata) );

         /* add the new variable to the pricer store */
         SCIP_CALL( SCIPaddPricedVar(scip, var, 1.0) );
         
         /* change the upper bound of the binary variable to lazy since the upper bound is already enforced 
          * due to the objective function the set covering constraint;
          * The reason for doing is that, is to avoid the bound of x <= 1 in the LP relaxation since this bound
          * constraint would produce a dual variable which might have a positive reduced cost 
          */
         SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );
         
         /* add the new variable to the corresponding set covering constarints */
         for( i = 0; i < nsolitems; ++i )
         {
            SCIP_CALL( SCIPaddCoefSetppc(scip, conss[solitems[i]], var) );
         }
         
         /* release the new variable such that SCIP takes care it */
         SCIP_CALL( SCIPreleaseVar(scip, &var) );      
      }

      (*result) = SCIP_SUCCESS;
   }
   else
   {
      /* since the knapsack problem was aborted, there is no guarantee that the current LP solution is optimal; to
       * indicate that to SCIP we set the result pointer to SCIP_DIDNOTRUN and continue with an incomplete pricing */
      (*result) = SCIP_DIDNOTRUN;
   }
      
   /* free buffer arrarys */
   SCIPfreeBufferArray(scip, &ids);
   SCIPfreeBufferArray(scip, &weights);
   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &profits);
