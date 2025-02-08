/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   datatree.c
 * @brief  unit tests for the methods of SCIP_DATATREE
 * @author Mohammed Ghannam
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "include/scip_test.h"
#include "scip/pub_datatree.h"
#include "scip/scip_datatree.h"
#include "scip/struct_datatree.h"

/** GLOBAL VARIABLES **/
static SCIP* scip;

/** TEST SUITES **/
static void setup(void)
{
   scip = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
}

static void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

TestSuite(datatree, .init = setup, .fini = teardown);

Test(datatree, test_create)
{
   SCIP_DATATREE* datatree = NULL;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 10) );

   cr_assert_not_null(datatree, "DataTree is NULL");
   cr_assert_eq(datatree->nitems, 0, "DataTree size is not 0");
   cr_assert_eq(datatree->itemssize, 10, "DataTree capacity is not 10");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_expands_beyond_capacity)
{
   SCIP_DATATREE *datatree = NULL;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 2) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "test_long", 10) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "test_long2", 20) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "test_long3", 30) );

   cr_assert_eq(datatree->nitems, 3, "DataTree size is not 3");
   cr_assert(datatree->itemssize >= 4, "DataTree capacity is not 4");

   cr_assert_str_eq(datatree->items[0].name, "test_long", "Name is not test_long");
   cr_assert_eq(datatree->items[0].value.type, SCIP_DATATREE_LONG, "Value type is not SCIP_DATATREE_LONG");
   cr_assert_eq(datatree->items[0].value.data.as_long, 10, "Value is not 10");

   cr_assert_str_eq(datatree->items[1].name, "test_long2", "Name is not test_long2");
   cr_assert_eq(datatree->items[1].value.type, SCIP_DATATREE_LONG, "Value type is not SCIP_DATATREE_LONG");
   cr_assert_eq(datatree->items[1].value.data.as_long, 20, "Value is not 20");

   cr_assert_str_eq(datatree->items[2].name, "test_long3", "Name is not test_long3");
   cr_assert_eq(datatree->items[2].value.type, SCIP_DATATREE_LONG, "Value type is not SCIP_DATATREE_LONG");
   cr_assert_eq(datatree->items[2].value.data.as_long, 30, "Value is not 30");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_create_without_size_hint)
{
   SCIP_DATATREE *datatree = NULL;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, -1) );

   cr_assert_not_null(datatree, "DataTree is NULL");
   cr_assert_eq(datatree->nitems, 0, "DataTree size is not 0");
   cr_assert(datatree->itemssize > 0, "DataTree capacity is not 0");

   SCIP_CALL(SCIPinsertDatatreeLong(scip, datatree, "test_long", 10));
   cr_assert_eq(datatree->nitems, 1, "DataTree size is not 1");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_insert_bool)
{
   SCIP_DATATREE *datatree = NULL;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeBool(scip, datatree, "test_bool", TRUE) );

   cr_assert_eq(datatree->nitems, 1, "DataTree size is not 1");
   cr_assert_str_eq(datatree->items[0].name, "test_bool", "Name is not test_bool");
   cr_assert_eq(datatree->items[0].value.type, SCIP_DATATREE_BOOL, "Value type is not SCIP_DATATREE_BOOL");
   cr_assert_eq(datatree->items[0].value.data.as_bool, TRUE, "Value is not TRUE");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_insert_long)
{
   SCIP_DATATREE *datatree = NULL;

   SCIP_CALL(SCIPcreateDatatree(scip, &datatree, 1));
   SCIP_CALL(SCIPinsertDatatreeLong(scip, datatree, "test_long", 10));

   cr_assert_eq(datatree->nitems, 1, "DataTree size is not 1");
   cr_assert_str_eq(datatree->items[0].name, "test_long", "Name is not test_long");
   cr_assert_eq(datatree->items[0].value.type, SCIP_DATATREE_LONG, "Value type is not SCIP_DATATREE_LONG");
   cr_assert_eq(datatree->items[0].value.data.as_long, 10, "Value is not 10");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_insert_real)
{
   SCIP_DATATREE *datatree = NULL;

   SCIP_CALL(SCIPcreateDatatree(scip, &datatree, 1));
   SCIP_CALL(SCIPinsertDatatreeReal(scip, datatree, "test_real", 10.0));

   cr_assert_eq(datatree->nitems, 1, "DataTree size is not 1");
   cr_assert_str_eq(datatree->items[0].name, "test_real", "Name is not test_real");
   cr_assert_eq(datatree->items[0].value.type, SCIP_DATATREE_REAL, "Value type is not SCIP_DATATREE_REAL");
   cr_assert_eq(datatree->items[0].value.data.as_real, 10.0, "Value is not 10.0");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_insert_string)
{
   SCIP_DATATREE *datatree = NULL;

   SCIP_CALL(SCIPcreateDatatree(scip, &datatree, 1));
   SCIP_CALL(SCIPinsertDatatreeString(scip, datatree, "test_string", "test"));

   cr_assert_eq(datatree->nitems, 1, "DataTree size is not 1");
   cr_assert_str_eq(datatree->items[0].name, "test_string", "Name is not test_string");
   cr_assert_eq(datatree->items[0].value.type, SCIP_DATATREE_STRING, "Value type is not SCIP_DATATREE_STRING");
   cr_assert_str_eq(datatree->items[0].value.data.as_string, "test", "Value is not test");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_insert_long_array)
{
   SCIP_DATATREE *datatree = NULL;
   SCIP_Longint long_values[3] = {10, 20, 30};

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeLongArray(scip, datatree, "int_array", long_values, 3) );

   cr_assert_eq(datatree->nitems, 1, "DataTree size is not 1");
   cr_assert_str_eq(datatree->items[0].name, "int_array", "Name is not int_array");
   cr_assert_eq(datatree->items[0].value.type, SCIP_DATATREE_LONGARRAY, "Value type is not SCIP_DATATREE_LONGARRAY");
   cr_assert_eq(datatree->items[0].value.nvalues, 3, "Array size is not 3");
   for (int i = 0; i < 3; ++i) {
      cr_assert_eq(datatree->items[0].value.data.as_longarray[i],
                   long_values[i], "Value is not %lld but %lld", long_values[i],
                   datatree->items[0].value.data.as_longarray[i]);
   }

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_insert_int_array)
{
   SCIP_DATATREE *datatree = NULL;
   int int_values[3] = {10, 20, 30};

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeIntArray(scip, datatree, "int_array", int_values, 3) );

   cr_assert_eq(datatree->nitems, 1, "DataTree size is not 1");
   cr_assert_str_eq(datatree->items[0].name, "int_array", "Name is not int_array");
   cr_assert_eq(datatree->items[0].value.type, SCIP_DATATREE_LONGARRAY, "Value type is not SCIP_DATATREE_LONGARRAY");
   cr_assert_eq(datatree->items[0].value.nvalues, 3, "Array size is not 3");
   for (int i = 0; i < 3; ++i) {
      cr_assert_eq(datatree->items[0].value.data.as_longarray[i],
                   int_values[i], "Value is not %d but %lld", int_values[i],
                   datatree->items[0].value.data.as_longarray[i]);
   }

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_insert_bool_array)
{
   SCIP_DATATREE *datatree = NULL;
   SCIP_Bool bool_values[3] = {TRUE, FALSE, TRUE};

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeBoolArray(scip, datatree, "bool_array", bool_values, 3) );

   cr_assert_eq(datatree->nitems, 1, "DataTree size is not 1");
   cr_assert_str_eq(datatree->items[0].name, "bool_array", "Name is not bool_array");
   cr_assert_eq(datatree->items[0].value.type, SCIP_DATATREE_BOOLARRAY, "Value type is not SCIP_DATATREE_BOOLARRAY");
   cr_assert_eq(datatree->items[0].value.nvalues, 3, "Array size is not 3");
   for (int i = 0; i < 3; ++i) {
      cr_assert_eq(datatree->items[0].value.data.as_boolarray[i],
                   bool_values[i], "Value is not %d but %d", bool_values[i],
                   datatree->items[0].value.data.as_boolarray[i]);
   }

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_insert_real_array)
{
   SCIP_DATATREE *datatree = NULL;
   SCIP_Real real_values[3] = {10.5, 20.5, 30.5};

   SCIP_CALL( SCIPcreateDatatree( scip, &datatree, 1 ));
   SCIP_CALL( SCIPinsertDatatreeRealArray( scip, datatree, "real_array", real_values, 3 ));

   cr_assert_eq( datatree->nitems, 1, "DataTree size is not 1" );
   cr_assert_str_eq( datatree->items[0].name, "real_array", "Name is not real_array" );
   cr_assert_eq( datatree->items[0].value.type, SCIP_DATATREE_REALARRAY, "Value type is not SCIP_DATATREE_REALARRAY" );
   cr_assert_eq( datatree->items[0].value.nvalues, 3, "Array size is not 3" );
   for ( int i = 0; i < 3; ++i )
   {
      cr_assert_eq( datatree->items[0].value.data.as_realarray[i],
                    real_values[i], "Value is not %f but %f", real_values[i],
                    datatree->items[0].value.data.as_realarray[i]);
   }

   SCIPfreeDatatree( scip, &datatree );
}


Test(datatree, test_insert_string_array)
{
   SCIP_DATATREE *datatree = NULL;
   const char* string_values[3] = {"test1", "test2", "test3"};

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeStringArray(scip, datatree, "string_array", string_values, 3) );

   cr_assert_eq(datatree->nitems, 1, "DataTree size is not 1");
   cr_assert_str_eq(datatree->items[0].name, "string_array", "Name is not string_array");
   cr_assert_eq(datatree->items[0].value.type, SCIP_DATATREE_STRINGARRAY, "Value type is not SCIP_DATATREE_STRINGARRAY");
   cr_assert_eq(datatree->items[0].value.nvalues, 3, "Array size is not 3");
   for (int i = 0; i < 3; ++i) {
      cr_assert_str_eq(datatree->items[0].value.data.as_stringarray[i],
                       string_values[i], "Value is not %s but %s", string_values[i],
                       datatree->items[0].value.data.as_stringarray[i]);
   }

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_insert_subtree)
{
   SCIP_DATATREE *datatree = NULL;
   SCIP_DATATREE *subdatatree = NULL;

   SCIP_CALL(SCIPcreateDatatree(scip, &datatree, 1));
   SCIP_CALL(SCIPcreateDatatree(scip, &subdatatree, 1));
   SCIP_CALL(SCIPinsertDatatreeLong(scip, subdatatree, "sub_long", 20));

   SCIP_CALL(SCIPinsertDatatreeTree(scip, datatree, "subdatatree", subdatatree));

   cr_assert_eq(datatree->nitems, 1, "DataTree size is not 1");
   cr_assert_str_eq(datatree->items[0].name, "subdatatree", "Name is not subdatatree");
   cr_assert_eq(datatree->items[0].value.type, SCIP_DATATREE_DATATREE, "Value type is not SCIP_DATATREE_STORE");
   cr_assert_eq(datatree->items[0].value.data.as_dtree, subdatatree, "Value is not subdatatree");

   cr_assert_eq(subdatatree->nitems, 1, "Subdatatree size is not 1");
   cr_assert_str_eq(subdatatree->items[0].name, "sub_long", "Subtree name is not sub_long");
   cr_assert_eq(subdatatree->items[0].value.type, SCIP_DATATREE_LONG, "Subtree value type is not SCIP_DATATREE_LONG");
   cr_assert_eq(subdatatree->items[0].value.data.as_long, 20, "Subtree value is not 20");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_long_getter)
{
   SCIP_DATATREE* datatree = NULL;
   SCIP_Longint long_value;
   SCIP_Real real_value;
   const char* string_value;
   SCIP_DATATREE* subdatatree_value;
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "test_long", 10) );

   // Correct usage
   SCIP_CALL( SCIPdatatreeGetLong(datatree, "test_long", &long_value) );
   cr_assert_eq(long_value, 10, "Value is not 10");

   // Wrong usage
   cr_redirect_stderr();
   retcode = SCIPdatatreeGetReal(datatree, "test_long", &real_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetString(datatree, "test_long", &string_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetTree(datatree, "test_long", &subdatatree_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_real_getter)
{
   SCIP_DATATREE* datatree = NULL;
   SCIP_Longint long_value;
   SCIP_Real real_value;
   const char* string_value;
   SCIP_DATATREE* subdatatree_value;
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "test_real", 10.0) );

   // Correct usage
   SCIP_CALL( SCIPdatatreeGetReal(datatree, "test_real", &real_value) );
   cr_assert_eq(real_value, 10.0, "Value is not 10.0");

   // Wrong usage
   cr_redirect_stderr();
   retcode = SCIPdatatreeGetLong(datatree, "test_real", &long_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetString(datatree, "test_real", &string_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetTree(datatree, "test_real", &subdatatree_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_string_getter)
{
   SCIP_DATATREE* datatree = NULL;
   SCIP_Longint long_value;
   SCIP_Real real_value;
   const char* string_value;
   SCIP_DATATREE* subdatatree_value;
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeString(scip, datatree, "test_string", "test") );

   // Correct usage
   SCIP_CALL( SCIPdatatreeGetString(datatree, "test_string", &string_value) );
   cr_assert_str_eq(string_value, "test", "Value is not test");

   // Wrong usage
   cr_redirect_stderr();
   retcode = SCIPdatatreeGetLong(datatree, "test_string", &long_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetReal(datatree, "test_string", &real_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetTree(datatree, "test_string", &subdatatree_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_bool_getter)
{
   SCIP_DATATREE* datatree = NULL;
   SCIP_Bool bool_value;
   SCIP_Longint long_value;
   SCIP_Real real_value;
   const char* string_value;
   SCIP_DATATREE* subdatatree_value;
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeBool(scip, datatree, "test_bool", TRUE) );

   // Correct usage
   SCIP_CALL( SCIPdatatreeGetBool(datatree, "test_bool", &bool_value) );
   cr_assert_eq(bool_value, TRUE, "Value is not TRUE");

   // Wrong usage
   cr_redirect_stderr();
   retcode = SCIPdatatreeGetLong(datatree, "test_bool", &long_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetReal(datatree, "test_bool", &real_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetString(datatree, "test_bool", &string_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetTree(datatree, "test_bool", &subdatatree_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_bool_array_getter)
{
   SCIP_DATATREE* datatree = NULL;
   SCIP_Bool bool_values[3] = {TRUE, FALSE, TRUE};
   SCIP_Bool* bool_array;
   SCIP_Longint long_value;
   SCIP_Real real_value;
   const char* string_value;
   SCIP_DATATREE* subdatatree_value;
   SCIP_RETCODE retcode;
   int nvalues;


   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeBoolArray(scip, datatree, "bool_array", bool_values, 3) );

   // Correct usage
   SCIP_CALL( SCIPdatatreeGetBoolArray(datatree, "bool_array", &bool_array, &nvalues) );
   for (int i = 0; i < 3; ++i) {
      cr_assert_eq(bool_array[i], bool_values[i], "Value is not %d", bool_values[i]);
   }

   // Wrong usage
   cr_redirect_stderr();
   retcode = SCIPdatatreeGetLong(datatree, "bool_array", &long_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetReal(datatree, "bool_array", &real_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetString(datatree, "bool_array", &string_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetTree(datatree, "bool_array", &subdatatree_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_long_array_getter)
{
   SCIP_DATATREE* datatree = NULL;
   SCIP_Longint long_values[3] = {10, 20, 30};
   SCIP_Longint* long_array;
   SCIP_Real real_value;
   const char* string_value;
   SCIP_DATATREE* subdatatree_value;
   SCIP_RETCODE retcode;
   int nvalues;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeLongArray(scip, datatree, "long_array", long_values, 3) );

   // Correct usage
   SCIP_CALL( SCIPdatatreeGetLongArray(datatree, "long_array", &long_array, &nvalues) );
   for (int i = 0; i < 3; ++i) {
      cr_assert_eq(long_array[i], long_values[i], "Value is not %lld", long_values[i]);
   }

   // Wrong usage
   cr_redirect_stderr();
   retcode = SCIPdatatreeGetReal(datatree, "long_array", &real_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetString(datatree, "long_array", &string_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetTree(datatree, "long_array", &subdatatree_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_real_array_getter)
{
   SCIP_DATATREE* datatree = NULL;
   SCIP_Real real_values[3] = {10.5, 20.5, 30.5};
   SCIP_Real* real_array;
   SCIP_Longint long_value;
   const char* string_value;
   SCIP_DATATREE* subdatatree_value;
   SCIP_RETCODE retcode;
   int nvalues;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeRealArray(scip, datatree, "real_array", real_values, 3) );

   // Correct usage
   SCIP_CALL( SCIPdatatreeGetRealArray(datatree, "real_array", &real_array, &nvalues) );
   for (int i = 0; i < 3; ++i) {
      cr_assert_float_eq(real_array[i], real_values[i], 1e-9, "Value is not %f", real_values[i]);
   }

   // Wrong usage
   cr_redirect_stderr();
   retcode = SCIPdatatreeGetLong(datatree, "real_array", &long_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetString(datatree, "real_array", &string_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetTree(datatree, "real_array", &subdatatree_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");

   SCIPfreeDatatree(scip, &datatree);
}


Test(datatree, test_string_array_getter)
{
   SCIP_DATATREE* datatree = NULL;
   const char* string_values[3] = {"test1", "test2", "test3"};
   char** string_array;
   SCIP_Longint long_value;
   SCIP_Real real_value;
   SCIP_DATATREE* subdatatree_value;
   SCIP_RETCODE retcode;
   int nvalues;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeStringArray(scip, datatree, "string_array", string_values, 3) );

   // Correct usage
   SCIP_CALL( SCIPdatatreeGetStringArray(datatree, "string_array", &string_array, &nvalues) );
   for (int i = 0; i < 3; ++i) {
      cr_assert_str_eq(string_array[i], string_values[i], "Value is not %s", string_values[i]);
   }

   // Wrong usage
   cr_redirect_stderr();
   retcode = SCIPdatatreeGetLong(datatree, "string_array", &long_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetReal(datatree, "string_array", &real_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");
   retcode = SCIPdatatreeGetTree(datatree, "string_array", &subdatatree_value);
   cr_assert_eq(retcode, SCIP_ERROR, "Function should return SCIP_ERROR");

   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_write)
{
   SCIP_DATATREE* datatree = NULL;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 3) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "test_long", 10) );
   SCIP_CALL( SCIPinsertDatatreeString(scip, datatree, "test_string", "test") );
   SCIP_Longint long_values[2] = {10, 20};
   SCIP_CALL( SCIPinsertDatatreeLongArray(scip, datatree, "test_long_array", long_values, 2) );
   SCIP_DATATREE* subdatatree = NULL;
   SCIP_CALL( SCIPcreateDatatreeInTree( scip, datatree, &subdatatree, "test_subtree", 3 ) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, subdatatree, "sub_long", 10) );
   SCIP_CALL( SCIPinsertDatatreeString(scip, subdatatree, "sub_string", "sub") );
   SCIP_CALL( SCIPinsertDatatreeLongArray(scip, subdatatree, "sub_long_array", long_values, 2) );


   // Write to file
   FILE *file = fopen("test_write.txt", "w");
   SCIP_CALL( SCIPwriteDatatreeJson(scip, file, datatree) );
   fclose(file); // Close the write file handle

   // Reopen file for reading and check content
   file = fopen("test_write.txt", "r");
   const char* expected = "{\"test_long\": 10, \"test_string\": \"test\", \"test_long_array\": [10, 20], \"test_subtree\": {\"sub_long\": 10, \"sub_string\": \"sub\", \"sub_long_array\": [10, 20]}}";
   char buffer[1024];
   fgets(buffer, 1024, file);
   cr_assert_str_eq(buffer, expected, "File content is not as expected");
   fclose(file);

   // Cleanup
   remove("test_write.txt");
   SCIPfreeDatatree(scip, &datatree);
}

Test(datatree, test_write_json_bool)
{
   SCIP_DATATREE *datatree = NULL;

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, 1) );
   SCIP_CALL( SCIPinsertDatatreeBool(scip, datatree, "test_bool", TRUE) );

   // Write to file
   FILE *file = fopen("test_write_bool.txt", "w");
   SCIP_CALL( SCIPwriteDatatreeJson(scip, file, datatree) );
   fclose(file); // Close the write file handle

   // Reopen file for reading and check content
   file = fopen("test_write_bool.txt", "r");
   const char* expected = "{\"test_bool\": true}";
   char buffer[1024];
   fgets(buffer, 1024, file);
   cr_assert_str_eq(buffer, expected, "File content is not as expected");
   fclose(file);

   // Cleanup
   remove("test_write_bool.txt");
   SCIPfreeDatatree(scip, &datatree);
}
