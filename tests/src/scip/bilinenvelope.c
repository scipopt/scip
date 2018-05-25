/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bilinenvelope.c
 * @brief  unit tests for computing envelopes of bilinear function
 * @author Benjamin Mueller
 */

#include "scip/scip.h"

#include "include/scip_test.h"

/* GLOBAL VARIABLES */
static SCIP* scip;
static SCIP_RANDNUMGEN* randnumgen;

/* TEST SUITE */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 0, TRUE) );
}

static
void teardown(void)
{
   SCIPfreeRandom(scip, &randnumgen);
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(bilinenvelope, .init = setup, .fini = teardown);

/* test xy for
 *    x <= y
 *    x in [0,1]
 *    y in [0,1]
 */
Test(bilinenvelope, gradcut1)
{
   SCIP_Real bilincoefs[6] = {1.0,   1.0,          1.0,             1.0,          -2.0,         -2.0};
   SCIP_Real refpointxs[6] = {0.5,   0.0,          0.02,            0.6,          0.5,          0.02};
   SCIP_Real refpointys[6] = {0.5,   0.0,          0.9,             0.4,          0.5,          0.9};
   SCIP_Bool successs[6]   = {TRUE,  FALSE,        TRUE,            FALSE,        FALSE,        TRUE};
   SCIP_Bool overests[6]   = {FALSE, FALSE,        FALSE,           FALSE,        FALSE,        TRUE};
   SCIP_Real coefxs[6]     = {0.75,  SCIP_INVALID, 0.3055555555555556,    SCIP_INVALID, SCIP_INVALID, -2.0 * 0.3055555555555556};
   SCIP_Real coefys[6]     = {0.25,  SCIP_INVALID, 0.027777777777777783,   SCIP_INVALID, SCIP_INVALID, -2.0 * 0.027777777777777783};
   SCIP_Real constants[6]  = {-0.25, SCIP_INVALID, -0.027777777777777783, SCIP_INVALID, SCIP_INVALID,  2.0 * 0.027777777777777783};
   SCIP_Real coefx;
   SCIP_Real coefy;
   SCIP_Real constant;
   SCIP_Bool success;
   int i;

   for( i = 0; i < 6; ++i )
   {
      /* reference point (0.5,0.5) */
      SCIPcomputeBilinEnvelope1(scip, bilincoefs[i], 0.0, 1.0, refpointxs[i], 0.0, 1.0, refpointys[i], overests[i], 1.0,
         1.0, 0.0, &coefx, &coefy, &constant, &success);

      cr_expect(success == successs[i], "%d: got %u expect %u", i, success, successs[i]);
      cr_expect(SCIPisEQ(scip, coefx, coefxs[i]));
      cr_expect(SCIPisEQ(scip, coefy, coefys[i]));
      cr_expect(SCIPisEQ(scip, constant, constants[i]));
   }
}

/* test xy for
 *    -x <= 2.16749 * y -1.30287
 *    x in [0.132186,0.434925]
 *    y in [0.471176,0.542986]
 */
Test(bilinenvelope, gradcut2)
{
   SCIP_Real bilincoefs[6] = {1.0, 1.0, 1.0, 1.0, -3.0, 1.0};
   SCIP_Real refpointxs[6] = {0.25, 0.2, 0.25, 0.25, 0.25, 0.35};
   SCIP_Real refpointys[6] = {0.5, 0.5088235701202773, 0.5, 0.5, 0.5, 0.5};
   SCIP_Bool successs[6]   = {TRUE, TRUE, TRUE, FALSE, TRUE, FALSE};
   SCIP_Bool overests[6]   = {TRUE, TRUE, TRUE, FALSE, FALSE, TRUE};
   SCIP_Real coefxs[6]     = {0.5269853305, 0.5347987912, 0.5269853305, SCIP_INVALID, -3.0 * 0.5269853305, SCIP_INVALID};
   SCIP_Real coefys[6]     = {0.2983037533, 0.2563010316, 0.2983037533, SCIP_INVALID, -3.0 * 0.2983037533, SCIP_INVALID};
   SCIP_Real constants[6]  = {-0.1550156705, -0.1356070501, -0.1550156705, SCIP_INVALID, 3.0 * 0.1550156705, SCIP_INVALID};
   SCIP_Real coefx;
   SCIP_Real coefy;
   SCIP_Real constant;
   SCIP_Bool success;
   int i;

   for( i = 0; i < 6; ++i )
   {
      SCIPcomputeBilinEnvelope1(scip, bilincoefs[i], 0.132186, 0.434925, refpointxs[i], 0.471176, 0.542986, refpointys[i],
	 overests[i], -1.0, 2.16749, -1.30287, &coefx, &coefy, &constant, &success);

      cr_expect(success == successs[i]);
      cr_expect(SCIPisEQ(scip, coefx, coefxs[i]));
      cr_expect(SCIPisEQ(scip, coefy, coefys[i]));
      cr_expect(SCIPisEQ(scip, constant, constants[i]));
   }
}

/* a numerically more challenging example */
Test(bilinenvelope, numerics)
{
   SCIP_Real lbx, ubx, lby, uby;
   SCIP_Real alpha, beta, gamma;
   SCIP_Real refpointx, refpointy;
   SCIP_Real coefx;
   SCIP_Real coefy;
   SCIP_Real constant;
   SCIP_Bool success;

   lbx = 1.3456657116093969;
   ubx = 5.9824198645533899;
   lby = 0.080398097539582644;
   uby = 0.35742536184127283;

   alpha = -1.0;
   beta = 16.737537240718815;
   gamma = -4.3942379318411717;

   /* note that the reference point needs to have a minimum distance to the variable bounds */
   refpointx = 5.2;
   refpointy = 0.3;

   SCIPcomputeBilinEnvelope1(scip, -325.08, lbx, ubx, refpointx, lby, uby, refpointy, FALSE, alpha, beta, gamma,
      &coefx, &coefy, &constant, &success);

   cr_expect(success);
   cr_expect(SCIPisEQ(scip, coefx, -71.50944001201063));
   cr_expect(SCIPisEQ(scip, coefy, -1449.1811847849826));
   cr_expect(SCIPisEQ(scip, constant, 250.66525223781193));
}

/* tests bilinear envelope computation for underestimating c*xy when given two linear inequalities */
Test(bilinenvelope, twoineqs_underestimate)
{
   SCIP_Real lbx = -1.0;
   SCIP_Real ubx = 2.0;
   SCIP_Real lby = -2.0;
   SCIP_Real uby = 3.0;
   SCIP_Real alpha1 = -1.2;
   SCIP_Real beta1 = -0.6;
   SCIP_Real gamma1 = 0.9;
   SCIP_Real alpha2 = 0.5;
   SCIP_Real beta2 = 0.3;
   SCIP_Real gamma2 = 0.6;
   int i;

   SCIP_Real bilincoefs[5] = {3.3, 100.3, -1.0, 1.0, 1.0};
   SCIP_Real refpointxs[5] = {0.3, -0.5, -0.5, 1.0, lbx};
   SCIP_Real refpointys[5] = {-0.15, -0.1, -0.1, 1.5, uby};
   SCIP_Real resx[5]       = {-0.0065211228254675, -0.7103545886154844, SCIP_INVALID, 1.45445115010332, SCIP_INVALID};
   SCIP_Real resy[5]       = {0.17701810527897127, -0.20848736066126078, SCIP_INVALID, 0.9772255750516624, SCIP_INVALID};
   SCIP_Real resconst[5]   = {-0.4315548420439882, -0.550126371865518, SCIP_INVALID, -1.9213268615442665, SCIP_INVALID};
   SCIP_Bool ressuccess[5] = {TRUE, TRUE, FALSE, TRUE, FALSE};

   for( i = 0; i < 5; ++i )
   {
      SCIP_Real constant;
      SCIP_Real coefx;
      SCIP_Real coefy;
      SCIP_Bool success;

      SCIPcomputeBilinEnvelope2(scip, bilincoefs[i], lbx, ubx, refpointxs[i], lby, uby, refpointys[i], FALSE,
            alpha1, beta1, gamma1, alpha2, beta2, gamma2, &coefx, &coefy, &constant, &success);

      /* check status */
      cr_expect(success == ressuccess[i], "%d: got %u expect %u\n", i, success, ressuccess[i]);

      /* check coefficients and the constant if successful */
      if( success )
      {
	 cr_expect(SCIPisEQ(scip, coefx, resx[i]*bilincoefs[i]));
	 cr_expect(SCIPisEQ(scip, coefy, resy[i]*bilincoefs[i]));
	 cr_expect(SCIPisEQ(scip, constant, resconst[i]*bilincoefs[i]));
      }
   }
}

/* tests bilinear envelope computation for overestimating c*xy when given two linear inequalities */
Test(bilinenvelope, twoineqs_overestimate)
{
   SCIP_Real lbx = -10.0;
   SCIP_Real ubx = 20.0;
   SCIP_Real lby = -20.0;
   SCIP_Real uby = 30.0;
   SCIP_Real alpha1 = 10.2;
   SCIP_Real beta1 = -6.0;
   SCIP_Real gamma1 = 120.0;
   SCIP_Real alpha2 = -5.0;
   SCIP_Real beta2 = 3.0;
   SCIP_Real gamma2 = 30;
   int i;

   SCIP_Real bilincoefs[4] = {3.0, 100.3, -1.0, 1.0};
   SCIP_Real refpointxs[4] = {10.0, -3.0, -3.0, lbx};
   SCIP_Real refpointys[4] = {-11.0, 15.0, 15.0, uby};
   SCIP_Real resx[4]       = {-11.49011102154311, 12.450432564701362, SCIP_INVALID, SCIP_INVALID};
   SCIP_Real resy[4]       = {9.752469181038805, -4.470333064628536, SCIP_INVALID, SCIP_INVALID};
   SCIP_Real resconst[4]   = {144.7533269821259, 89.40702892160698, SCIP_INVALID, SCIP_INVALID};
   SCIP_Bool ressuccess[4] = {TRUE, TRUE, FALSE, FALSE};

   for( i = 0; i < 4; ++i )
   {
      SCIP_Real constant;
      SCIP_Real coefx;
      SCIP_Real coefy;
      SCIP_Bool success;

      SCIPcomputeBilinEnvelope2(scip, bilincoefs[i], lbx, ubx, refpointxs[i], lby, uby, refpointys[i], TRUE,
            alpha1, beta1, gamma1, alpha2, beta2, gamma2, &coefx, &coefy, &constant, &success);

      /* check status */
      cr_expect(success == ressuccess[i]);

      /* check coefficients and the constant if successful */
      if( success )
      {
	 cr_expect(SCIPisEQ(scip, coefx, resx[i]*bilincoefs[i]));
	 cr_expect(SCIPisEQ(scip, coefy, resy[i]*bilincoefs[i]));
	 cr_expect(SCIPisEQ(scip, constant, resconst[i]*bilincoefs[i]));
      }
   }
}
