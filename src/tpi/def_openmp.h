/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   def_openmp.h
 * @ingroup PARAINTERFACE
 * @brief  wrappers for OpenMP defines
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DEF_OPENMP_H__
#define __DEF_OPENMP_H__

#define STR(x)                #x
#define STRINGIFY(x)          STR(x)
#define CONCATENATE(x, y)     x y
#define CONCATPARENTH(x, y)   x ( y )

#define SPI_NONE  NULL

#define SPI_PRAGMA_CLAUSE(directive, clause)    _Pragma( STRINGIFY( CONCATENATE( directive, clause ) ) )
#define SPI_PRAGMA(directive)               _Pragma( STRINGIFY( directive ) )
#define SPI_PRAGMA_PARENTH(directive, var)  _Pragma( STRINGIFY( CONCATPARENTH( directive, var ) ) )



#define SPI_FOR_CLAUSE(clause)      SPI_PRAGMA_CLAUSE( SPI_DIR_FOR, clause )
#define SPI_FOR                     SPI_PRAGMA( SPI_DIR_FOR )
#define SPI_PARA_CLAUSE(clause)     SPI_PRAGMA_CLAUSE( SPI_DIR_PARA, clause )

#define SPI_PARA_CLAUSE_SHARED(priv, clause)       SPI_PRAGMA_CLAUSE( SPI_DIR_PARA,                           \
                                                      SPI_CLAUSE_DEFAULT( shared )                            \
                                                      SPI_CLAUSE_PRIVATE( (priv) ) clause )

#define SPI_PARA_SHARED                            SPI_PRAGMA_CLAUSE( SPI_DIR_PARA,                           \
                                                      SPI_CLAUSE_DEFAULT( shared ) )

#define SPI_PARA_SHARED_PRIVATE(priv)              SPI_PRAGMA_CLAUSE( SPI_DIR_PARA,                           \
                                                      SPI_CLAUSE_DEFAULT( shared )                            \
                                                      SPI_CLAUSE_PRIVATE( ( priv ) ) )

#define SPI_PARA_CLAUSE_NONE(share, priv, clause)  SPI_PRAGMA_CLAUSE( SPI_DIR_PARA,                           \
                                                      SPI_CLAUSE_DEFAULT( (none) )                            \
                                                      SPI_CLAUSE_SHARED( (share) )                            \
                                                      SPI_CLAUSE_PRIVATE( (priv) ) clause )

#define SPI_PARA                    SPI_PRAGMA( SPI_DIR_PARA )
#define SPI_CRITICAL(var)           SPI_PRAGMA_PARENTH( SPI_DIR_CRITICAL, var)
#define SPI_MASTER                  SPI_PRAGMA( SPI_DIR_MASTER )
#define SPI_WAIT                    SPI_PRAGMA( SPI_DIR_WAIT )
#define SPI_ORDERED                 SPI_PRAGMA( SPI_DIR_ORDERED )
#define SPI_SINGLE                  SPI_PRAGMA( SPI_DIR_SINGLE )
#define SPI_CLAUSE_SINGLE(clause)   SPI_PRAGMA_CLAUSE( SPI_DIR_SINGLE, clause )
#define SPI_TASK                    SPI_PRAGMA( SPI_DIR_TASK )
#define SPI_TASK_SHARED             SPI_PRAGMA_CLAUSE( SPI_DIR_TASK,                                           \
                                       SPI_CLAUSE_DEFAULT(shared) )
#define SPI_CLAUSE_TASK(clause)     SPI_PRAGMA_CLAUSE( SPI_DIR_TASK, clause )
#define SPI_TASKWAIT                SPI_PRAGMA( SPI_DIR_TASKWAIT )


/* OpenMP pragma directives */
#define SPI_DIR_PARA             omp parallel
#define SPI_DIR_FOR              omp for
#define SPI_DIR_CRITICAL         omp critical
#define SPI_DIR_MASTER           omp master
#define SPI_DIR_WAIT             omp barrier
#define SPI_DIR_ORDERED          omp ordered
#define SPI_DIR_TASK             omp task
#define SPI_DIR_SINGLE           omp single
#define SPI_DIR_TASKWAIT         omp taskwait


/* OpenMP clauses */
#define SPI_CLAUSE_PRIVATE(var)                 CONCATENATE( private, var )
#define SPI_CLAUSE_FSTPRIVATE(var)              CONCATENATE( firstprivate, var )
#define SPI_CLAUSE_LSTPRIVATE(var)              CONCATENATE( lastprivate, var )
#define SPI_CLAUSE_CPYPRIVATE(var)              CONCATENATE( copyprivate, var )
#define SPI_CLAUSE_NOWAIT                       nowait
#define SPI_CLAUSE_SHARED(var)                  CONCATENATE( shared, var )
#define SPI_CLAUSE_DEFAULT(var)                 CONCATPARENTH( default, var )
/* The reduce clause requires op as either an operator or intrinsic procedure.
 * Operators: +, *, .and., .or., .eqv., .neqv.
 * intrinsic procedures: max, min, iand, ior, or ieor*/
#define SPI_CLAUSE_REDUCE(op, var)              CONCATENATE( reduction, CONCATENATE( CONCATENATE( op, : ), var ) )
#define SPI_CLAUSE_ORDERED                      ordered
#define SPI_CLAUSE_IF(var)                      CONCATENATE( if, var )
#define SPI_CLAUSE_NUMTHREADS(var)              CONCATENATE( num_threads, var )
#define SPI_CLAUSE_SCHEDULE(type)               CONCATENATE( schedule, type )
#define SPI_CLAUSE_SCHEDULE_CHUNK(type, chunk)  CONCATENATE( schedule, CONCATPARENTH( type, chunk ) )
#define SPI_CLAUSE_COPYIN(var)                  CONCATENATE( copyin, var )
#define SPI_CLAUSE_FINAL(var)                   CONCATENATE( final, var )
#define SPI_CLAUSE_UNTIED                       untied
#define SPI_CLAUSE_MERGEABLE                    mergeable
#define SPI_CLAUSE_DEPEND(type, var)            CONCATENATE( depend, CONCATENATE( CONCATENATE( type, : ), var ) )
#define SPI_CLAUSE_PRIORITY(var)                CONCATENATE( priority, var )



#define SPI_SHARED_DATA(name, members)       struct SPI_Shared_Data {                                         \
                                                members                                                       \
                                             } name;

#define SPI_PRIVATE_DATA(name, members)      struct SPI_Private_Data {                                        \
                                                members                                                       \
                                             } name;

#define SPI_FSTPRIVATE_DATA(name, members)   struct SPI_FirstPrivate_Data {                                   \
                                                members                                                       \
                                             } name;

#define SPI_LSTPRIVATE_DATA(name, members)   struct SPI_LastPrivate_Data {                                    \
                                                members                                                       \
                                             } name;

#define SPI_COPYIN_DATA(name, members)       struct SPI_CopyIn_Data {                                         \
                                                members                                                       \
                                             } name;


#endif
