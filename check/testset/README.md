|      Testset               | n instances | Description |              Regular testing (see scip/scripts/jenkins/*_testruns.sh)              | Remarks |
|:--------------------------:|:-----------:|:-----------:|:----------------------------------------------------------------------------------:|:-------:|
| mipdebug                   |        4656 | -           | debug: weekdays 60s soplex (different settings) and cplex;                         | -       |
| MINLP                      |        1682 | -           | debug: weekdays 60s soplex (different settings) and cplex;                         | -       |
| mipdev-solvable            |         425 | -           | debug: weekends 7200s soplex; performance: saturdays 7200s 5 seeds exclusive M620v3; | -       |
| minlpdev-solvable          |         113 | -           | debug: weekends 7200s soplex; performance: saturdays 3600s 5 perms exclusive M640; | -       |
| sapdev-solvable            |         422 | -           | performance: sundays 3600s M630v2;                                                 | -       |
| miplib2017_benchmark       |         241 | -           | debug: fiberscip: sundays                                                          | -       |
| MMM                        |         168 | -           | debug: fiberscip: sundays                                                          | -       |
| clusterbench               |          15 | -           | fridays                                                                            | -       |
| allProblems                |        5182 | -           | -                                                                                  | -       |
| bugs                       |         140 | -           | -                                                                                  | -       |
| challenge                  |         164 | -           | -                                                                                  | -       |
| conflict                   |          86 | -           | -                                                                                  | -       |
| coral                      |         349 | -           | -                                                                                  | -       |
| coverage                   |          49 | -           | -                                                                                  | -       |
| fctp                       |          32 | -           | -                                                                                  | -       |
| indicator_zimpl            |           7 | -           | -                                                                                  | -       |
| infeasible                 |         185 | -           | -                                                                                  | -       |
| IP_0s_1s                   |         408 | -           | -                                                                                  | -       |
| IP_1h_INF                  |          36 | -           | -                                                                                  | -       |
| IP_1m_1h                   |         155 | -           | -                                                                                  | -       |
| IP_1s_1m                   |         187 | -           | -                                                                                  | -       |
| largeTestset               |        4820 | -           | -                                                                                  | -       |
| lotsize                    |          65 | -           | -                                                                                  | -       |
| mcf_ind_frac               |         133 | -           | -                                                                                  | -       |
| mcf_ind_int                |         132 | -           | -                                                                                  | -       |
| MINLP_0s_1m                |         604 | -           | -                                                                                  | -       |
| MINLP_0s_1s                |         381 | -           | -                                                                                  | -       |
| MINLP_boxqp                |          90 | -           | -                                                                                  | -       |
| minlpdev-complete          |         143 | -           | -                                                                                  | -       |
| minlpdev-hard              |          30 | -           | -                                                                                  | -       |
| MINLP_minlplib             |        1598 | -           | -                                                                                  | -       |
| MINLP_qplib                |         453 | -           | -                                                                                  | -       |
| MinULR-GuieuChinneck       |          20 | -           | -                                                                                  | -       |
| mipdev-1000seconds         |         101 | -           | -                                                                                  | -       |
| mipdev-100seconds          |          91 | -           | -                                                                                  | -       |
| mipdev-10seconds           |          99 | -           | -                                                                                  | -       |
| mipdev-7200seconds         |         130 | -           | -                                                                                  | -       |
| mipdev-complete            |         666 | -           | -                                                                                  | -       |
| mipdev-hard                |         241 | -           | -                                                                                  | -       |
| miplib2003                 |          60 | -           | -                                                                                  | -       |
| miplib2010_inf             |          19 | -           | -                                                                                  | -       |
| miplib2010                 |          87 | -           | -                                                                                  | -       |
| miplib3                    |          64 | -           | -                                                                                  | -       |
| miplib                     |          64 | -           | -                                                                                  | -       |
| mittelmann_current         |          75 | -           | -                                                                                  | -       |
| mittelmann_feas            |          32 | -           | -                                                                                  | -       |
| mittelmann_old             |          27 | -           | -                                                                                  | -       |
| mittelmann                 |          95 | -           | -                                                                                  | -       |
| mmm_0sec_to_10sec          |          48 | -           | -                                                                                  | -       |
| mmm_100sec_to_1h           |          36 | -           | -                                                                                  | -       |
| mmm_10sec_to_100sec        |          38 | -           | -                                                                                  | -       |
| MMMc                       |         496 | -           | -                                                                                  | -       |
| mmm_gt_1hour               |          31 | -           | -                                                                                  | -       |
| mmm_old                    |         163 | -           | -                                                                                  | -       |
| MMMshort                   |         123 | -           | -                                                                                  | -       |
| numerics                   |          27 | -           | -                                                                                  | -       |
| PBALL-DEC-SMALLINT-LIN     |         725 | -           | -                                                                                  | -       |
| PBALL-DEC-SMALLINT-NLC     |         100 | -           | -                                                                                  | -       |
| PBALL-OPT-SMALLINT-LIN     |        1994 | -           | -                                                                                  | -       |
| PBALL-OPT-SMALLINT-NLC     |         409 | -           | -                                                                                  | -       |
| PBALL-PARTIAL-SMALLINT-LIN |         536 | -           | -                                                                                  | -       |
| PBALL-PURE-SAT             |         334 | -           | -                                                                                  | -       |
| PBALL-SOFT-SMALLINT-LIN    |         201 | -           | -                                                                                  | -       |
| pure_subproblems           |          41 | -           | -                                                                                  | -       |
| sapdev-complete            |        4382 | -           | -                                                                                  | -       |
| sapdev-deco                |        4342 | -           | -                                                                                  | -       |
| sapdev-pure                |          40 | -           | -                                                                                  | -       |
| SAP-MMP                    |           4 | -           | -                                                                                  | -       |
| semicon                    |           4 | -           | -                                                                                  | -       |
| shortmiplib                |          23 | -           | -                                                                                  | -       |
| short                      |          46 | -           | -                                                                                  | -       |
| sos1                       |           3 | -           | -                                                                                  | -       |
| sos2                       |          13 | -           | -                                                                                  | -       |
| stableset                  |         119 | -           | -                                                                                  | -       |
| stochastic                 |           8 | -           | -                                                                                  | -       |
| stp                        |          75 | -           | -                                                                                  | -       |
| ufcn                       |          84 | -           | -                                                                                  | -       |
| xor                        |          31 | -           | -                                                                                  | -       |
| zib                        |          93 | -           | -                                                                                  | -       |
