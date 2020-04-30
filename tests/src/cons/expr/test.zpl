var X1;
var X2;
var X3;
maximize obj: 25 * X1 + 30 * X2 + 10 * X3;
subto quad: X2^2 + 4*X1^2 == 0;
subto poly: X2^2 + 2 * X1 * X2^3 - X1^2 >= -0.2;
subto polyfun: 2 * (log(X2^4))^2 + 2 * X1 * X2^3 + X1^4 <= 1.0;


comment

This test problem has 3 constraints and 3 variables
the polyfun constraint is transformed by ZIMPL into:

<polyfun_1_a_0>: -@@polyfun_1_t_0 + X2^4 == 0.0
<polyfun_1_b_1>: 0.434294 * LN(@@polyfun_1_t_0) - @@polyfun_1_r_1 == 0.0
<polyfun_1>:     2 * @@polyfun_1_r_1^2 + 2 * X2^3 * X1 + X1^4 <= 1.0

Therefore we get 5 variables and 5 constraints, 3 of the
constraints representing reformulated constraint <polyfun>.
