your_version = 'GemmFiveLoops\_Pack\_12x4';
 
% number of repeats:% 3
% enter first, last, inc:% 48 960 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   960 4.8564e-02 3.6436e+01    6.5820e-02 2.6883e+01 9.6634e-13
   912 4.2014e-02 3.6109e+01    4.6898e-02 3.2349e+01 1.1084e-12
   864 3.5574e-02 3.6261e+01    4.1301e-02 3.1232e+01 8.8107e-13
   816 3.0215e-02 3.5964e+01    3.2978e-02 3.2952e+01 8.2423e-13
   768 2.4640e-02 3.6769e+01    2.6824e-02 3.3774e+01 7.3896e-13
   720 2.0327e-02 3.6725e+01    2.2370e-02 3.3371e+01 7.1054e-13
   672 1.6576e-02 3.6614e+01    1.8445e-02 3.2904e+01 6.5370e-13
   624 1.3368e-02 3.6351e+01    1.4787e-02 3.2862e+01 5.4001e-13
   576 1.0639e-02 3.5925e+01    1.1601e-02 3.2946e+01 5.1159e-13
   528 8.2871e-03 3.5525e+01    9.5720e-03 3.0756e+01 4.5475e-13
   480 6.1025e-03 3.6245e+01    6.7304e-03 3.2863e+01 3.2685e-13
   432 4.5016e-03 3.5819e+01    6.0965e-03 2.6448e+01 2.8422e-13
   384 3.2327e-03 3.5031e+01    3.4984e-03 3.2370e+01 2.1316e-13
   336 2.2062e-03 3.4387e+01    2.4111e-03 3.1465e+01 1.7053e-13
   288 1.4279e-03 3.3458e+01    1.6095e-03 2.9683e+01 1.1369e-13
   240 8.2190e-04 3.3639e+01    9.6837e-04 2.8551e+01 4.2633e-14
   192 4.3431e-04 3.2594e+01    4.9726e-04 2.8468e+01 2.8422e-14
   144 1.9714e-04 3.0293e+01    2.1242e-04 2.8114e+01 2.8422e-14
    96 6.6790e-05 2.6493e+01    6.4939e-05 2.7248e+01 1.4211e-14
    48 1.3770e-05 1.6063e+01    9.2090e-06 2.4018e+01 5.3291e-15
];

% Maximum difference between reference and your implementation: 1.108447e-12.
