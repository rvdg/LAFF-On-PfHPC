your_version = 'GemmFiveLoops\_Pack\_12x4';
 
% number of repeats:% 3
% enter first, last, inc:% 48 960 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   960 4.7569e-02 3.7198e+01    5.3040e-02 3.3361e+01 9.6634e-13
   912 4.1035e-02 3.6971e+01    5.4781e-02 2.7694e+01 1.1084e-12
   864 3.4793e-02 3.7075e+01    4.6307e-02 2.7857e+01 8.8107e-13
   816 2.9416e-02 3.6942e+01    4.2473e-02 2.5585e+01 8.2423e-13
   768 2.4181e-02 3.7466e+01    3.2284e-02 2.8062e+01 7.3896e-13
   720 1.9926e-02 3.7463e+01    2.6920e-02 2.7730e+01 7.1054e-13
   672 1.6277e-02 3.7288e+01    2.1993e-02 2.7596e+01 6.5370e-13
   624 1.3132e-02 3.7005e+01    1.7550e-02 2.7688e+01 5.4001e-13
   576 1.0392e-02 3.6780e+01    1.4101e-02 2.7105e+01 5.1159e-13
   528 8.1048e-03 3.6323e+01    1.0779e-02 2.7311e+01 4.5475e-13
   480 5.9830e-03 3.6969e+01    8.0604e-03 2.7441e+01 3.2685e-13
   432 4.4097e-03 3.6565e+01    6.0802e-03 2.6519e+01 2.8422e-13
   384 3.1812e-03 3.5599e+01    4.2483e-03 2.6657e+01 2.1316e-13
   336 2.1646e-03 3.5049e+01    3.0722e-03 2.4694e+01 1.7053e-13
   288 1.4052e-03 3.3999e+01    1.8203e-03 2.6246e+01 1.1369e-13
   240 7.9691e-04 3.4694e+01    1.0867e-03 2.5442e+01 4.2633e-14
   192 4.2294e-04 3.3470e+01    5.5511e-04 2.5501e+01 2.8422e-14
   144 1.9465e-04 3.0680e+01    2.4663e-04 2.4214e+01 2.8422e-14
    96 6.7823e-05 2.6090e+01    7.5928e-05 2.3305e+01 1.4211e-14
    48 1.3809e-05 1.6017e+01    1.0105e-05 2.1889e+01 5.3291e-15
];

% Maximum difference between reference and your implementation: 1.108447e-12.
