your_version = 'JI\_4x4Kernel';
 
% number of repeats:% 3
% enter first, last, inc:% 48 960 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   960 4.7511e-02 3.7243e+01    2.0618e-01 8.5820e+00 9.6634e-13
   912 4.0973e-02 3.7026e+01    1.3544e-01 1.1201e+01 1.1084e-12
   864 3.4671e-02 3.7206e+01    1.3243e-01 9.7405e+00 8.8107e-13
   816 2.9460e-02 3.6886e+01    9.7767e-02 1.1115e+01 8.2423e-13
   768 2.4139e-02 3.7532e+01    1.0357e-01 8.7471e+00 7.3896e-13
   720 1.9680e-02 3.7931e+01    6.5976e-02 1.1315e+01 7.1054e-13
   672 1.6250e-02 3.7350e+01    6.6025e-02 9.1925e+00 6.5370e-13
   624 1.3084e-02 3.7141e+01    5.0887e-02 9.5494e+00 5.4001e-13
   576 1.0402e-02 3.6742e+01    4.1739e-02 9.1571e+00 5.1159e-13
   528 8.0984e-03 3.6352e+01    2.6444e-02 1.1133e+01 4.5475e-13
   480 5.9630e-03 3.7092e+01    2.2316e-02 9.9116e+00 3.2685e-13
   432 4.4042e-03 3.6612e+01    1.7935e-02 8.9905e+00 2.8422e-13
   384 3.1883e-03 3.5520e+01    1.2515e-02 9.0492e+00 2.1316e-13
   336 2.1588e-03 3.5143e+01    6.6367e-03 1.1431e+01 1.7053e-13
   288 1.4109e-03 3.3863e+01    5.2383e-03 9.1205e+00 1.1369e-13
   240 7.9288e-04 3.4870e+01    2.3969e-03 1.1535e+01 4.2633e-14
   192 4.2311e-04 3.3456e+01    1.1953e-03 1.1843e+01 2.8422e-14
   144 1.9378e-04 3.0818e+01    3.9628e-04 1.5070e+01 2.8422e-14
    96 6.6843e-05 2.6472e+01    1.1348e-04 1.5592e+01 1.4211e-14
    48 1.4002e-05 1.5797e+01    1.3534e-05 1.6343e+01 5.3291e-15
];

% Maximum difference between reference and your implementation: 1.108447e-12.
