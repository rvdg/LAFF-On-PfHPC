your_version = 'JI\_8x4Kernel';
 
% number of repeats:% 3
% enter first, last, inc:% 48 960 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   960 4.8442e-02 3.6527e+01    1.6525e-01 1.0708e+01 9.6634e-13
   912 4.1680e-02 3.6399e+01    1.1882e-01 1.2768e+01 1.1084e-12
   864 3.5411e-02 3.6428e+01    1.1038e-01 1.1686e+01 8.8107e-13
   816 2.9888e-02 3.6358e+01    8.5860e-02 1.2656e+01 8.2423e-13
   768 2.4642e-02 3.6765e+01    8.6491e-02 1.0475e+01 7.3896e-13
   720 2.0137e-02 3.7071e+01    5.8307e-02 1.2803e+01 7.1054e-13
   672 1.6560e-02 3.6650e+01    5.0955e-02 1.1911e+01 6.5370e-13
   624 1.3309e-02 3.6512e+01    4.0050e-02 1.2133e+01 5.4001e-13
   576 1.0599e-02 3.6061e+01    3.5485e-02 1.0771e+01 5.1159e-13
   528 8.2464e-03 3.5700e+01    2.4217e-02 1.2157e+01 4.5475e-13
   480 6.0624e-03 3.6485e+01    1.7807e-02 1.2421e+01 3.2685e-13
   432 4.4734e-03 3.6045e+01    1.3368e-02 1.2062e+01 2.8422e-13
   384 3.2550e-03 3.4791e+01    1.0269e-02 1.1028e+01 2.1316e-13
   336 2.1955e-03 3.4555e+01    6.0039e-03 1.2636e+01 1.7053e-13
   288 1.4346e-03 3.3302e+01    4.3100e-03 1.1085e+01 1.1369e-13
   240 8.0749e-04 3.4240e+01    2.1572e-03 1.2816e+01 4.2633e-14
   192 4.2791e-04 3.3081e+01    1.0430e-03 1.3572e+01 2.8422e-14
   144 1.9694e-04 3.0323e+01    2.4769e-04 2.4111e+01 2.8422e-14
    96 6.8574e-05 2.5804e+01    7.3328e-05 2.4131e+01 1.4211e-14
    48 1.4214e-05 1.5561e+01    8.2940e-06 2.6668e+01 5.3291e-15
];

% Maximum difference between reference and your implementation: 1.108447e-12.
