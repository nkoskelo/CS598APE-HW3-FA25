# 598APE-HW3

## Reproducability Instructions

Use the same VM and docker image as the other assignments. 

### Timing tests
To build the executable use the command:
```
make all
```

Our test cases were
```
./main.exe 100 2.269 100000 64
```
```
./main.exe 1024 2.269 100000 64
```
```
./main.exe 1024 2.269 1000000000 64
```

Please checkout the following commits to see reproduce the times for those optimizations.

Initial - `3d0550609610e188a2ea1bbfd10b9345b8e73dc0`

Optimization 1 Predict $\Delta E$ - `3d0550609610e188a2ea1bbfd10b9345b8e73dc0`

Optimization 2 Reduce Data Storage Size - `83305cea551edf6f0b791d77ab53b8a8fff7cbef`

Optimization 3 Return $\Delta E$ directly - `88dd48dc02ae052be21034ce08460ac560b19518`

Optimization 4 Naive Parallelism - `8c6b7940ed3d1c518fc6dc71622d6ee4b9c967cd`

Optimization 5 Parallel Rand Without Load Balancing - `5b09082cc7e0cd91cf1b01c31f14f123fa596b69`

Optimization 6 Parallel Rand With Load Balancing - `fc4d7eeafd514e4c15a9bef1ffcbbb0e293ccf3e`

Optimization 7 Rand Seed False Sharing - 
`df639a32a226f36612aa113b300e33d9c71bde79`

Optimization 8 Partitioned Rows - 
`603d8da01142deca9f939f689fdbc618769994f0`

Optimization 9 Precompute Exponentials - 
`1c3e8ea0bb88107212c6bd9814c88c43f4d09700`

Optimization 10 Unroll Inner Loop - 
`a4241e4b31442f7750bf09d4c1ca751ce1fd18fa`


### Convergence Claims
To reproduce results about convergence speed, checkout branch `Convergence`. 
Convergence results can be obtained by saving the output of each execution to a file, and then
plotted alongside other executions with:
```
python3 plot_convergence.py file1 file2 file...
```
(This can be done on your local machine with the files downloaded. Requires `matplotlib`: `pip install matplotlib`)
