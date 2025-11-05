# 598APE-HW3

## Reproducability Instructions

Use the same VM and docker image as the other assignments. 


### Convergence Claims
To reproduce results about convergence speed, checkout brain `Convergence`. 
Convergence results can be obtained by saving the output of each execution to a file, and then
plotted alongside other executions with:
```
python3 plot_convergence.py file1 file2 file...
```
(This can be done on your local machine with the files downloaded. Requires `matplotlib`: `pip install matplotlib`)