Generating sequences:
> python3 Assignment12.py -n toy.nwk -m 0.1 -l 100
Output will be written to Assignment12.fasta


Test of proper parsing:
> python3 JC_Tree.py 
root : 0
        node 2 : 0.5
                node 0 : 0.25
                        f1 : 0.25
                        f2 : 0.25
                node 1 : 0.45
                        f3 : 0.05
                        f4 : 0.05
        node 5 : 0.4
                node 3 : 0.2
                        f5 : 0.4
                        f6 : 0.4
                node 4 : 0.3
                        f7 : 0.3
                        f8 : 0.3

Example Distance Matrices:
['> f8', '> f7', '> f6', '> f5', '> f4', '> f3', '> f2', '> f1']
[[ 0.  9. 16. 10. 22. 23. 25. 24.]
 [ 0.  0. 18. 12. 24. 25. 26. 24.]
 [ 0.  0.  0.  8. 25. 26. 26. 27.]
 [ 0.  0.  0.  0. 20. 21. 22. 22.]
 [ 0.  0.  0.  0.  0.  1. 12. 11.]
 [ 0.  0.  0.  0.  0.  0. 13. 12.]
 [ 0.  0.  0.  0.  0.  0.  0.  5.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.]]
[[-0.      0.0959  0.18    0.1073  0.2604  0.2747  0.3041  0.2892]
 [-0.     -0.      0.2058  0.1308  0.2892  0.3041  0.3193  0.2892]
 [-0.     -0.     -0.      0.0846  0.3041  0.3193  0.3193  0.3347]
 [-0.     -0.     -0.     -0.      0.2326  0.2464  0.2604  0.2604]
 [-0.     -0.     -0.     -0.     -0.      0.0101  0.1308  0.119 ]
 [-0.     -0.     -0.     -0.     -0.     -0.      0.1428  0.1308]
 [-0.     -0.     -0.     -0.     -0.     -0.     -0.      0.0517]
 [-0.     -0.     -0.     -0.     -0.     -0.     -0.     -0.    ]]

 The matrices are as expected. We can clearly see the higher distances between 1, 2, 3 & 4 and 5, 6, 7 & 8.
 The distance between f3 and f4 are as expected very low, between f5 and f6 still pretty low, between f1 and f5,f6,f7,f8 high.