%nprocshared=4
%chk=N-Me-2-piridazone.chk
#p opt b3lyp/6-31g(d,p) nosymm pop=(npa,mk,hirshfeld)

Title Card Required

0 1
 C                 -3.14688800    1.79873400    0.00065300
 C                 -0.34036800    1.76702700    0.00012000
 O                  0.87682700    1.84695600   -0.00028600
 C                 -2.45699100    2.97909100    0.00077900
 H                 -2.94190400    3.95060900    0.00109100
 N                 -1.10292000    2.98324200    0.00052300
 C                 -0.33015400    4.22451700    0.00062800
 H                  0.31337200    4.25800800    0.88222700
 H                  0.31294300    4.25843700   -0.88126900
 H                 -1.01501400    5.07399300    0.00099900
 H                 -4.22831300    1.77544100    0.00087000
 N                 -2.36257100    0.61518500    0.00021000
 C                 -1.05393400    0.58037600   -0.00007200
 H                 -0.53963082   -0.35791612   -0.00044800

--Link1--
%nprocshared=4
%chk=N-Me-2-piridazone.chk
#p ub3lyp/6-31g(d,p) nosymm pop=(npa,mk,hirshfeld) SP scf=tight geom=checkpoint

Title Card Required

-1 2

--Link1--
%nprocshared=4
%chk=N-Me-2-piridazone.chk
#p ub3lyp/6-31g(d,p) nosymm pop=(npa,mk,hirshfeld) SP scf=tight geom=checkpoint

Title Card Required

1 2


