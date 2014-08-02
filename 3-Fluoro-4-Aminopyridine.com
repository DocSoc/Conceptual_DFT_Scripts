%mem=200MW
%nproc=4
%chk=3-Fluoro-4-Aminopyridine
#p opt b3lyp/6-31+G(d) nosymm pop=(npa,hirshfeld)

3-Fluoro-4-Aminopyridine: Parent Compound Optimisation [B3LYP/6-31+G(d)]

0 1
C	-0.709	0.087	-0.062
N	2.125	-0.114	0.037
C	0.111	1.226	-0.037
C	-0.033	-1.142	-0.035
C	1.347	-1.206	0.015
C	1.494	1.067	0.010
H	-0.331	2.218	-0.060
H	1.835	-2.179	0.040
H	2.133	1.949	0.029
N	-2.085	0.126	-0.160
H	-2.566	-0.725	0.099
H	-2.539	0.968	0.166
F	-0.780	-2.276	-0.042

--link1--
%mem=200MW
%nproc=4
%chk=3-Fluoro-4-Aminopyridine
#p rob3lyp/6-31+G(d) scf=tight nosymm guess=read geom=checkpoint pop=(npa,hirshfeld)

3-Fluoro-4-Aminopyridine: Radical Anion Single Point.

-1 2

--link1--
%mem=200MW
%nproc=4
%chk=3-Fluoro-4-Aminopyridine
#p rob3lyp/6-31+G(d) scf=tight nosymm guess=read geom=checkpoint pop=(npa,hirshfeld)

3-Fluoro-4-Aminopyridine: Radical Cation Single Point.

1 2
