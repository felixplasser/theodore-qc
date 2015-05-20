1
 program ciudg      
 multireference single and double excitation configuration
 interaction based on the graphical unitary group approach.


 references:  h. lischka, r. shepard, f. b. brown, and i. shavitt,
                  int. j. quantum chem. s 15, 91 (1981).
              r. shepard, r. a. bair, r. a. eades, a. f. wagner,
                  m. j. davis, l. b. harding, and t. h. dunning,
                  int j. quantum chem. s 17, 613 (1983).
              r. ahlrichs, h.-j. boehm, c. ehrhardt, p. scharf,
                  h. schiffer, h. lischka, and m. schindler,
                  j. comp. chem. 6, 200 (1985).
              r. shepard, i. shavitt, r. m. pitzer, d. c. comeau, m. pepper
                  h. lischka, p. g. szalay, r. ahlrichs, f. b. brown, and
                  j.-g. zhao, int. j. quantum chem. symp. 22, 149 (1988).

 This Version of Program CIUDG is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              CIUDG       **
     **    PROGRAM VERSION:      2009-03.    **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************


================ Computing sorting integral file structure ================

                    -----z----- -----y----- -----x----- -----w----- ---total---

                CSFs        57           0           0           0          57
      internal walks        57           0           0           0          57
valid internal walks        57           0           0           0          57
  MR-CIS calculation - skip 3- and 4-external 2e integrals
 getinfoarray: info=                     6 :                     1
                  8192                  6552                  8192
                  5460                     0
 icd(3)=                  2461 ci%nnlev=                   820  l2rec=
                  8192  n2max=                  5460
 lcore1,lcore2=             917465227             917442361
 lencor,maxblo             917504000                 60000
========================================
 current settings:
 minbl3         543
 minbl4         543
 locmaxbl3     9104
 locmaxbuf     4552
 maxbl3       60000
 maxbl3       60000
 maxbl4       60000
 maxbuf       30006
========================================
 Orig.  diagonal integrals:  1electron:        40
                             0ext.    :       210
                             2ext.    :       728
                             4ext.    :       702


 Orig. off-diag. integrals:  4ext.    :     35217
                             3ext.    :     68250
                             2ext.    :     55419
                             1ext.    :     21840
                             0ext.    :      3444
                             2ext. SO :         0
                             1ext. SO :         0
                             0ext. SO :         0
                             1electron:       182


 Sorted integrals            3ext.  w :         0 x :         0
                             4ext.  w :         0 x :         0


Cycle #  1 sortfile size=    360437(      11 records of    32767) #buckets=   2
distributed memory consumption per node=         0 available core 917465227
Cycle #  2 sortfile size=    360437(      11 records of    32767) #buckets=   2
distributed memory consumption per node=         0 available core 917465227
 minimum size of srtscr:    294903 WP (     9 records)
 maximum size of srtscr:    360437 WP (    11 records)
 compressed index vector length=                     1
 echo of the input for program ciudg:
 ------------------------------------------------------------------------
  &input
  NTYPE = 0,
  GSET = 0,
   DAVCOR =10,
  NCOREL = 24
  NROOT = 3
  IVMODE = 3
  NBKITR = 1
  NVBKMN = 3
  RTOLBK = 1e-3,1e-3,1e-3,
  NITER = 60
  NVCIMN = 5
  RTOLCI = 1e-3,1e-3,1e-3,
  NVCIMX = 8
  NVRFMX = 8
  NVBKMX = 8
   iden=1
  CSFPRN = 10,
 /&end
 ------------------------------------------------------------------------
lodens (list->root)=  3
invlodens (root->list)= -1 -1  1
 USING SEGMENTS OF EQUAL SIZE

****************  list of control variables  ****************
 lvlprt =    0      nroot  =    3      noldv  =   0      noldhv =   0
 nunitv =    3      nbkitr =    1      niter  =  60      davcor =  10
 csfprn =   10      ivmode =    3      istrt  =   0      vout   =   0
 iortls =    0      nvbkmx =    8      ibktv  =  -1      ibkthv =  -1
 nvcimx =    8      icitv  =   -1      icithv =  -1      frcsub =   0
 nvbkmn =    3      nvcimn =    5      maxseg =   4      nrfitr =  30
 ncorel =   24      nvrfmx =    8      nvrfmn =   5      iden   =   1
 itran  =    0      froot  =    0      rtmode =   0      ncouple=   1
 skipso =    F      dalton2=    0      molcas =   0      finalv =   0
 finalw =    0      cosmocalc=   0    with_tsklst=   0
 nsegwx =    1     1     1     1
 nseg0x =    1     1     1     1
 nseg1x =    1     1     1     1
 nseg2x =    1     1     1     1
 nseg3x =    1     1     1     1
 nseg4x =    1     1     1     1
 no0ex  =      0    no1ex  =      0    no2ex  =     0    no3ex  =     0
 no4ex  =      0    nodiag =      0
 cdg4ex =    1      c3ex1ex=    1      c2ex0ex=   1
 fileloc=    0     0     0     0     0     0     0     1     1     1
 directhd=   1      noaqccshift_zyxw=      0
 critical_crit=-1.00000    critical_delta= 0.05000

 ctol   = 0.010000    lrtshift=1.000000    smalld =0.001000


 convergence tolerances of bk and full diagonalization steps
 root #       rtolbk        rtol
 ------      --------      ------
    1        1.000E-03    1.000E-03
    2        1.000E-03    1.000E-03
    3        1.000E-03    1.000E-03
 Computing density:                    .drt1.state3
 using                      1  nodes and                      1  cores.
 szdg/szodg per processor=                  3064                 82728
 Main memory management:
 global                1 DP per process
 vdisk                 0 DP per process
 stack                 0 DP per process
 core          917503999 DP per process

********** Integral sort section *************


 workspace allocation information: lencor= 917503999

 echo of the input for program cisrt:
 ------------------------------------------------------------------------
  &input
  maxbl3=60000
  maxbl4=60000
  &end
 ------------------------------------------------------------------------
 
 ( 6) listing file:                    ciudgls             
 ( 5) input file:                      cisrtin   
 (17) cidrt file:                      cidrtfl             
 (11) transformed integrals file:      moints    
 (12) diagonal integral file:          diagint             
 (13) off-diagonal integral file:      ofdgint             
 (31) 4-external w integrals file:     fil4w               
 (32) 4-external x integrals file:     fil4x               
 (33) 3-external w integrals file:     fil3w               
 (34) 3-external x integrals file:     fil3x               
 (21) scratch da sorting file:         srtscr              
 (12) 2-e integral file [fsplit=2]:    moints2   

 input integral file header information:
 Hermit Integral Program : SIFS version  vanadium.itc.univ 11:15:16.184 20-May-15
  cidrt_title                                                                    
 MO-coefficients from mcscf.x                                                    
  with dummy occupation 1.0 for active orbitals                                  
  total ao core energy =  100.585552274                                          
 MCSCF energy =    -227.500709909                                                
 SIFS file created by program tran.      vanadium.itc.univ 11:15:17.007 20-May-15

 input energy(*) values:
 energy( 1)=  1.005855522736E+02, ietype=   -1,    core energy of type: Nuc.Rep.

 total core energy =   1.005855522736E+02

 nsym = 2 nmot=  40

 symmetry  =    1    2
 slabel(*) =  A'   A" 
 nmpsy(*)  =   20   20

 info(*) =          1      8192      6552      8192      5460         0

 orbital labels, i:molab(i)=
   1:tout:001   2:tout:002   3:tout:003   4:tout:004   5:tout:005   6:tout:006   7:tout:007   8:tout:008   9:tout:009  10:tout:010
  11:tout:011  12:tout:012  13:tout:013  14:tout:014  15:tout:015  16:tout:016  17:tout:017  18:tout:018  19:tout:019  20:tout:020
  21:tout:021  22:tout:022  23:tout:023  24:tout:024  25:tout:025  26:tout:026  27:tout:027  28:tout:028  29:tout:029  30:tout:030
  31:tout:031  32:tout:032  33:tout:033  34:tout:034  35:tout:035  36:tout:036  37:tout:037  38:tout:038  39:tout:039  40:tout:040

 input parameters:
 prnopt=  0
 ldamin=    4095 ldamax=   32767 ldainc=      64
 maxbuf=   30006 maxbl3=   60000 maxbl4=   60000 intmxo=     766
  Using 32 bit compression 

 drt information:
  cidrt_title                                                                    
 nmotd =  44 nfctd =   4 nfvtc =   0 nmot  =  40
 nlevel =  40 niot  =  14 lowinl=  27
 orbital-to-level map(*)
   -1  -1  27  28  29  30  35  36  37   1   2   3   4   5   6   7   8   9  10  11
   12  13  -1  -1  31  32  33  34  38  39  40  14  15  16  17  18  19  20  21  22
   23  24  25  26
 compressed map(*)
   27  28  29  30  35  36  37   1   2   3   4   5   6   7   8   9  10  11  12  13
   31  32  33  34  38  39  40  14  15  16  17  18  19  20  21  22  23  24  25  26
 levsym(*)
    1   1   1   1   1   1   1   1   1   1   1   1   1   2   2   2   2   2   2   2
    2   2   2   2   2   2   1   1   1   1   2   2   2   2   1   1   1   2   2   2
 repartitioning mu(*)=
   2.  2.  2.  2.  2.  2.  2.  2.  0.  0.  0.  0.  0.  0.

 new core energy added to the energy(*) list.
 from the integral file: h1_core= -2.266432920827E+02

 indxdg: diagonal integral statistics.
 total number of integrals contributing to diagonal matrix elements:      1640
 number with all external indices:       702
 number with half external - half internal indices:       728
 number with all internal indices:       210

 indxof: off-diagonal integral statistics.
    4-external integrals: num=      35217 strt=          1
    3-external integrals: num=      68250 strt=      35218
    2-external integrals: num=      55419 strt=     103468
    1-external integrals: num=      21840 strt=     158887
    0-external integrals: num=       3444 strt=     180727

 total number of off-diagonal integrals:      184170


 indxof(2nd)  ittp=   3 numx(ittp)=       55419
 indxof(2nd)  ittp=   4 numx(ittp)=       21840
 indxof(2nd)  ittp=   5 numx(ittp)=        3444

 intermediate da file sorting parameters:
 nbuk=   2 lendar=   32767 nipbk=   21844 nipsg= 917328555
 pro2e        1     821    1641    2461    3281    3386    3491    4311   47999   91687
   124454  132646  138106  159945

 pro2e:    162205 integrals read in    30 records.
 pro1e        1     821    1641    2461    3281    3386    3491    4311   47999   91687
   124454  132646  138106  159945
 pro1e: eref =   -8.846144849560261E+01
 total size of srtscr:                    11  records of                  32767 
 WP =               2883496 Bytes

 new core energy added to the energy(*) list.
 from the hamiltonian repartitioning, eref= -8.846144849560E+01
 putdg        1     821    1641    2461    3227   35994   57839    4311   47999   91687
   124454  132646  138106  159945

 putf:       4 buffers of length     766 written to file 12
 diagonal integral file completed.
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd     55419         3    103468    103468     55419    184170
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd     21840         4    158887    158887     21840    184170
ptofdgf: num,ittp,ipos,istrtx,numx,maxrd      3444         5    180727    180727      3444    184170

 putf:     108 buffers of length     766 written to file 13
 off-diagonal files sort completed.
 executing brd_struct for cisrtinfo
cisrtinfo:
bufszi   766
 diagfile 4ext:     702 2ext:     728 0ext:     210
 fil4w,fil4x  :   35217 fil3w,fil3x :   68250
 ofdgint  2ext:   55419 1ext:   21840 0ext:    3444so0ext:       0so1ext:       0so2ext:       0
buffer minbl4     543 minbl3     543 maxbl2     546nbas:  13  13   0   0   0   0   0   0 maxbuf 30006
 CIUDG version 5.9.7 ( 5-Oct-2004)

 workspace allocation information: lcore= 917503999

 core energy values from the integral file:
 energy( 1)=  1.005855522736E+02, ietype=   -1,    core energy of type: Nuc.Rep.
 energy( 2)= -2.266432920827E+02, ietype=    6,   fcore energy of type: H1(*)   
 energy( 3)= -8.846144849560E+01, ietype=    5,   fcore energy of type: Vref(*) 

 total core repulsion energy = -2.145191883047E+02
 nmot  =    44 niot  =    14 nfct  =     4 nfvt  =     0
 nrow  =   103 nsym  =     2 ssym  =     1 lenbuf=  1600
 nwalk,xbar:         57       57        0        0        0
 nvalwt,nvalw:       57       57        0        0        0
 ncsft:              57
 total number of valid internal walks:      57
 nvalz,nvaly,nvalx,nvalw =       57       0       0       0

 cisrt info file parameters:
 file number  12 blocksize    766
 mxbld    766
 nd4ext,nd2ext,nd0ext   702   728   210
 n4ext,n3ext,n2ext,n1ext,n0ext,n2int,n1int,n0int    35217    68250    55419    21840     3444        0        0        0
 minbl4,minbl3,maxbl2   543   543   546
 maxbuf 30006
 number of external orbitals per symmetry block:  13  13
 nmsym   2 number of internal orbitals  14
 executing brd_struct for drt
 executing brd_struct for orbinf
 executing brd_struct for momap
 calcthrxt: niot,maxw1=                    14                     1
 block size     0
 pthz,pthy,pthx,pthw:    57     0     0     0 total internal walks:      57
 maxlp3,n2lp,n1lp,n0lp     1     0     0     0
 orbsym(*)= 1 1 1 1 2 2 2 2 1 1 1 2 2 2

 setref:       57 references kept,
                0 references were marked as invalid, out of
               57 total.
 nmb.of records onel     1
 nmb.of records 2-ext    73
 nmb.of records 1-ext    29
 nmb.of records 0-ext     5
 nmb.of records 2-int     0
 nmb.of records 1-int     0
 nmb.of records 0-int     0
 ---------memory usage in DP -----------------
 < n-ex core usage >
     routines:
    fourex            61129
    threx             60019
    twoex              3613
    onex               1157
    allin               766
    diagon             1343
               =======
   maximum            61129
 
  __ static summary __ 
   reflst                57
   hrfspc                57
               -------
   static->              57
 
  __ core required  __ 
   totstc                57
   max n-ex           61129
               -------
   totnec->           61186
 
  __ core available __ 
   totspc         917503999
   totnec -           61186
               -------
   totvec->       917442813

 number of external paths / symmetry
 vertex x     156     169
 vertex w     182     169
segment: free space=   917442813
 reducing frespc by                  1279 to              917441534 
  for index/conft/indsym storage .
 resegmenting ...



                   segmentation summary for type all-internal
 -------------------------------------------------------------------------------
 seg.      no. of|    no. of|  starting|  internal|  starting|  starting|
  no.    internal|        ci|       csf|     walks|      walk|       DRT|
            paths|  elements|    number|     /seg.|    number|    record|
 -------------------------------------------------------------------------------
  Z 1          57|        57|         0|        57|         0|         1|
 -------------------------------------------------------------------------------
max. additional memory requirements:index=           4DP  conft+indsym=         228DP  drtbuffer=        1047 DP

dimension of the ci-matrix ->>>        57

 executing brd_struct for civct
 gentasklist: ntask=                     2
                    TASKLIST
----------------------------------------------------------------------------------------------------
TASK# BRA# KET#  T-TYPE    DESCR.   SEGMENTTYPE    SEGEL              SEGCI          VWALKS   
----------------------------------------------------------------------------------------------------
     1  1   1     1      allint zz    OX  1 1      57      57         57         57      57      57
     2  1   1    75      dg-024ext z  OX  1 1      57      57         57         57      57      57
----------------------------------------------------------------------------------------------------
REDTASK #   1 TIME=   1.000 N=  1 (task/type/sgbra)=(   1/ 1/0) (
REDTASK #   2 TIME=   0.000 N=  1 (task/type/sgbra)=(   2/75/1) (
 initializing v-file: 1:                    57

    ---------trial vector generation----------

    trial vectors will be created by: 

    (ivmode= 3) diagonalizing h in the reference space.                     

      3 vectors will be written to unit 11 beginning with logical record   1

            3 vectors will be created
 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         888 2x:           0 4x:           0
All internal counts: zz :        1548 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 reference space has dimension      57
 dsyevx: computed roots 1 to    6(converged:   6)

    root           eigenvalues
    ----           ------------
       1        -227.6783581879
       2        -227.5324063198
       3        -227.4155723468
       4        -227.3853631405
       5        -227.3442327433
       6        -227.2711953660

 strefv generated    3 initial ci vector(s).
    ---------end of vector generation---------

 ufvoutnew: ... writing  recamt=                    57

         vector  1 from unit 11 written to unit 49 filename cirefv              
 ufvoutnew: ... writing  recamt=                    57

         vector  2 from unit 11 written to unit 49 filename cirefv              
 ufvoutnew: ... writing  recamt=                    57

         vector  3 from unit 11 written to unit 49 filename cirefv              

 ************************************************************************
 beginning the bk-type iterative procedure (nzcsf=    57)...
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:                57
 number of initial trial vectors:                         3
 number of initial matrix-vector products:                0
 maximum dimension of the subspace vectors:               8
 number of roots to converge:                             3
 number of iterations:                                    1
 residual norm convergence criteria:               0.001000  0.001000  0.001000

          starting bk iteration   1

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         888 2x:           0 4x:           0
All internal counts: zz :        1548 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         888 2x:           0 4x:           0
All internal counts: zz :        1548 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         888 2x:           0 4x:           0
All internal counts: zz :        1548 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3
   ht   1   -13.15916988
   ht   2    -0.00000000   -13.01321802
   ht   3    -0.00000000     0.00000000   -12.89638404

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3
 refs   1   1.000000        0.00000      -1.181374E-16
 refs   2    0.00000       1.000000       1.649093E-17
 refs   3    0.00000        0.00000       1.000000    

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3
 civs   1   1.000000      -2.536727E-14  -7.896196E-14
 civs   2  -2.560281E-14  -1.000000       2.780241E-14
 civs   3   7.907732E-14   2.781890E-14   1.000000    

          calcsovref: sovref block   1

              v      1       v      2       v      3
 ref    1   1.000000      -2.536727E-14  -7.908010E-14
 ref    2  -2.560281E-14  -1.000000       2.781890E-14
 ref    3   7.907732E-14   2.781890E-14   1.000000    

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3
 ref:   1     1.00000000    -0.00000000    -0.00000000
 ref:   2    -0.00000000    -1.00000000     0.00000000
 ref:   3     0.00000000     0.00000000     1.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -227.6783581879  3.1974E-14  6.4363E-27  7.4719E-14  1.0000E-03
 mr-sdci #  1  2   -227.5324063198  6.5725E-14  0.0000E+00  1.0568E-13  1.0000E-03
 mr-sdci #  1  3   -227.4155723468 -3.5527E-14  0.0000E+00  1.2992E-13  1.0000E-03
 
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000000
time for cinew                         0.000000
time for eigenvalue solver             0.000000
time for vector access                 0.000000

 mr-sdci  convergence criteria satisfied after  1 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -227.6783581879  3.1974E-14  6.4363E-27  7.4719E-14  1.0000E-03
 mr-sdci #  1  2   -227.5324063198  6.5725E-14  0.0000E+00  1.0568E-13  1.0000E-03
 mr-sdci #  1  3   -227.4155723468 -3.5527E-14  0.0000E+00  1.2992E-13  1.0000E-03
 
diagon:itrnv=   2
 expansion vectors are not transformed.
 matrix-vector products are not transformed.

    4 expansion eigenvectors written to unit nvfile (= 11)
    3 matrix-vector products written to unit nhvfil (= 10)

 ************************************************************************
 beginning the ci iterative diagonalization procedure... 
 ************************************************************************

               initial diagonalization conditions:

 number of configuration state functions:                57
 number of initial trial vectors:                         4
 number of initial matrix-vector products:                3
 maximum dimension of the subspace vectors:               8
 number of roots to converge:                             3
 number of iterations:                                   60
 residual norm convergence criteria:               0.001000  0.001000  0.001000

          starting ci iteration   1

 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         888 2x:           0 4x:           0
All internal counts: zz :        1548 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


 Final subspace hamiltonian 

                ht   1         ht   2         ht   3         ht   4
   ht   1   -13.15916988
   ht   2    -0.00000000   -13.01321802
   ht   3    -0.00000000     0.00000000   -12.89638404
   ht   4    -0.00000000    -0.00000000     0.00000000    -0.00000000

                v:   1         v:   2         v:   3         v:   4

   eig(s)   9.569831E-27    1.00000       1.000000       1.000000    
 
   x:   1  -4.832361E-14   3.249079E-02  -0.704735       0.708726    
   x:   2  -9.641580E-15   0.999185       5.910195E-03  -3.992966E-02
   x:   3   5.028423E-15  -2.395114E-02  -0.709446      -0.704353    
   x:   4   1.000000       1.132423E-14  -3.043098E-14   3.740499E-14
 bummer (warning):overlap matrix: # small eigenvalues=                      1

          calcsovref: tciref block   1

              civs   1       civs   2       civs   3       civs   4
 refs   1   1.000000        0.00000      -1.181374E-16   4.832361E-14
 refs   2    0.00000       1.000000       1.649093E-17   9.641580E-15
 refs   3    0.00000        0.00000       1.000000      -5.028423E-15

          calcsovref: scrb block   1

                ci   1         ci   2         ci   3         ci   4
 civs   1   -1.00000       1.374113E-14  -6.926282E-14  -0.493978    
 civs   2   1.122037E-14   -1.00000       3.269373E-14  -9.855895E-02
 civs   3  -7.247654E-14   3.146649E-14   1.000000       5.140195E-02
 civs   4   -1.01405      -0.315732       3.883434E-02   1.022228E+13

          calcsovref: sovref block   1

              v      1       v      2       v      3       v      4
 ref    1  -1.000000      -1.516181E-15  -6.750434E-14  -9.892087E-14
 ref    2   1.443290E-15  -1.000000       3.308465E-14  -3.108624E-14
 ref    3  -6.737745E-14   3.305412E-14   1.000000      -3.858025E-15

          reference overlap matrix  block   1

                ci   1         ci   2         ci   3         ci   4
 ref:   1    -1.00000000    -0.00000000    -0.00000000    -0.00000000
 ref:   2     0.00000000    -1.00000000     0.00000000    -0.00000000
 ref:   3    -0.00000000     0.00000000     1.00000000    -0.00000000

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -227.6783581879 -3.5527E-15  1.7507E-27  4.0059E-14  1.0000E-03
 mr-sdci #  1  2   -227.5324063198  0.0000E+00  0.0000E+00  1.0452E-13  1.0000E-03
 mr-sdci #  1  3   -227.4155723468 -3.5527E-15  0.0000E+00  1.3003E-13  1.0000E-03
 mr-sdci #  1  4   -226.9705698603  1.2451E+01  0.0000E+00  4.5899E-01  1.0000E-04
 
================ TIMING STATISTICS FOR JOB     ================
time for subspace matrix construction  0.000000
time for cinew                         0.000000
time for eigenvalue solver             0.000000
time for vector access                 0.000000

 mr-sdci  convergence criteria satisfied after  1 iterations.

 final mr-sdci  convergence information:

   iter      root         energy      deltae       apxde    residual       rtol
        ---- ----   --------------  ----------  ----------  ----------  ----------
 mr-sdci #  1  1   -227.6783581879 -3.5527E-15  1.7507E-27  4.0059E-14  1.0000E-03
 mr-sdci #  1  2   -227.5324063198  0.0000E+00  0.0000E+00  1.0452E-13  1.0000E-03
 mr-sdci #  1  3   -227.4155723468 -3.5527E-15  0.0000E+00  1.3003E-13  1.0000E-03
 mr-sdci #  1  4   -226.9705698603  1.2451E+01  0.0000E+00  4.5899E-01  1.0000E-04

####################CIUDGINFO####################

   ci vector at position   1 energy= -227.678358187948
   ci vector at position   2 energy= -227.532406319758
   ci vector at position   3 energy= -227.415572346828

################END OF CIUDGINFO################

 
diagon:itrnv=   0
    3 of the   5 expansion vectors are transformed.
    3 of the   4 matrix-vector products are transformed.

    3 expansion eigenvectors written to unit nvfile (= 11)
    3 matrix-vector products written to unit nhvfil (= 10)
maximum overlap with reference    1(overlap= 1.00000)

 information on vector: 1 from unit 11 written to unit 48 filename civout              
maximum overlap with reference    2(overlap= 1.00000)

 information on vector: 2 from unit 11 written to unit 48 filename civout              
maximum overlap with reference    3(overlap= 1.00000)

 information on vector: 3 from unit 11 written to unit 48 filename civout              


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 1) =      -227.6783581879

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7    8    9   10   11   12   13   14

                                          orbital     3    4    5    6   25   26   27   28    7    8    9   29   30   31

                                         symmetry   a'   a'   a'   a'   a"   a"   a"   a"   a'   a'   a'   a"   a"   a" 

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       1  0.010960                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-           
 z*  1  1       4  0.115906                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-      
 z*  1  1       7 -0.897431                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-   +-      
 z*  1  1       9  0.010065                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-        +- 
 z*  1  1      10 -0.242104                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +    +-    - 
 z*  1  1      12  0.080062                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-             +-   +- 
 z*  1  1      14  0.010811                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +-   +     - 
 z*  1  1      19 -0.017959                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +    +-    -    - 
 z*  1  1      22  0.011098                        +-   +-   +-   +-   +-   +-   +-   +-   +-        +-   +-   +-      
 z*  1  1      29 -0.232651                        +-   +-   +-   +-   +-   +-   +-   +-   +    +-    -   +-   +-      
 z*  1  1      32  0.052539                        +-   +-   +-   +-   +-   +-   +-   +-   +    +-    -   +    +-    - 
 z*  1  1      34  0.047557                        +-   +-   +-   +-   +-   +-   +-   +-   +    +-    -        +-   +- 
 z*  1  1      36 -0.187307                        +-   +-   +-   +-   +-   +-   +-   +-   +    +-   +     -   +-    - 
 z*  1  1      49  0.080204                        +-   +-   +-   +-   +-   +-   +-   +-        +-   +-   +-   +-      
 z*  1  1      52  0.051756                        +-   +-   +-   +-   +-   +-   +-   +-        +-   +-   +    +-    - 
 z*  1  1      54 -0.039589                        +-   +-   +-   +-   +-   +-   +-   +-        +-   +-        +-   +- 
 z*  1  1      55  0.101370                        +-   +-   +-   +-   +-   +-   +-   +-        +-        +-   +-   +- 

 ci coefficient statistics:
           rq > 0.1                6
      0.1> rq > 0.01              11
     0.01> rq > 0.001             15
    0.001> rq > 0.0001             1
   0.0001> rq > 0.00001            0
  0.00001> rq > 0.000001           0
 0.000001> rq                     24
           all                    57
 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         888 2x:           0 4x:           0
All internal counts: zz :        1548 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


  iref  icsf         v(icsf)             hv(icsf)
     1     1      0.010959598616     -2.495263419387
     2     2      0.000000008293     -0.000001888202
     3     3     -0.002871218754      0.653714371946
     4     4      0.115906334628    -26.389363971777
     5     5      0.000000004184     -0.000000952542
     6     6     -0.002274390333      0.517829456802
     7     7     -0.897431015040    204.325620091162
     8     8     -0.000000043675      0.000009943919
     9     9      0.010064951311     -2.291571589775
    10    10     -0.242103603576     55.121750973651
    11    11     -0.000000000887      0.000000201871
    12    12      0.080062449718    -18.228487104276
    13    13     -0.000000034679      0.000007895625
    14    14      0.010810632227     -2.461346996305
    15    15      0.000000000872     -0.000000198536
    16    16     -0.000000003555      0.000000809301
    17    17      0.001650594446     -0.375804633428
    18    18      0.000000004380     -0.000000997279
    19    19     -0.017959107412      4.088900089986
    20    20     -0.000000008770      0.000001996759
    21    21      0.002961303061     -0.674224618966
    22    22      0.011098075747     -2.526791665043
    23    23      0.000000000953     -0.000000217089
    24    24     -0.000463183086      0.105456764663
    25    25      0.002824439348     -0.643063713494
    26    26     -0.000000000192      0.000000043644
    27    27     -0.002253780808      0.513137114164
    28    28      0.009787877119     -2.228487792564
    29    29     -0.232651179725     52.969638630276
    30    30     -0.000000006155      0.000001401331
    31    31      0.002304942023     -0.524785415415
    32    32      0.052539048323    -11.962004262895
    33    33      0.000000004593     -0.000001045747
    34    34      0.047557030908    -10.827706717343
    35    35     -0.000000009281      0.000002113008
    36    36     -0.187307113496     42.645776077762
    37    37     -0.000000004742      0.000001079722
    38    38      0.000000001211     -0.000000275684
    39    39      0.001894170677     -0.431261669824
    40    40     -0.000000000186      0.000000042394
    41    41      0.000000004599     -0.000001047137
    42    42     -0.002180768744      0.496513847294
    43    43     -0.000000001802      0.000000410251
    44    44      0.000000006643     -0.000001512396
    45    45      0.003386642965     -0.771065309929
    46    46     -0.000000003389      0.000000771582
    47    47      0.003796559496     -0.864394432751
    48    48     -0.002341268019      0.533056058721
    49    49      0.080203763844    -18.260661272476
    50    50      0.000000005422     -0.000001234431
    51    51     -0.002127778992      0.484449227518
    52    52      0.051755895850    -11.783697393679
    53    53     -0.000000001399      0.000000318618
    54    54     -0.039589347237      9.013637580720
    55    55      0.101370381503    -23.079842029375
    56    56      0.000000002575     -0.000000586325
    57    57     -0.002119474231      0.482558413205

 number of reference csfs (nref) is    57.  root number (iroot) is  1.
 c0**2 =   1.00000000  c**2 (all zwalks) =   1.00000000

 pople ci energy extrapolation is computed with 24 correlated electrons.

 eref      =   -227.678358187948   "relaxed" cnot**2         =   1.000000000000
 eci       =   -227.678358187948   deltae = eci - eref       =   0.000000000000
 eci+dv1   =   -227.678358187948   dv1 = (1-cnot**2)*deltae  =   0.000000000000
 eci+dv2   =   -227.678358187948   dv2 = dv1 / cnot**2       =   0.000000000000
 eci+dv3   =   -227.678358187948   dv3 = dv1 / (2*cnot**2-1) =   0.000000000000
 eci+pople =   -227.678358187948   ( 24e- scaled deltae )    =   0.000000000000


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 2) =      -227.5324063198

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7    8    9   10   11   12   13   14

                                          orbital     3    4    5    6   25   26   27   28    7    8    9   29   30   31

                                         symmetry   a'   a'   a'   a'   a"   a"   a"   a"   a'   a'   a'   a"   a"   a" 

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       2 -0.139058                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -      
 z*  1  1       5 -0.056463                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +     - 
 z*  1  1       8  0.627784                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-   +     - 
 z*  1  1      13  0.687246                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +-   +-      
 z*  1  1      15 -0.014498                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +-        +- 
 z*  1  1      16  0.090248                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +    +-    - 
 z*  1  1      18 -0.082216                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -        +-   +- 
 z*  1  1      20  0.148666                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +     -   +-    - 
 z*  1  1      23 -0.015816                        +-   +-   +-   +-   +-   +-   +-   +-   +-        +-   +-   +     - 
 z*  1  1      30  0.072897                        +-   +-   +-   +-   +-   +-   +-   +-   +    +-    -   +-   +     - 
 z*  1  1      33 -0.071513                        +-   +-   +-   +-   +-   +-   +-   +-   +    +-    -   +     -   +- 
 z*  1  1      35  0.153071                        +-   +-   +-   +-   +-   +-   +-   +-   +    +-   +    +-    -    - 
 z*  1  1      37  0.064980                        +-   +-   +-   +-   +-   +-   +-   +-   +    +-   +     -    -   +- 
 z*  1  1      41 -0.080226                        +-   +-   +-   +-   +-   +-   +-   +-   +     -   +-   +    +-    - 
 z*  1  1      43  0.026244                        +-   +-   +-   +-   +-   +-   +-   +-   +     -   +-        +-   +- 
 z*  1  1      44 -0.112226                        +-   +-   +-   +-   +-   +-   +-   +-   +     -        +-   +-   +- 
 z*  1  1      46  0.072729                        +-   +-   +-   +-   +-   +-   +-   +-   +    +    +-    -   +-    - 
 z*  1  1      50 -0.080896                        +-   +-   +-   +-   +-   +-   +-   +-        +-   +-   +-   +     - 
 z*  1  1      53  0.027471                        +-   +-   +-   +-   +-   +-   +-   +-        +-   +-   +     -   +- 
 z*  1  1      56 -0.056246                        +-   +-   +-   +-   +-   +-   +-   +-        +     -   +-   +-   +- 

 ci coefficient statistics:
           rq > 0.1                6
      0.1> rq > 0.01              14
     0.01> rq > 0.001              3
    0.001> rq > 0.0001             1
   0.0001> rq > 0.00001            0
  0.00001> rq > 0.000001           0
 0.000001> rq                     33
           all                    57
 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         888 2x:           0 4x:           0
All internal counts: zz :        1548 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


  iref  icsf         v(icsf)             hv(icsf)
     1     1     -0.000000011458      0.000002607062
     2     2     -0.139058129599     31.640230846055
     3     3     -0.000000006287      0.000001430512
     4     4      0.000000012543     -0.000002854051
     5     5     -0.056462600808     12.847071428899
     6     6     -0.000000001846      0.000000419964
     7     7     -0.000000040824      0.000009288803
     8     8      0.627784304668   -142.841273490992
     9     9      0.000000021322     -0.000004851357
    10    10     -0.000000029288      0.000006664050
    11    11     -0.003965769267      0.902341024165
    12    12      0.000000004311     -0.000000980800
    13    13      0.687245591182   -156.370643094179
    14    14      0.000000047555     -0.000010820194
    15    15     -0.014497837203      3.298727785306
    16    16      0.090248384039    -20.534431986765
    17    17     -0.000000006846      0.000001557763
    18    18     -0.082216059862     18.706817938641
    19    19      0.000000016502     -0.000003754823
    20    20      0.148665921107    -33.826314767216
    21    21      0.000000006715     -0.000001527884
    22    22      0.000000016875     -0.000003839703
    23    23     -0.015815716154      3.598587954101
    24    24     -0.000000000391      0.000000089032
    25    25     -0.000000003468      0.000000789049
    26    26      0.003153016595     -0.717413453058
    27    27     -0.000000000227      0.000000051658
    28    28     -0.000000008570      0.000001949888
    29    29     -0.000000031115      0.000007079766
    30    30      0.072897075433    -16.586446986867
    31    31     -0.000000002187      0.000000497625
    32    32     -0.000000001387      0.000000315652
    33    33     -0.071513140806     16.271557011085
    34    34      0.000000008151     -0.000001854610
    35    35      0.153071492362    -34.828724996148
    36    36     -0.000000019177      0.000004363289
    37    37      0.064980297355    -14.785123420668
    38    38     -0.000935671307      0.212895543968
    39    39     -0.000000008126      0.000001848858
    40    40      0.003039336741     -0.691547602313
    41    41     -0.080226486882     18.254125610899
    42    42      0.000000000412     -0.000000093672
    43    43      0.026244497917     -5.971473763597
    44    44     -0.112226178499     25.535092446009
    45    45      0.000000006235     -0.000001418778
    46    46      0.072729411949    -16.548298111052
    47    47      0.000000002732     -0.000000621704
    48    48     -0.000000004616      0.000001050388
    49    49      0.000000005273     -0.000001199674
    50    50     -0.080896475217     18.406569668838
    51    51     -0.000000000859      0.000000195454
    52    52      0.000000008496     -0.000001933032
    53    53      0.027470988431     -6.250540101683
    54    54     -0.000000004892      0.000001113168
    55    55      0.000000009888     -0.000002249756
    56    56     -0.056245782374     12.797738208986
    57    57     -0.000000001287      0.000000292846

 number of reference csfs (nref) is    57.  root number (iroot) is  2.
 c0**2 =   1.00000000  c**2 (all zwalks) =   1.00000000

 pople ci energy extrapolation is computed with 24 correlated electrons.

 eref      =   -227.532406319758   "relaxed" cnot**2         =   1.000000000000
 eci       =   -227.532406319758   deltae = eci - eref       =  -0.000000000000
 eci+dv1   =   -227.532406319758   dv1 = (1-cnot**2)*deltae  =  -0.000000000000
 eci+dv2   =   -227.532406319758   dv2 = dv1 / cnot**2       =  -0.000000000000
 eci+dv3   =   -227.532406319758   dv3 = dv1 / (2*cnot**2-1) =  -0.000000000000
 eci+pople =   -227.532406319758   ( 24e- scaled deltae )    =   0.000000000000


 --- list of ci coefficients ( ctol =   1.00E-02 )  total energy( 3) =      -227.4155723468

                                                       internal orbitals

                                          level       1    2    3    4    5    6    7    8    9   10   11   12   13   14

                                          orbital     3    4    5    6   25   26   27   28    7    8    9   29   30   31

                                         symmetry   a'   a'   a'   a'   a"   a"   a"   a"   a'   a'   a'   a"   a"   a" 

 path  s ms    csf#    c(i)    ext. orb.(sym)
 z*  1  1       1 -0.483398                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-           
 z*  1  1       3 -0.086704                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +-   +          - 
 z*  1  1       9  0.385682                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +-        +-        +- 
 z*  1  1      14 -0.391227                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +-   +     - 
 z*  1  1      17  0.052078                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +     -   +     -   +- 
 z*  1  1      19 -0.225404                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +    +-    -    - 
 z*  1  1      21 -0.029624                        +-   +-   +-   +-   +-   +-   +-   +-   +-   +    +     -    -   +- 
 z*  1  1      22  0.479324                        +-   +-   +-   +-   +-   +-   +-   +-   +-        +-   +-   +-      
 z*  1  1      25 -0.086017                        +-   +-   +-   +-   +-   +-   +-   +-   +-        +-   +    +-    - 
 z*  1  1      28 -0.382836                        +-   +-   +-   +-   +-   +-   +-   +-   +-             +-   +-   +- 
 z*  1  1      31 -0.085875                        +-   +-   +-   +-   +-   +-   +-   +-   +    +-    -   +-        +- 
 z*  1  1      39  0.054154                        +-   +-   +-   +-   +-   +-   +-   +-   +     -   +-   +-   +     - 
 z*  1  1      45 -0.031793                        +-   +-   +-   +-   +-   +-   +-   +-   +    +    +-   +-    -    - 
 z*  1  1      48 -0.085070                        +-   +-   +-   +-   +-   +-   +-   +-   +          -   +-   +-   +- 
 z*  1  1      51  0.011966                        +-   +-   +-   +-   +-   +-   +-   +-        +-   +-   +-        +- 
 z*  1  1      57 -0.011871                        +-   +-   +-   +-   +-   +-   +-   +-             +-   +-   +-   +- 

 ci coefficient statistics:
           rq > 0.1                6
      0.1> rq > 0.01              10
     0.01> rq > 0.001              4
    0.001> rq > 0.0001             7
   0.0001> rq > 0.00001            3
  0.00001> rq > 0.000001           2
 0.000001> rq                     25
           all                    57
 ========= Executing OUT-OF-CORE method ========


====================================================================================================
Diagonal     counts:  0x:         888 2x:           0 4x:           0
All internal counts: zz :        1548 yy:           0 xx:           0 ww:           0
One-external counts: yz :           0 yx:           0 yw:           0
Two-external counts: yy :           0 ww:           0 xx:           0 xz:           0 wz:           0 wx:           0
Three-ext.   counts: yx :           0 yw:           0

SO-0ex       counts: zz :           0 yy:           0 xx:           0 ww:           0
SO-1ex       counts: yz :           0 yx:           0 yw:           0
SO-2ex       counts: yy :           0 xx:           0 wx:           0
====================================================================================================


  iref  icsf         v(icsf)             hv(icsf)
     1     1     -0.483398416101    109.932327469054
     2     2      0.000000012995     -0.000002955338
     3     3     -0.086704382245     19.717926713278
     4     4      0.000410842541     -0.093431991587
     5     5      0.000000002205     -0.000000501352
     6     6     -0.001895087193      0.430972338710
     7     7     -0.000021294934      0.004842799529
     8     8     -0.000000009560      0.000002174137
     9     9      0.385681798975    -87.710047057599
    10    10     -0.000095892840      0.021807525160
    11    11     -0.000000009102      0.000002069879
    12    12     -0.000338052505      0.076878403968
    13    13      0.000000006060     -0.000001378122
    14    14     -0.391226556217     88.971011199322
    15    15     -0.000000004905      0.000001115561
    16    16      0.000000005181     -0.000001178257
    17    17      0.052077949023    -11.843336583761
    18    18     -0.000000001176      0.000000267371
    19    19     -0.225404252201     51.260437023714
    20    20      0.000000003480     -0.000000791443
    21    21     -0.029624148780      6.736992750124
    22    22      0.479324076275   -109.005759145642
    23    23      0.000000006623     -0.000001506275
    24    24      0.000004604554     -0.001047147171
    25    25     -0.086017249615     19.561662052835
    26    26     -0.000000000923      0.000000209992
    27    27      0.001986708385     -0.451808424362
    28    28     -0.382835729399     87.062806515978
    29    29      0.000143384926     -0.032607964904
    30    30      0.000000007146     -0.000001625137
    31    31     -0.085875205034     19.529358903247
    32    32      0.000313752855     -0.071352285075
    33    33      0.000000001709     -0.000000388646
    34    34      0.000002809445     -0.000638911438
    35    35      0.000000002868     -0.000000652303
    36    36      0.000140458195     -0.031942380806
    37    37      0.000000000065     -0.000000014675
    38    38     -0.000000012618      0.000002869520
    39    39      0.054153809765    -12.315419642368
    40    40      0.000000000702     -0.000000159706
    41    41      0.000000001671     -0.000000380091
    42    42     -0.005583813729      1.269846195082
    43    43     -0.000000000038      0.000000008532
    44    44      0.000000009560     -0.000002174135
    45    45     -0.031792978733      7.230218455189
    46    46      0.000000000412     -0.000000093753
    47    47     -0.003296703484      0.749721709728
    48    48     -0.085070335440     19.346319023919
    49    49     -0.000212189720      0.048255246534
    50    50     -0.000000000843      0.000000191812
    51    51      0.011965686298     -2.721183398091
    52    52      0.000000306385     -0.000069676683
    53    53     -0.000000000315      0.000000071552
    54    54     -0.000011947552      0.002717059404
    55    55      0.000217451426     -0.049451840570
    56    56      0.000000001516     -0.000000344755
    57    57     -0.011870617472      2.699563266516

 number of reference csfs (nref) is    57.  root number (iroot) is  3.
 c0**2 =   1.00000000  c**2 (all zwalks) =   1.00000000

 pople ci energy extrapolation is computed with 24 correlated electrons.

 eref      =   -227.415572346828   "relaxed" cnot**2         =   1.000000000000
 eci       =   -227.415572346828   deltae = eci - eref       =  -0.000000000000
 eci+dv1   =   -227.415572346828   dv1 = (1-cnot**2)*deltae  =   0.000000000000
 eci+dv2   =   -227.415572346828   dv2 = dv1 / cnot**2       =   0.000000000000
 eci+dv3   =   -227.415572346828   dv3 = dv1 / (2*cnot**2-1) =   0.000000000000
 eci+pople =   -227.415572346828   ( 24e- scaled deltae )    =   0.000000000000
 passed aftci ... 
 readint2: molcas,dalton2=                     0                     0
 files%faoints=aoints              
lodens (list->root)=  1  2  3
                       Size (real*8) of d2temp for two-external contributions      51779
 
                       Size (real*8) of d2temp for all-internal contributions       3444
                       Size (real*8) of d2temp for one-external contributions      21840
                       Size (real*8) of d2temp for two-external contributions      51779
size_thrext:  lsym   l1    ksym   k1strt   k1       cnt3 
                1    1    1    1   14     1326
                1    1    2    1   14     3549
                1    2    1    1   14     1326
                1    2    2    1   14     3549
                1    3    1    1   14     1326
                1    3    2    1   14     3549
                1    4    1    1   14     1326
                1    4    2    1   14     3549
                2    5    2    1   14     4875
                2    6    2    1   14     4875
                2    7    2    1   14     4875
                2    8    2    1   14     4875
                1    9    1    1   14     1326
                1    9    2    1   14     3549
                1   10    1    1   14     1326
                1   10    2    1   14     3549
                1   11    1    1   14     1326
                1   11    2    1   14     3549
                2   12    2    1   14     4875
                2   13    2    1   14     4875
                2   14    2    1   14     4875
                       Size (real*8) of d2temp for three-external contributions      68250
                       Size (real*8) of d2temp for four-external contributions      35763
 serial operation: forcing vdisk for temporary dd012 matrices
location of d2temp files... fileloc(dd012)=       1
location of d2temp files... fileloc(d3)=       1
location of d2temp files... fileloc(d4)=       1
 files%dd012ext =  unit=  22  vdsk=   1  filestart=       1
 files%d3ext =     unit=  23  vdsk=   1  filestart=   78899
 files%d4ext =     unit=  24  vdsk=   1  filestart=  198899
            0xdiag    0ext      1ext      2ext      3ext      4ext
d2off                   767      4597     26811         1         1
d2rec                     5        29        68         2         1
recsize                 766       766       766     60000     60000
d2bufferlen=          60000
maxbl3=               60000
maxbl4=               60000
  allocated                 258899  DP for d2temp 
sifcfg setup: record length 4096 DP
# d1 elements per record  3272
# d2 elements per record  2730
  The MR-CISD density will be calculated.
 item #                     1 suffix=:.drt1.state1:
 read_civout: repnuc=  -214.519188304727     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref#100% root-following 0
 MR-CISD energy:  -227.67835819   -13.15916988
 residuum:     0.00000000
 deltae:    -0.00000000
 apxde:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -1.00000000    -0.00000000    -0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   2     0.00000000    -1.00000000     0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   3    -0.00000000     0.00000000     1.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -1.00000000    -0.00000000    -0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   2     0.00000000    -1.00000000     0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   3    -0.00000000     0.00000000     1.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    1
--------------------------------------------------------------------------------
================================================================================
   DYZ=       0  DYX=       0  DYW=       0
   D0Z=     144  D0Y=       0  D0X=       0  D0W=       0
  DDZI=     288 DDYI=       0 DDXI=       0 DDWI=       0
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
Trace of MO density:    24.000000
   24  correlated and     8  frozen core electrons

           D1 (irrep)  block   1

                MO   1         MO   2         MO   3         MO   4         MO   5         MO   6         MO   7         MO   8
   MO   1     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   2     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   3     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   4     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   5     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000
   MO   6     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000
   MO   7     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.86379522     0.00000001
   MO   8     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000001     1.99903045
   MO   9     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000    -0.25757951     0.00000005

                MO   9         MO  10         MO  11         MO  12         MO  13         MO  14         MO  15         MO  16
   MO   7    -0.25757951     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   8     0.00000005     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   9     0.14354649     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                MO  17         MO  18         MO  19         MO  20         MO  21         MO  22

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99903045     1.90153561
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.10580609     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

           D1 (irrep)  block   1

                MO   1         MO   2         MO   3         MO   4         MO   5         MO   6         MO   7         MO   8
   MO   1     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   2     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   3     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   4     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   5     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000
   MO   6     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000
   MO   7     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.85344895     0.00000001
   MO   8     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000001     1.99902510
   MO   9     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000    -0.26920031     0.00000006

                MO   9         MO  10         MO  11         MO  12         MO  13         MO  14         MO  15         MO  16
   MO   7    -0.26920031     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   8     0.00000006     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   9     0.14115379     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                MO  17         MO  18         MO  19         MO  20         MO  21         MO  22

Natural orbital populations,block 2
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99902510     1.89477421
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.09982853     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000


 total number of electrons =   32.0000000000

 test slabel:                     2
  a'  a" 


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a'  partial gross atomic populations
   ao class       1a'        2a'        3a'        4a'        5a'        6a' 
     3_ s       0.882842   0.268246   1.999829   0.000630   0.053888   0.000000
     3_ p       0.151100   0.129436   0.000007   0.000713   0.463295   1.089076
     4_ s       0.262726   1.459845  -0.000185   1.998648   0.215344   0.000000
     4_ p       0.208274   0.138935  -0.000578   0.000005   1.043039   0.129044
     4_ s       0.247529   0.001769   0.000464   0.000002   0.112217   0.390940
     4_ s       0.247529   0.001769   0.000464   0.000002   0.112217   0.390940
 
   ao class       7a'        8a'        9a'       10a'       11a'       12a' 
     3_ s       0.000000   0.000265  -0.000008   0.000000   0.000000   0.000000
     3_ p       0.009449   0.850334   0.059772   0.000000   0.000000   0.000000
     4_ s      -0.000000  -0.000057   0.000000   0.000000   0.000000   0.000000
     4_ p       1.776987   1.050843   0.046011   0.000000   0.000000   0.000000
     4_ s       0.106297   0.000075   0.000015   0.000000   0.000000   0.000000
     4_ s       0.106297   0.000075   0.000015   0.000000   0.000000   0.000000
 
   ao class      13a'       14a'       15a'       16a'       17a'       18a' 
 
   ao class      19a'       20a'       21a'       22a' 

                        a"  partial gross atomic populations
   ao class       1a"        2a"        3a"        4a"        5a"        6a" 
     3_ s       0.889376   0.258172   1.999309   0.000725   0.053802   0.000000
     3_ p       0.143137   0.125521   0.000006   0.000767   0.471388   1.097909
     4_ s       0.266831   1.471273   0.000013   1.998505   0.214572   0.000000
     4_ p       0.213648   0.140990  -0.000399   0.000005   1.032626   0.129341
     4_ s       0.243505   0.002023   0.000535  -0.000001   0.113806   0.386375
     4_ s       0.243505   0.002023   0.000535  -0.000001   0.113806   0.386375
 
   ao class       7a"        8a"        9a"       10a"       11a"       12a" 
     3_ s       0.000000   0.000351  -0.000012   0.000000   0.000000   0.000000
     3_ p       0.010169   0.803970   0.058312   0.000000   0.000000   0.000000
     4_ s       0.000000   0.000016   0.000001   0.000000   0.000000   0.000000
     4_ p       1.779966   1.090279   0.041502   0.000000   0.000000   0.000000
     4_ s       0.104445   0.000079   0.000013   0.000000   0.000000   0.000000
     4_ s       0.104445   0.000079   0.000013   0.000000   0.000000   0.000000
 
   ao class      13a"       14a"       15a"       16a"       17a"       18a" 
 
   ao class      19a"       20a"       21a"       22a" 


                        gross atomic populations
     ao            3_         4_         4_         4_
      s         6.407412   7.887532   1.710089   1.710089
      p         5.464361   8.820516   0.000000   0.000000
    total      11.871774  16.708049   1.710089   1.710089
 

 Total number of electrons:   32.00000000

 item #                     2 suffix=:.drt1.state2:
 read_civout: repnuc=  -214.519188304727     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref#100% root-following 0
 MR-CISD energy:  -227.67835819   -13.15916988
 residuum:     0.00000000
 deltae:    -0.00000000
 apxde:     0.00000000
================================================================================
  Reading record                      2  of civout
 INFO:ref#  2vector#  2 method:  0 last record  0max overlap with ref#100% root-following 0
 MR-CISD energy:  -227.53240632   -13.01321802
 residuum:     0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -1.00000000    -0.00000000    -0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   2     0.00000000    -1.00000000     0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   3    -0.00000000     0.00000000     1.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -1.00000000    -0.00000000    -0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   2     0.00000000    -1.00000000     0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   3    -0.00000000     0.00000000     1.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    2
--------------------------------------------------------------------------------
================================================================================
   DYZ=       0  DYX=       0  DYW=       0
   D0Z=     144  D0Y=       0  D0X=       0  D0W=       0
  DDZI=     288 DDYI=       0 DDXI=       0 DDWI=       0
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
Trace of MO density:    24.000000
   24  correlated and     8  frozen core electrons

           D1 (irrep)  block   1

                MO   1         MO   2         MO   3         MO   4         MO   5         MO   6         MO   7         MO   8
   MO   1     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   2     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   3     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   4     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   5     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000
   MO   6     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000
   MO   7     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.91597414    -0.00000002
   MO   8     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000    -0.00000002     1.46177430
   MO   9     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000    -0.02679282     0.00000001

                MO   9         MO  10         MO  11         MO  12         MO  13         MO  14         MO  15         MO  16
   MO   7    -0.02679282     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   8     0.00000001     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   9     0.63578533     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                MO  17         MO  18         MO  19         MO  20         MO  21         MO  22

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.91653464     1.46177430
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.63522483     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

           D1 (irrep)  block   1

                MO   1         MO   2         MO   3         MO   4         MO   5         MO   6         MO   7         MO   8
   MO   1     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   2     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   3     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   4     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   5     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000
   MO   6     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000
   MO   7     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.90730120    -0.00000001
   MO   8     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000    -0.00000001     1.53726669
   MO   9     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000    -0.04657561     0.00000002

                MO   9         MO  10         MO  11         MO  12         MO  13         MO  14         MO  15         MO  16
   MO   7    -0.04657561     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   8     0.00000002     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   9     0.54189833     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                MO  17         MO  18         MO  19         MO  20         MO  21         MO  22

Natural orbital populations,block 2
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.90888811     1.53726669
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.54031143     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000


 total number of electrons =   32.0000000000

 test slabel:                     2
  a'  a" 


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a'  partial gross atomic populations
   ao class       1a'        2a'        3a'        4a'        5a'        6a' 
     3_ s       0.882842   0.268246   1.999829   0.000630   0.053888   0.000000
     3_ p       0.151100   0.129436   0.000007   0.000713   0.463295   1.089076
     4_ s       0.262726   1.459845  -0.000185   1.998648   0.215344   0.000000
     4_ p       0.208274   0.138935  -0.000578   0.000005   1.043039   0.129044
     4_ s       0.247529   0.001769   0.000464   0.000002   0.112217   0.390940
     4_ s       0.247529   0.001769   0.000464   0.000002   0.112217   0.390940
 
   ao class       7a'        8a'        9a'       10a'       11a'       12a' 
     3_ s       0.000243   0.000000  -0.000040   0.000000   0.000000   0.000000
     3_ p       0.604936   0.006910   0.442410   0.000000   0.000000   0.000000
     4_ s      -0.000063  -0.000000   0.000004   0.000000   0.000000   0.000000
     4_ p       1.311331   1.299407   0.192645   0.000000   0.000000   0.000000
     4_ s       0.000044   0.077729   0.000102   0.000000   0.000000   0.000000
     4_ s       0.000044   0.077729   0.000102   0.000000   0.000000   0.000000
 
   ao class      13a'       14a'       15a'       16a'       17a'       18a' 
 
   ao class      19a'       20a'       21a'       22a' 

                        a"  partial gross atomic populations
   ao class       1a"        2a"        3a"        4a"        5a"        6a" 
     3_ s       0.889376   0.258172   1.999309   0.000725   0.053802   0.000000
     3_ p       0.143137   0.125521   0.000006   0.000767   0.471388   1.097909
     4_ s       0.266831   1.471273   0.000013   1.998505   0.214572   0.000000
     4_ p       0.213648   0.140990  -0.000399   0.000005   1.032626   0.129341
     4_ s       0.243505   0.002023   0.000535  -0.000001   0.113806   0.386375
     4_ s       0.243505   0.002023   0.000535  -0.000001   0.113806   0.386375
 
   ao class       7a"        8a"        9a"       10a"       11a"       12a" 
     3_ s       0.000313   0.000000  -0.000056   0.000000   0.000000   0.000000
     3_ p       0.574655   0.007820   0.382213   0.000000   0.000000   0.000000
     4_ s       0.000022   0.000000   0.000005   0.000000   0.000000   0.000000
     4_ p       1.333793   1.368809   0.157998   0.000000   0.000000   0.000000
     4_ s       0.000053   0.080319   0.000076   0.000000   0.000000   0.000000
     4_ s       0.000053   0.080319   0.000076   0.000000   0.000000   0.000000
 
   ao class      13a"       14a"       15a"       16a"       17a"       18a" 
 
   ao class      19a"       20a"       21a"       22a" 


                        gross atomic populations
     ao            3_         4_         4_         4_
      s         6.407277   7.887540   1.657487   1.657487
      p         5.691298   8.698910   0.000000   0.000000
    total      12.098576  16.586450   1.657487   1.657487
 

 Total number of electrons:   32.00000000

 item #                     3 suffix=:.drt1.state3:
 read_civout: repnuc=  -214.519188304727     
================================================================================
  Reading record                      1  of civout
 INFO:ref#  1vector#  1 method:  0 last record  0max overlap with ref#100% root-following 0
 MR-CISD energy:  -227.67835819   -13.15916988
 residuum:     0.00000000
 deltae:    -0.00000000
 apxde:     0.00000000
================================================================================
  Reading record                      2  of civout
 INFO:ref#  2vector#  2 method:  0 last record  0max overlap with ref#100% root-following 0
 MR-CISD energy:  -227.53240632   -13.01321802
 residuum:     0.00000000
================================================================================
  Reading record                      3  of civout
 INFO:ref#  3vector#  3 method:  0 last record  1max overlap with ref#100% root-following 0
 MR-CISD energy:  -227.41557235   -12.89638404
 residuum:     0.00000000
 deltae:    -0.00000000

          sovref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -1.00000000    -0.00000000    -0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   2     0.00000000    -1.00000000     0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   3    -0.00000000     0.00000000     1.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40

          tciref  block   1

                ci   1         ci   2         ci   3         ci   4         ci   5         ci   6         ci   7         ci   8
 ref:   1    -1.00000000    -0.00000000    -0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   2     0.00000000    -1.00000000     0.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 ref:   3    -0.00000000     0.00000000     1.00000000    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                ci   9         ci  10         ci  11         ci  12         ci  13         ci  14         ci  15         ci  16

                ci  17         ci  18         ci  19         ci  20         ci  21         ci  22         ci  23         ci  24

                ci  25         ci  26         ci  27         ci  28         ci  29         ci  30         ci  31         ci  32

                ci  33         ci  34         ci  35         ci  36         ci  37         ci  38         ci  39         ci  40
 computing final density
 =========== Executing IN-CORE method ==========
--------------------------------------------------------------------------------
  1e-density for root #    3
--------------------------------------------------------------------------------
================================================================================
   DYZ=       0  DYX=       0  DYW=       0
   D0Z=     144  D0Y=       0  D0X=       0  D0W=       0
  DDZI=     288 DDYI=       0 DDXI=       0 DDWI=       0
  DDZE=       0 DDYE=       0 DDXE=       0 DDWE=       0
================================================================================
Trace of MO density:    24.000000
   24  correlated and     8  frozen core electrons

           D1 (irrep)  block   1

                MO   1         MO   2         MO   3         MO   4         MO   5         MO   6         MO   7         MO   8
   MO   1     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   2     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   3     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   4     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   5     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000
   MO   6     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000
   MO   7     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.98083451    -0.00000002
   MO   8     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000    -0.00000002     1.00636842
   MO   9     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.12452006     0.00000000

                MO   9         MO  10         MO  11         MO  12         MO  13         MO  14         MO  15         MO  16
   MO   7     0.12452006     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   8     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   9     1.18730572     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                MO  17         MO  18         MO  19         MO  20         MO  21         MO  22

Natural orbital populations,block 1
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99991531     1.16822491
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     1.00636842     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

           D1 (irrep)  block   1

                MO   1         MO   2         MO   3         MO   4         MO   5         MO   6         MO   7         MO   8
   MO   1     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   2     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   3     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   4     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   5     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000     0.00000000
   MO   6     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     2.00000000     0.00000000     0.00000000
   MO   7     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     1.98143586    -0.00000002
   MO   8     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000    -0.00000002     0.99363263
   MO   9     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.14551463    -0.00000000

                MO   9         MO  10         MO  11         MO  12         MO  13         MO  14         MO  15         MO  16
   MO   7     0.14551463     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   8    -0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
   MO   9     0.85042286     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

                MO  17         MO  18         MO  19         MO  20         MO  21         MO  22

Natural orbital populations,block 2
              MO     1       MO     2       MO     3       MO     4       MO     5       MO     6       MO     7       MO     8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99985753     0.99363263
              MO     9       MO    10       MO    11       MO    12       MO    13       MO    14       MO    15       MO    16
  occ(*)=     0.83200119     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
              MO    17       MO    18       MO    19       MO    20       MO    21       MO    22       MO
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000


 total number of electrons =   32.0000000000

 test slabel:                     2
  a'  a" 


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        a'  partial gross atomic populations
   ao class       1a'        2a'        3a'        4a'        5a'        6a' 
     3_ s       0.882842   0.268246   1.999829   0.000630   0.000000   0.053888
     3_ p       0.151100   0.129436   0.000007   0.000713   1.089076   0.463295
     4_ s       0.262726   1.459845  -0.000185   1.998648   0.000000   0.215344
     4_ p       0.208274   0.138935  -0.000578   0.000005   0.129044   1.043039
     4_ s       0.247529   0.001769   0.000464   0.000002   0.390940   0.112217
     4_ s       0.247529   0.001769   0.000464   0.000002   0.390940   0.112217
 
   ao class       7a'        8a'        9a'       10a'       11a'       12a' 
     3_ s       0.000200  -0.000042  -0.000000   0.000000   0.000000   0.000000
     3_ p       0.308203   1.002331   0.004757   0.000000   0.000000   0.000000
     4_ s      -0.000069   0.000010   0.000000   0.000000   0.000000   0.000000
     4_ p       1.691554   0.165512   0.894585   0.000000   0.000000   0.000000
     4_ s       0.000014   0.000207   0.053513   0.000000   0.000000   0.000000
     4_ s       0.000014   0.000207   0.053513   0.000000   0.000000   0.000000
 
   ao class      13a'       14a'       15a'       16a'       17a'       18a' 
 
   ao class      19a'       20a'       21a'       22a' 

                        a"  partial gross atomic populations
   ao class       1a"        2a"        3a"        4a"        5a"        6a" 
     3_ s       0.889376   0.258172   1.999309   0.000725   0.053802   0.000000
     3_ p       0.143137   0.125521   0.000006   0.000767   0.471388   1.097909
     4_ s       0.266831   1.471273   0.000013   1.998505   0.214572   0.000000
     4_ p       0.213648   0.140990  -0.000399   0.000005   1.032626   0.129341
     4_ s       0.243505   0.002023   0.000535  -0.000001   0.113806   0.386375
     4_ s       0.243505   0.002023   0.000535  -0.000001   0.113806   0.386375
 
   ao class       7a"        8a"        9a"       10a"       11a"       12a" 
     3_ s       0.000248   0.000000  -0.000053   0.000000   0.000000   0.000000
     3_ p       0.306948   0.005054   0.711320   0.000000   0.000000   0.000000
     4_ s       0.000031   0.000000   0.000004   0.000000   0.000000   0.000000
     4_ p       1.692576   0.884748   0.120473   0.000000   0.000000   0.000000
     4_ s       0.000027   0.051915   0.000129   0.000000   0.000000   0.000000
     4_ s       0.000027   0.051915   0.000129   0.000000   0.000000   0.000000
 
   ao class      13a"       14a"       15a"       16a"       17a"       18a" 
 
   ao class      19a"       20a"       21a"       22a" 


                        gross atomic populations
     ao            3_         4_         4_         4_
      s         6.407170   7.887547   1.604969   1.604969
      p         6.010968   8.484376   0.000000   0.000000
    total      12.418138  16.371923   1.604969   1.604969
 

 Total number of electrons:   32.00000000

 DA ...
