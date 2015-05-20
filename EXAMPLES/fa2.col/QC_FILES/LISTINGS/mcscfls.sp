

     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.5         **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 This program allows the csf mixing coefficient and orbital expansion coefficient
 optimization using the graphical unitary group approach and the exponential
 operator mcscf method.
 references:  r. shepard and j. simons, ' int. j. quantum chem. symp. 14, 211 (1980).
              r. shepard, i. shavitt, and j. simons, j. chem. phys. 76, 543 (1982).
              r. shepard in "ab initio methods in quantum chemistry ii" advances in chemical
                  physics 69, edited by k. p. lawley (wiley, new york, 1987) pp. 63-200.
 Original autor: Ron Shepard, ANL
 Later revisions: Michal Dallos, University Vienna

 This Version of Program MCSCF is Maintained by:
     Thomas Mueller
     Juelich Supercomputing Centre (JSC)
     Institute of Advanced Simulation (IAS)
     D-52425 Juelich, Germany 
     Email: th.mueller@fz-juelich.de



     ******************************************
     **    PROGRAM:              MCSCF       **
     **    PROGRAM VERSION:      5.4.0.2     **
     **    DISTRIBUTION VERSION: 5.9.a       **
     ******************************************

 Workspace allocation information:
       917504000 of real*8 words ( 7000.00 MB) of work space has been allocated.

 user input information:

 ======== echo of the mcscf input ========
 ------------------------------------------------------------------------
  &input
   niter=100,
   nmiter=50,
   nciitr=300,
   tol(3)=1.e-4,
   tol(2)=1.e-4,
   tol(1)=1.e-8,
   NSTATE=0,
   npath=1,3,9,10,13,17,19,21,-11,12, 2, 34,30
   ncoupl=5,
   tol(9)=1.e-3,
   FCIORB=  1,7,40,1,8,40,1,9,40,2,7,40,2,8,40,2,9,40
   NAVST(1) = 3,
   WAVST(1,1)=1 ,
   WAVST(1,2)=1 ,
   WAVST(1,3)=1 ,
   NAVST(2) = 2,
   WAVST(2,1)=1 ,
   WAVST(2,2)=1 ,
  &end
 ------------------------------------------------------------------------


 ***  Integral file informations  ***


 input integral file : /public/vanadium/scratch/plasserf/temp/WORK/aoints       
    

 Integral file header information:
 Hermit Integral Program : SIFS version  vanadium.itc.univ 11:15:16.184 20-May-15

 Core type energy values:
 energy( 1)=  1.005855522736E+02, ietype=   -1,    core energy of type: Nuc.Rep.
 total ao core energy =  100.585552274


   ******  Basis set information:  ******

 Number of irreps:                  2
 Total number of basis functions:  44

 irrep no.              1    2
 irrep label           A'   A" 
 no. of bas.fcions.    22   22


 ***  MCSCF optimization procedure parmeters:  ***


 maximum number of mcscf iterations:        niter=   100

 maximum number of psci micro-iterations:   nmiter=   50
 maximum r,s subspace dimension allowed:    nvrsmx=   30

 tol(1)=  1.0000E-08. . . . delta-emc convergence criterion.
 tol(2)=  1.0000E-04. . . . wnorm convergence criterion.
 tol(3)=  1.0000E-04. . . . knorm convergence criterion.
 tol(4)=  1.0000E-08. . . . apxde convergence criterion.
 tol(5)=  1.0000E-04. . . . small diagonal matrix element tolerance.
 tol(6)=  1.0000E-06. . . . minimum ci-psci residual norm.
 tol(7)=  1.0000E-05. . . . maximum ci-psci residual norm.
 tol(8)=  1.0000E+00. . . . maximum abs(k(xy)) allowed.
 tol(9)=  1.0000E-03. . . . wnorm coupling tolerance.
 tol(10)= 0.0000E+00. . . . maximum psci emergency shift parameter.
 tol(11)= 0.0000E+00. . . . minimum psci emergency shift parameter.
 tol(12)= 0.0000E+00. . . . increment of psci emergency shift parameter.


 *** State averaging informations: ***


 MCSCF calculation performed for  2 DRTs.

 DRT  first state   no.of aver.states   weights
  1   ground state          3             0.200 0.200 0.200
  2   ground state          2             0.200 0.200

 The number of hmc(*) eigenvalues and eigenvectors calculated each iteration per DRT:
 DRT.   no.of eigenv.(=ncol)
    1        4
    2        3

 Orbitals included in invariant subspaces:
   symmetry   orbital   mask
       1       7(  7)    40
       1       8(  8)    40
       1       9(  9)    40
       2       7( 29)    40
       2       8( 30)    40
       2       9( 31)    40

 npath(*) options:
  2:  orbital-state coupling terms will be included beginning on iteration ncoupl=  5
  3:  print intermediate timing information.
  9:  suppress the drt listing.
 10:  suppress the hmc(*) eigenvector listing.
 12:  diagonalize the hmc(*) matrix iteratively.
        nunitv= 1 nciitr=** mxvadd=20 nvcimx=20
       rtolci(*),wnorm=     1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02 1.0000E-02
                            1.0000E-02 1.0000E-02 0.0000E+00
   noldv =   0  0
 13:  get initial orbitals from the formatted file, mocoef.
 17:  print the final natural orbitals and occupations.
 19:  transform the virtual orbitals to diagonalize qvv(*).
 21:  write out the one- and two- electron density for further use (files:mcd1fl, mcd2fl).
 30:  Compute mcscf (transition) density matrices


   ******  DRT info section  ******


 Informations for the DRT no.  1

 DRT file header:
  title                                                                          
 Molecular symmetry group:    a' 
 Total number of electrons:   32
 Spin multiplicity:            1
 Number of active orbitals:    6
 Number of active electrons:   8
 Total number of CSFs:        57

 Informations for the DRT no.  2

 DRT file header:
  title                                                                          
 Molecular symmetry group:    a" 
 Total number of electrons:   32
 Spin multiplicity:            1
 Number of active orbitals:    6
 Number of active electrons:   8
 Total number of CSFs:        48
 

 faar:   0 active-active rotations allowed out of:   6 possible.


 Number of active-double rotations:        36
 Number of active-active rotations:         0
 Number of double-virtual rotations:      156
 Number of active-virtual rotations:       78
 lenbfsdef=                131071  lenbfs=                   338
  number of integrals per class 1:11 (cf adda 
 class  1 (pq|rs):         #         123
 class  2 (pq|ri):         #         756
 class  3 (pq|ia):         #        3276
 class  4 (pi|qa):         #        5616
 class  5 (pq|ra):         #        1638
 class  6 (pq|ij)/(pi|qj): #        2340
 class  7 (pq|ab):         #        3705
 class  8 (pa|qb):         #        7098
 class  9 p(bp,ai)         #       12168
 class 10p(ai,jp):        #        5616
 class 11p(ai,bj):        #       13182

 Size of orbital-Hessian matrix B:                    38079
 Size of the orbital-state Hessian matrix C:          72090
 Total size of the state Hessian matrix M:                0
 Size of HESSIAN-matrix for quadratic conv.:         110169


 Source of the initial MO coeficients:

 Input MO coefficient file: /public/vanadium/scratch/plasserf/temp/WORK/mocoef          
 

               starting mcscf iteration...   1

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     2, naopsy(1) =    22, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 917497823

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 917330879
 address segment size,           sizesg = 917179393
 number of in-core blocks,       nincbk =         4
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:      91275 transformed 1/r12    array elements were written in      17 records.


 mosort: allocated sort2 space, avc2is=   917359676 available sort2 space, avcisx=   917359928

 trial vectors are generated internally.

 trial vector  1 is unit matrix column     7
 ciiter=  35 noldhv= 10 noldv= 10

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.6647233699     -328.2502756435        0.0000001551        0.0000010000
    2      -227.4760263671     -328.0615786407        0.0000004365        0.0000010000
    3      -227.3257747168     -327.9113269904        0.0000009058        0.0000010000
    4      -227.2857391078     -327.8712913814        0.0016340608        0.0100000000

 trial vectors are generated internally.

 trial vector  1 is unit matrix column     1
 ciiter=  15 noldhv=  9 noldv=  9

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.4769067603     -328.0624590339        0.0000002396        0.0000010000
    2      -227.2735994919     -327.8591517655        0.0000005610        0.0000010000
    3      -227.2200344057     -327.8055866793        0.0001759608        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  4.858892827992890E-002
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.93252643 pnorm= 0.0000E+00 rznorm= 6.8490E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -41.354450  -22.629446   -2.939532   -1.745191   -1.399367   -1.350899
 i,qaaresolved                     1  -1.13091884388698     
 i,qaaresolved                     2 -0.722780135308617     
 i,qaaresolved                     3 -1.348741110696561E-002

 qvv(*) eigenvalues. symmetry block  1
     0.466660    0.638311    0.685309    1.342354    1.583549    1.674826    1.982560    2.265821    2.268428    2.378083
     2.425038    2.694710    3.721831

 fdd(*) eigenvalues. symmetry block  2
   -41.354457  -22.629470   -2.937787   -1.737187   -1.394294   -1.347034
 i,qaaresolved                     1  -1.09559914858769     
 i,qaaresolved                     2 -0.717873094133593     
 i,qaaresolved                     3  9.617597945215772E-002

 qvv(*) eigenvalues. symmetry block  2
     0.516892    0.687769    0.796341    1.589895    1.715331    1.947172    2.136414    2.317494    2.363153    2.450049
     2.571401    2.760035    3.840972

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    1 emc=   -227.4434061412 demc= 2.2744E+02 wnorm= 3.8871E-01 knorm= 3.6110E-01 apxde= 4.4576E-02    *not conv.*     

               starting mcscf iteration...   2

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     2, naopsy(1) =    22, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 917497823

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 917330879
 address segment size,           sizesg = 917179393
 number of in-core blocks,       nincbk =         4
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     128342 transformed 1/r12    array elements were written in      24 records.


 mosort: allocated sort2 space, avc2is=   917359676 available sort2 space, avcisx=   917359928

   4 trial vectors read from nvfile (unit= 29).
 ciiter=  15 noldhv=  7 noldv=  7

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.6778758228     -328.2634280964        0.0000010371        0.0000100000
    2      -227.5280877959     -328.1136400695        0.0000088762        0.0000100000
    3      -227.4106367311     -327.9961890047        0.0000036918        0.0000100000
    4      -227.3771745486     -327.9627268222        0.0022271244        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   8 noldhv= 19 noldv= 19

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.5286085905     -328.1141608641        0.0000057295        0.0000100000
    2      -227.3407414222     -327.9262936958        0.0000094670        0.0000100000
    3      -227.2294633570     -327.8150156305        0.0024166993        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  7.962801270777765E-003
 Total number of micro iterations:    8

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99265319 pnorm= 0.0000E+00 rznorm= 3.2278E-06 rpnorm= 0.0000E+00 noldr=  8 nnewr=  8 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -41.245666  -22.598661   -2.910493   -1.723226   -1.350117   -1.328539
 i,qaaresolved                     1  -1.11354270838814     
 i,qaaresolved                     2 -0.705165880556201     
 i,qaaresolved                     3 -6.286898548410358E-003

 qvv(*) eigenvalues. symmetry block  1
     0.483307    0.653120    0.692112    1.342517    1.597155    1.694342    2.000922    2.290662    2.297362    2.381924
     2.440894    2.711070    3.736402

 fdd(*) eigenvalues. symmetry block  2
   -41.245673  -22.598685   -2.908694   -1.715056   -1.344466   -1.324469
 i,qaaresolved                     1  -1.07802775907265     
 i,qaaresolved                     2 -0.700498729181113     
 i,qaaresolved                     3  9.638426694682756E-002

 qvv(*) eigenvalues. symmetry block  2
     0.534338    0.703134    0.805139    1.596192    1.729924    1.967044    2.152410    2.346127    2.384350    2.470215
     2.571518    2.776555    3.856185

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    2 emc=   -227.4971900725 demc= 5.3784E-02 wnorm= 6.3702E-02 knorm= 1.2099E-01 apxde= 3.2618E-03    *not conv.*     

               starting mcscf iteration...   3

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     2, naopsy(1) =    22, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 917497823

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 917330879
 address segment size,           sizesg = 917179393
 number of in-core blocks,       nincbk =         4
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     129163 transformed 1/r12    array elements were written in      24 records.


 mosort: allocated sort2 space, avc2is=   917359676 available sort2 space, avcisx=   917359928

   4 trial vectors read from nvfile (unit= 29).
 ciiter=  14 noldhv=  7 noldv=  7

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.6783629936     -328.2639152672        0.0000008711        0.0000020584
    2      -227.5322596428     -328.1178119164        0.0000011404        0.0000020584
    3      -227.4156895197     -328.0012417933        0.0000015368        0.0000020584
    4      -227.3850544795     -327.9706067531        0.0024412671        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=   9 noldhv=  4 noldv=  4

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.5327018980     -328.1182541716        0.0000005353        0.0000020584
    2      -227.3444030521     -327.9299553257        0.0000015479        0.0000020584
    3      -227.2207646471     -327.8063169207        0.0012902111        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  1.194373560077489E-003
 Total number of micro iterations:    6

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-0.99996405 pnorm= 0.0000E+00 rznorm= 7.5351E-06 rpnorm= 0.0000E+00 noldr=  6 nnewr=  6 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -41.237027  -22.585187   -2.909135   -1.718115   -1.326702   -1.324426
 i,qaaresolved                     1  -1.11362824686424     
 i,qaaresolved                     2 -0.720081725975270     
 i,qaaresolved                     3  9.436761626764323E-004

 qvv(*) eigenvalues. symmetry block  1
     0.488335    0.657823    0.693209    1.337705    1.599178    1.700360    2.002965    2.296687    2.306296    2.378353
     2.440957    2.709263    3.734823

 fdd(*) eigenvalues. symmetry block  2
   -41.237034  -22.585210   -2.907339   -1.709987   -1.322539   -1.318226
 i,qaaresolved                     1  -1.07846891473025     
 i,qaaresolved                     2 -0.715969072245921     
 i,qaaresolved                     3  0.100260068560946     

 qvv(*) eigenvalues. symmetry block  2
     0.539522    0.708002    0.806777    1.595461    1.731945    1.970592    2.156191    2.354584    2.387216    2.474003
     2.565550    2.775462    3.854979

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    3 emc=   -227.5006834212 demc= 3.4933E-03 wnorm= 9.5550E-03 knorm= 8.4795E-03 apxde= 2.3679E-05    *not conv.*     

               starting mcscf iteration...   4

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     2, naopsy(1) =    22, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 917497823

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 917330879
 address segment size,           sizesg = 917179393
 number of in-core blocks,       nincbk =         4
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     129649 transformed 1/r12    array elements were written in      24 records.


 mosort: allocated sort2 space, avc2is=   917359676 available sort2 space, avcisx=   917359928

   4 trial vectors read from nvfile (unit= 29).
 ciiter=  17 noldhv= 12 noldv= 12

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.6783317067     -328.2638839803        0.0000002707        0.0000010000
    2      -227.5323881875     -328.1179404611        0.0000003791        0.0000010000
    3      -227.4156124383     -328.0011647119        0.0000003271        0.0000010000
    4      -227.3853460344     -327.9708983080        0.0025628966        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=  10 noldhv=  5 noldv=  5

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.5328247515     -328.1183770251        0.0000003559        0.0000010000
    2      -227.3443905919     -327.9299428655        0.0000004487        0.0000010000
    3      -227.2200178816     -327.8055701552        0.0010948384        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  2.229380538375214E-004
 Total number of micro iterations:    7

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0= 0.99999972 pnorm= 0.0000E+00 rznorm= 2.5853E-07 rpnorm= 0.0000E+00 noldr=  7 nnewr=  7 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -41.239069  -22.582782   -2.909547   -1.717662   -1.326883   -1.322278
 i,qaaresolved                     1  -1.11411655334845     
 i,qaaresolved                     2 -0.722011325126235     
 i,qaaresolved                     3  1.520539747110083E-003

 qvv(*) eigenvalues. symmetry block  1
     0.488813    0.658288    0.693326    1.337261    1.599380    1.700977    2.003158    2.297190    2.307180    2.378165
     2.440659    2.708602    3.734278

 fdd(*) eigenvalues. symmetry block  2
   -41.239076  -22.582805   -2.907753   -1.709552   -1.322707   -1.315984
 i,qaaresolved                     1  -1.07905536820899     
 i,qaaresolved                     2 -0.717990746398623     
 i,qaaresolved                     3  0.100029175335898     

 qvv(*) eigenvalues. symmetry block  2
     0.540024    0.708483    0.806932    1.595677    1.732134    1.970919    2.156584    2.355391    2.387261    2.474292
     2.565184    2.774909    3.854477

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    4 emc=   -227.5007095352 demc= 2.6114E-05 wnorm= 1.7835E-03 knorm= 7.4758E-04 apxde= 3.1590E-07    *not conv.*     

               starting mcscf iteration...   5

 orbital-state coupling will not be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     2, naopsy(1) =    22, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 917497823

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 917330879
 address segment size,           sizesg = 917179393
 number of in-core blocks,       nincbk =         4
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     124256 transformed 1/r12    array elements were written in      23 records.


 mosort: allocated sort2 space, avc2is=   917359676 available sort2 space, avcisx=   917359928

   4 trial vectors read from nvfile (unit= 29).
 ciiter=  17 noldhv= 12 noldv= 12

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.6783518491     -328.2639041227        0.0000005227        0.0000010000
    2      -227.5324037367     -328.1179560103        0.0000008179        0.0000010000
    3      -227.4155798541     -328.0011321277        0.0000003683        0.0000010000
    4      -227.3853570493     -327.9709093229        0.0025416563        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=  10 noldhv=  5 noldv=  5

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.5328397920     -328.1183920656        0.0000003444        0.0000010000
    2      -227.3443742637     -327.9299265373        0.0000004327        0.0000010000
    3      -227.2199684435     -327.8055207171        0.0011544535        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  4.327999073898266E-005
 Total number of micro iterations:    5

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-1.00000000 pnorm= 0.0000E+00 rznorm= 2.7016E-07 rpnorm= 0.0000E+00 noldr=  5 nnewr=  5 nolds=  0 nnews=  0
 

 fdd(*) eigenvalues. symmetry block  1
   -41.239564  -22.582443   -2.909649   -1.717620   -1.326944   -1.322073
 i,qaaresolved                     1  -1.11423485676550     
 i,qaaresolved                     2 -0.722258518541627     
 i,qaaresolved                     3  1.576783577161178E-003

 qvv(*) eigenvalues. symmetry block  1
     0.488869    0.658339    0.693337    1.337254    1.599406    1.701051    2.003168    2.297234    2.307280    2.378089
     2.440580    2.708494    3.734177

 fdd(*) eigenvalues. symmetry block  2
   -41.239571  -22.582467   -2.907854   -1.709513   -1.322767   -1.315771
 i,qaaresolved                     1  -1.07918993067357     
 i,qaaresolved                     2 -0.718244323702756     
 i,qaaresolved                     3  9.998716263448060E-002

 qvv(*) eigenvalues. symmetry block  2
     0.540083    0.708535    0.806947    1.595740    1.732158    1.970948    2.156633    2.355480    2.387239    2.474307
     2.565085    2.774815    3.854383

 restrt: restart information saved on the restart file (unit= 13).

 not all mcscf convergence criteria are satisfied.
 iter=    5 emc=   -227.5007098991 demc= 3.6393E-07 wnorm= 3.4624E-04 knorm= 8.4386E-05 apxde= 8.2385E-09    *not conv.*     

               starting mcscf iteration...   6

 orbital-state coupling will be calculated this iteration.

 *** Starting integral transformation ***

 module tranlib input parameters:

 prnopt    =     1, chkopt    =     0,ortopt    =     0, denopt    =     0
 mapin(1 ) =     1, nsymao    =     2, naopsy(1) =    22, freeze(1) =     1
 mapout(1) =     1, nsymmo    =    -1, nmopsy(1) =    -1, fsplit    =     1
 outlab    =     0, seward    =     0, lumorb    =     0, DALTON2   =     0
 nextint   =     2
 LDAMIN    =   127, LDAMAX    = 64959, LDAINC    =    64
 LRC1MX    =    -1, LRC2MX    =    -1, LRCSCR    = 65000

 THRESH    =  5.0000E-12  [cutoff threshold]

 module tranlib: workspace lcore= 917497823

 inoutp: segmentation information:
 in-core transformation space,   avcinc = 917330879
 address segment size,           sizesg = 917179393
 number of in-core blocks,       nincbk =         4
 number of out-of-core blocks,   noutbk =         0
 number of in-core segments,     incseg =         1
 number of out-of-core segments, outseg =         0
 trmain:     121959 transformed 1/r12    array elements were written in      23 records.


 mosort: allocated sort2 space, avc2is=   917359676 available sort2 space, avcisx=   917359928

   4 trial vectors read from nvfile (unit= 29).
 ciiter=  17 noldhv= 12 noldv= 12

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.6783581879     -328.2639104615        0.0000005155        0.0000010000
    2      -227.5324063198     -328.1179585933        0.0000009356        0.0000010000
    3      -227.4155723468     -328.0011246204        0.0000003608        0.0000010000
    4      -227.3853557946     -327.9709080682        0.0025169916        0.0100000000

   3 trial vectors read from nvfile (unit= 29).
 ciiter=  10 noldhv=  5 noldv=  5

 Eigenvalues of the hmc(*) matrix
             total energy     electronic energy      residual norm          rtolci(*)
    1*     -227.5328423294     -328.1183946030        0.0000003404        0.0000010000
    2      -227.3443703603     -327.9299226339        0.0000004301        0.0000010000
    3      -227.2199654270     -327.8055177006        0.0011626873        0.0100000000
 
  tol(10)=  0.000000000000000E+000  eshsci=  8.385926128170196E-006
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 performing all-state projection
 Total number of micro iterations:    5

 ***  micro: final psci convergence values:  ***
    imxov=  1 z0=-1.00000000 pnorm= 2.4935E-04 rznorm= 3.9996E-07 rpnorm= 8.1973E-07 noldr=  5 nnewr=  5 nolds=  4 nnews=  4
 

 fdd(*) eigenvalues. symmetry block  1
   -41.239663  -22.582388   -2.909670   -1.717615   -1.326958   -1.322048
 i,qaaresolved                     1  -1.11425943070208     
 i,qaaresolved                     2 -0.722296436358318     
 i,qaaresolved                     3  1.581769308999670E-003

 qvv(*) eigenvalues. symmetry block  1
     0.488877    0.658347    0.693339    1.337259    1.599410    1.701062    2.003168    2.297238    2.307295    2.378072
     2.440563    2.708475    3.734158

 fdd(*) eigenvalues. symmetry block  2
   -41.239670  -22.582411   -2.907875   -1.709509   -1.322780   -1.315746
 i,qaaresolved                     1  -1.07921734674178     
 i,qaaresolved                     2 -0.718282479783208     
 i,qaaresolved                     3  9.998305680110708E-002

 qvv(*) eigenvalues. symmetry block  2
     0.540091    0.708543    0.806949    1.595751    1.732162    1.970952    2.156641    2.355493    2.387233    2.474307
     2.565063    2.774799    3.854365

 restrt: restart information saved on the restart file (unit= 13).

 all mcscf convergence criteria are satisfied.

 final mcscf convergence values:
 iter=    6 emc=   -227.5007099089 demc= 9.7374E-09 wnorm= 6.7087E-05 knorm= 1.5020E-05 apxde= 3.4669E-10    *converged*     




   ---------Individual total energies for all states:----------
   DRT #1 state # 1 wt 0.200 total energy=     -227.678358188, rel. (eV)=   0.000000
   DRT #1 state # 2 wt 0.200 total energy=     -227.532406320, rel. (eV)=   3.971554
   DRT #1 state # 3 wt 0.200 total energy=     -227.415572347, rel. (eV)=   7.150770
   DRT #2 state # 1 wt 0.200 total energy=     -227.532842329, rel. (eV)=   3.959690
   DRT #2 state # 2 wt 0.200 total energy=     -227.344370360, rel. (eV)=   9.088275
   ------------------------------------------------------------


 MO-coefficient print-out skipped (no flag 32)
 They may be found in the MOCOEF directory.

          natural orbitals of the final iteration,block  1    -  A' 
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.84083382     1.49985181
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.75038980     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          natural orbitals of the final iteration,block  2    -  A" 
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.82524697     1.49938167
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.58429593     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
 d1(*), fmc(*), and qmc(*) written to the 1-particle density matrix file.
        555 d2(*) elements written to the 2-particle density matrix file: mcd2fl                                                      
 Computing the requested mcscf (transition) density matrices (flag 30)
 Reading mcdenin ...
 Number of density matrices (ndens):                     6
 Number of unique bra states (ndbra):                     3
 qind: F
 (Transition) density matrices:
 d1(*) written to the 1-particle density matrix file.
        555 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st01                                           
 d1(*) written to the 1-particle density matrix file.
 d1(*) written to the 1-particle density matrix file.
        555 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st02-st01                                      
 d1(*) written to the 1-particle density matrix file.
 d1(*) written to the 1-particle density matrix file.
        555 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st03-st01                                      
 d1(*) written to the 1-particle density matrix file.
        555 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st02                                           
 d1(*) written to the 1-particle density matrix file.
 d1(*) written to the 1-particle density matrix file.
        555 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st03-st02                                      
 d1(*) written to the 1-particle density matrix file.
        555 d2(*) elements written to the 2-particle density matrix file: mcsd2fl.drt1.st03                                           

          state spec. NOs: DRT 1, State  1
 *** warning *** large active-orbital occupation. i=  3 nocc= 1.9990E+00
 *** warning *** large active-orbital occupation. i=  3 nocc= 1.9990E+00

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99903045     1.90153558
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.10580613     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          block  2
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99902510     1.89477416
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.09982857     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          NOs of symm. trans. dens. mat.: DRT 1, State  2 - DRT  1, state   1
 *** warning *** small active-orbital occupation. i=  1 nocc=-4.6982E-01
 *** warning *** small active-orbital occupation. i=  2 nocc= 1.4249E-07
 *** warning *** small active-orbital occupation. i=  1 nocc=-4.2947E-01
 *** warning *** small active-orbital occupation. i=  2 nocc= 1.4583E-07

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.46981510     0.00000014
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=    -0.46981524     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          block  2
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.42947437     0.00000015
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=    -0.42947453     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          NOs of symm. trans. dens. mat.: DRT 1, State  3 - DRT  1, state   1
 *** warning *** small active-orbital occupation. i=  1 nocc=-2.9101E-03
 *** warning *** small active-orbital occupation. i=  2 nocc=-1.7598E-04
 *** warning *** small active-orbital occupation. i=  3 nocc= 9.3596E-05
 *** warning *** small active-orbital occupation. i=  1 nocc=-8.7911E-05
 *** warning *** small active-orbital occupation. i=  2 nocc= 1.1958E-04

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00009360    -0.00017598
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=    -0.00291010     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          block  2
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00296081     0.00011958
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=    -0.00008791     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          state spec. NOs: DRT 1, State  2

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.91653460     1.46177415
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.63522503     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          block  2
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.90888804     1.53726683
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.54031134     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          NOs of symm. trans. dens. mat.: DRT 1, State  3 - DRT  1, state   2
 *** warning *** small active-orbital occupation. i=  1 nocc=-5.8783E-02
 *** warning *** small active-orbital occupation. i=  2 nocc= 5.9266E-09
 *** warning *** small active-orbital occupation. i=  1 nocc=-2.2921E-02
 *** warning *** small active-orbital occupation. i=  2 nocc= 2.8377E-08

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.05878290     0.00000001
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=    -0.05878291     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          block  2
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.02292079     0.00000003
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=    -0.02292081     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          state spec. NOs: DRT 1, State  3
 *** warning *** large active-orbital occupation. i=  3 nocc= 1.9999E+00
 *** warning *** large active-orbital occupation. i=  3 nocc= 1.9999E+00

          block  1
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99991532     1.16822491
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     1.00636842     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000

          block  2
               MO    1        MO    2        MO    3        MO    4        MO    5        MO    6        MO    7        MO    8
  occ(*)=     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     2.00000000     1.99985754     0.99363263
               MO    9        MO   10        MO   11        MO   12        MO   13        MO   14        MO   15        MO   16
  occ(*)=     0.83200119     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000
               MO   17        MO   18        MO   19        MO   20        MO   21        MO   22
  occ(*)=     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A'  partial gross atomic populations
   ao class       1A'        2A'        3A'        4A'        5A'        6A' 
     3_ s       0.000630   1.999828   0.268234   0.882840   0.053902   0.000000
     3_ p       0.000713   0.000007   0.129436   0.151074   0.463321   1.089076
     4_ s       1.998648  -0.000185   1.459864   0.262734   0.215317   0.000000
     4_ p       0.000005  -0.000578   0.138928   0.208325   1.042994   0.129044
     4_ s       0.000002   0.000464   0.001769   0.247513   0.112234   0.390940
     4_ s       0.000002   0.000464   0.001769   0.247513   0.112234   0.390940
 
   ao class       7A'        8A'        9A'       10A'       11A'       12A' 
     3_ s       0.000237   0.000000  -0.000049   0.000000   0.000000   0.000000
     3_ p       0.609346   0.007090   0.511081   0.000000   0.000000   0.000000
     4_ s      -0.000060  -0.000000   0.000005   0.000000   0.000000   0.000000
     4_ p       1.231220   1.333255   0.239113   0.000000   0.000000   0.000000
     4_ s       0.000045   0.079754   0.000120   0.000000   0.000000   0.000000
     4_ s       0.000045   0.079754   0.000120   0.000000   0.000000   0.000000
 
   ao class      13A'       14A'       15A'       16A'       17A'       18A' 
 
   ao class      19A'       20A'       21A'       22A' 

                        A"  partial gross atomic populations
   ao class       1A"        2A"        3A"        4A"        5A"        6A" 
     3_ s       0.000725   1.999309   0.258160   0.889373   0.053816   0.000000
     3_ p       0.000767   0.000006   0.125521   0.143110   0.471415   1.097909
     4_ s       1.998505   0.000013   1.471292   0.266839   0.214544   0.000000
     4_ p       0.000005  -0.000399   0.140983   0.213701   1.032579   0.129341
     4_ s      -0.000001   0.000535   0.002022   0.243489   0.113823   0.386375
     4_ s      -0.000001   0.000535   0.002022   0.243489   0.113823   0.386375
 
   ao class       7A"        8A"        9A"       10A"       11A"       12A" 
     3_ s       0.000291  -0.000000  -0.000058   0.000000   0.000000   0.000000
     3_ p       0.511085   0.007627   0.425617   0.000000   0.000000   0.000000
     4_ s       0.000022   0.000000   0.000005   0.000000   0.000000   0.000000
     4_ p       1.313757   1.335075   0.158565   0.000000   0.000000   0.000000
     4_ s       0.000046   0.078340   0.000083   0.000000   0.000000   0.000000
     4_ s       0.000046   0.078340   0.000083   0.000000   0.000000   0.000000
 
   ao class      13A"       14A"       15A"       16A"       17A"       18A" 
 
   ao class      19A"       20A"       21A"       22A" 


                        gross atomic populations
     ao            3_         4_         4_         4_
      s         6.407239   7.887543   1.657553   1.657553
      p         5.744199   8.645913   0.000000   0.000000
    total      12.151438  16.533457   1.657553   1.657553
 

 Total number of electrons:   32.00000000

 Mulliken population for:
DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A'  partial gross atomic populations
   ao class       1A'        2A'        3A'        4A'        5A'        6A' 
     3_ s       0.000630   1.999829   0.268246   0.882842   0.053888   0.000000
     3_ p       0.000713   0.000007   0.129436   0.151100   0.463295   1.089076
     4_ s       1.998648  -0.000185   1.459845   0.262726   0.215344   0.000000
     4_ p       0.000005  -0.000578   0.138935   0.208274   1.043039   0.129044
     4_ s       0.000002   0.000464   0.001769   0.247529   0.112217   0.390940
     4_ s       0.000002   0.000464   0.001769   0.247529   0.112217   0.390940
 
   ao class       7A'        8A'        9A'       10A'       11A'       12A' 
     3_ s       0.000000   0.000265  -0.000008   0.000000   0.000000   0.000000
     3_ p       0.009449   0.850334   0.059772   0.000000   0.000000   0.000000
     4_ s      -0.000000  -0.000057   0.000000   0.000000   0.000000   0.000000
     4_ p       1.776987   1.050843   0.046011   0.000000   0.000000   0.000000
     4_ s       0.106297   0.000075   0.000015   0.000000   0.000000   0.000000
     4_ s       0.106297   0.000075   0.000015   0.000000   0.000000   0.000000
 
   ao class      13A'       14A'       15A'       16A'       17A'       18A' 
 
   ao class      19A'       20A'       21A'       22A' 

                        A"  partial gross atomic populations
   ao class       1A"        2A"        3A"        4A"        5A"        6A" 
     3_ s       0.000725   1.999309   0.258172   0.889376   0.053802   0.000000
     3_ p       0.000767   0.000006   0.125521   0.143137   0.471388   1.097909
     4_ s       1.998505   0.000013   1.471273   0.266831   0.214572   0.000000
     4_ p       0.000005  -0.000399   0.140990   0.213648   1.032626   0.129341
     4_ s      -0.000001   0.000535   0.002023   0.243505   0.113806   0.386375
     4_ s      -0.000001   0.000535   0.002023   0.243505   0.113806   0.386375
 
   ao class       7A"        8A"        9A"       10A"       11A"       12A" 
     3_ s       0.000000   0.000351  -0.000012   0.000000   0.000000   0.000000
     3_ p       0.010169   0.803970   0.058312   0.000000   0.000000   0.000000
     4_ s       0.000000   0.000016   0.000001   0.000000   0.000000   0.000000
     4_ p       1.779966   1.090279   0.041502   0.000000   0.000000   0.000000
     4_ s       0.104445   0.000079   0.000013   0.000000   0.000000   0.000000
     4_ s       0.104445   0.000079   0.000013   0.000000   0.000000   0.000000
 
   ao class      13A"       14A"       15A"       16A"       17A"       18A" 
 
   ao class      19A"       20A"       21A"       22A" 


                        gross atomic populations
     ao            3_         4_         4_         4_
      s         6.407412   7.887532   1.710089   1.710089
      p         5.464361   8.820516   0.000000   0.000000
    total      11.871774  16.708049   1.710089   1.710089
 

 Total number of electrons:   32.00000000

 Off-diagonal Mulliken population for:
DRT 1,state 02 - DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A'  partial gross atomic populations
   ao class       1A'        2A'        3A'        4A'        5A'        6A' 
 
   ao class       7A'        8A'        9A'       10A'       11A'       12A' 
     3_ s      -0.000019   0.000000   0.000019   0.000000   0.000000   0.000000
     3_ p       0.112146   0.000000  -0.112146   0.000000   0.000000   0.000000
     4_ s       0.000000  -0.000000  -0.000000   0.000000   0.000000   0.000000
     4_ p       0.332644   0.000000  -0.332644   0.000000   0.000000   0.000000
     4_ s       0.012502   0.000000  -0.012542   0.000000   0.000000   0.000000
     4_ s       0.012542   0.000000  -0.012502   0.000000   0.000000   0.000000
 
   ao class      13A'       14A'       15A'       16A'       17A'       18A' 
 
   ao class      19A'       20A'       21A'       22A' 

                        A"  partial gross atomic populations
   ao class       1A"        2A"        3A"        4A"        5A"        6A" 
 
   ao class       7A"        8A"        9A"       10A"       11A"       12A" 
     3_ s      -0.000029   0.000000   0.000029   0.000000   0.000000   0.000000
     3_ p       0.103195   0.000000  -0.103195   0.000000   0.000000   0.000000
     4_ s       0.000003   0.000000  -0.000003   0.000000   0.000000   0.000000
     4_ p       0.303818   0.000000  -0.303818   0.000000   0.000000   0.000000
     4_ s       0.011274   0.000000  -0.011214   0.000000   0.000000   0.000000
     4_ s       0.011214   0.000000  -0.011274   0.000000   0.000000   0.000000
 
   ao class      13A"       14A"       15A"       16A"       17A"       18A" 
 
   ao class      19A"       20A"       21A"       22A" 


                        gross atomic populations
     ao            3_         4_         4_         4_
      s         0.000000  -0.000000   0.000019  -0.000019
      p         0.000000  -0.000000   0.000000   0.000000
    total       0.000000  -0.000000   0.000019  -0.000019
 

 Total number of electrons:    0.00000000

 Off-diagonal Mulliken population for:
DRT 1,state 03 - DRT 1, state 01


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A'  partial gross atomic populations
   ao class       1A'        2A'        3A'        4A'        5A'        6A' 
 
   ao class       7A'        8A'        9A'       10A'       11A'       12A' 
     3_ s       0.000000  -0.000000   0.000000   0.000000   0.000000   0.000000
     3_ p      -0.000002  -0.000183  -0.000014   0.000000   0.000000   0.000000
     4_ s      -0.000000  -0.000000   0.000000   0.000000   0.000000   0.000000
     4_ p       0.000096   0.000007  -0.002587   0.000000   0.000000   0.000000
     4_ s       0.000000  -0.000000  -0.000155   0.000000   0.000000   0.000000
     4_ s       0.000000  -0.000000  -0.000155   0.000000   0.000000   0.000000
 
   ao class      13A'       14A'       15A'       16A'       17A'       18A' 
 
   ao class      19A'       20A'       21A'       22A' 

                        A"  partial gross atomic populations
   ao class       1A"        2A"        3A"        4A"        5A"        6A" 
 
   ao class       7A"        8A"        9A"       10A"       11A"       12A" 
     3_ s      -0.000000   0.000000  -0.000000   0.000000   0.000000   0.000000
     3_ p       0.000015   0.000121   0.000000   0.000000   0.000000   0.000000
     4_ s       0.000000  -0.000000  -0.000000   0.000000   0.000000   0.000000
     4_ p       0.002636  -0.000001  -0.000088   0.000000   0.000000   0.000000
     4_ s       0.000155   0.000000  -0.000000   0.000000   0.000000   0.000000
     4_ s       0.000155   0.000000  -0.000000   0.000000   0.000000   0.000000
 
   ao class      13A"       14A"       15A"       16A"       17A"       18A" 
 
   ao class      19A"       20A"       21A"       22A" 


                        gross atomic populations
     ao            3_         4_         4_         4_
      s        -0.000000  -0.000000  -0.000000  -0.000000
      p        -0.000062   0.000063   0.000000   0.000000
    total      -0.000062   0.000063  -0.000000  -0.000000
 

 Total number of electrons:    0.00000000

 Mulliken population for:
DRT 1, state 02


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A'  partial gross atomic populations
   ao class       1A'        2A'        3A'        4A'        5A'        6A' 
     3_ s       0.000630   1.999829   0.268246   0.882842   0.053888   0.000000
     3_ p       0.000713   0.000007   0.129436   0.151100   0.463295   1.089076
     4_ s       1.998648  -0.000185   1.459845   0.262726   0.215344   0.000000
     4_ p       0.000005  -0.000578   0.138935   0.208274   1.043039   0.129044
     4_ s       0.000002   0.000464   0.001769   0.247529   0.112217   0.390940
     4_ s       0.000002   0.000464   0.001769   0.247529   0.112217   0.390940
 
   ao class       7A'        8A'        9A'       10A'       11A'       12A' 
     3_ s       0.000243  -0.000000  -0.000040   0.000000   0.000000   0.000000
     3_ p       0.604936   0.006910   0.442410   0.000000   0.000000   0.000000
     4_ s      -0.000063   0.000000   0.000004   0.000000   0.000000   0.000000
     4_ p       1.311331   1.299406   0.192645   0.000000   0.000000   0.000000
     4_ s       0.000044   0.077729   0.000102   0.000000   0.000000   0.000000
     4_ s       0.000044   0.077729   0.000102   0.000000   0.000000   0.000000
 
   ao class      13A'       14A'       15A'       16A'       17A'       18A' 
 
   ao class      19A'       20A'       21A'       22A' 

                        A"  partial gross atomic populations
   ao class       1A"        2A"        3A"        4A"        5A"        6A" 
     3_ s       0.000725   1.999309   0.258172   0.889376   0.053802   0.000000
     3_ p       0.000767   0.000006   0.125521   0.143137   0.471388   1.097909
     4_ s       1.998505   0.000013   1.471273   0.266831   0.214572   0.000000
     4_ p       0.000005  -0.000399   0.140990   0.213648   1.032626   0.129341
     4_ s      -0.000001   0.000535   0.002023   0.243505   0.113806   0.386375
     4_ s      -0.000001   0.000535   0.002023   0.243505   0.113806   0.386375
 
   ao class       7A"        8A"        9A"       10A"       11A"       12A" 
     3_ s       0.000313  -0.000000  -0.000056   0.000000   0.000000   0.000000
     3_ p       0.574655   0.007820   0.382213   0.000000   0.000000   0.000000
     4_ s       0.000022   0.000000   0.000005   0.000000   0.000000   0.000000
     4_ p       1.333793   1.368809   0.157998   0.000000   0.000000   0.000000
     4_ s       0.000053   0.080319   0.000076   0.000000   0.000000   0.000000
     4_ s       0.000053   0.080319   0.000076   0.000000   0.000000   0.000000
 
   ao class      13A"       14A"       15A"       16A"       17A"       18A" 
 
   ao class      19A"       20A"       21A"       22A" 


                        gross atomic populations
     ao            3_         4_         4_         4_
      s         6.407277   7.887540   1.657487   1.657487
      p         5.691298   8.698910   0.000000   0.000000
    total      12.098576  16.586450   1.657487   1.657487
 

 Total number of electrons:   32.00000000

 Off-diagonal Mulliken population for:
DRT 1,state 03 - DRT 1, state 02


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A'  partial gross atomic populations
   ao class       1A'        2A'        3A'        4A'        5A'        6A' 
 
   ao class       7A'        8A'        9A'       10A'       11A'       12A' 
     3_ s      -0.000002   0.000000   0.000002   0.000000   0.000000   0.000000
     3_ p       0.021272   0.000000  -0.021272   0.000000   0.000000   0.000000
     4_ s       0.000000  -0.000000  -0.000000   0.000000   0.000000   0.000000
     4_ p       0.034377   0.000000  -0.034377   0.000000   0.000000   0.000000
     4_ s       0.001565   0.000000  -0.001570   0.000000   0.000000   0.000000
     4_ s       0.001570   0.000000  -0.001565   0.000000   0.000000   0.000000
 
   ao class      13A'       14A'       15A'       16A'       17A'       18A' 
 
   ao class      19A'       20A'       21A'       22A' 

                        A"  partial gross atomic populations
   ao class       1A"        2A"        3A"        4A"        5A"        6A" 
 
   ao class       7A"        8A"        9A"       10A"       11A"       12A" 
     3_ s       0.000001  -0.000000  -0.000001   0.000000   0.000000   0.000000
     3_ p       0.012033  -0.000000  -0.012033   0.000000   0.000000   0.000000
     4_ s      -0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
     4_ p       0.009686   0.000000  -0.009686   0.000000   0.000000   0.000000
     4_ s       0.000591   0.000000  -0.000610   0.000000   0.000000   0.000000
     4_ s       0.000610   0.000000  -0.000591   0.000000   0.000000   0.000000
 
   ao class      13A"       14A"       15A"       16A"       17A"       18A" 
 
   ao class      19A"       20A"       21A"       22A" 


                        gross atomic populations
     ao            3_         4_         4_         4_
      s         0.000000   0.000000  -0.000024   0.000024
      p         0.000000  -0.000000   0.000000   0.000000
    total       0.000000  -0.000000  -0.000024   0.000024
 

 Total number of electrons:    0.00000000

 Mulliken population for:
DRT 1, state 03


          Mulliken population analysis


  NOTE: For HERMIT use spherical harmonics basis sets !!!
 

                        A'  partial gross atomic populations
   ao class       1A'        2A'        3A'        4A'        5A'        6A' 
     3_ s       0.000630   1.999829   0.268246   0.882842   0.053888   0.000000
     3_ p       0.000713   0.000007   0.129436   0.151100   0.463295   1.089076
     4_ s       1.998648  -0.000185   1.459845   0.262726   0.215344   0.000000
     4_ p       0.000005  -0.000578   0.138935   0.208274   1.043039   0.129044
     4_ s       0.000002   0.000464   0.001769   0.247529   0.112217   0.390940
     4_ s       0.000002   0.000464   0.001769   0.247529   0.112217   0.390940
 
   ao class       7A'        8A'        9A'       10A'       11A'       12A' 
     3_ s       0.000200  -0.000042  -0.000000   0.000000   0.000000   0.000000
     3_ p       0.308203   1.002331   0.004757   0.000000   0.000000   0.000000
     4_ s      -0.000069   0.000010   0.000000   0.000000   0.000000   0.000000
     4_ p       1.691554   0.165511   0.894585   0.000000   0.000000   0.000000
     4_ s       0.000014   0.000207   0.053513   0.000000   0.000000   0.000000
     4_ s       0.000014   0.000207   0.053513   0.000000   0.000000   0.000000
 
   ao class      13A'       14A'       15A'       16A'       17A'       18A' 
 
   ao class      19A'       20A'       21A'       22A' 

                        A"  partial gross atomic populations
   ao class       1A"        2A"        3A"        4A"        5A"        6A" 
     3_ s       0.000725   1.999309   0.258172   0.889376   0.053802   0.000000
     3_ p       0.000767   0.000006   0.125521   0.143137   0.471388   1.097909
     4_ s       1.998505   0.000013   1.471273   0.266831   0.214572   0.000000
     4_ p       0.000005  -0.000399   0.140990   0.213648   1.032626   0.129341
     4_ s      -0.000001   0.000535   0.002023   0.243505   0.113806   0.386375
     4_ s      -0.000001   0.000535   0.002023   0.243505   0.113806   0.386375
 
   ao class       7A"        8A"        9A"       10A"       11A"       12A" 
     3_ s       0.000248  -0.000000  -0.000053   0.000000   0.000000   0.000000
     3_ p       0.306948   0.005054   0.711320   0.000000   0.000000   0.000000
     4_ s       0.000031   0.000000   0.000004   0.000000   0.000000   0.000000
     4_ p       1.692576   0.884748   0.120473   0.000000   0.000000   0.000000
     4_ s       0.000027   0.051915   0.000129   0.000000   0.000000   0.000000
     4_ s       0.000027   0.051915   0.000129   0.000000   0.000000   0.000000
 
   ao class      13A"       14A"       15A"       16A"       17A"       18A" 
 
   ao class      19A"       20A"       21A"       22A" 


                        gross atomic populations
     ao            3_         4_         4_         4_
      s         6.407170   7.887547   1.604969   1.604969
      p         6.010968   8.484376   0.000000   0.000000
    total      12.418138  16.371923   1.604969   1.604969
 

 Total number of electrons:   32.00000000

 !timer: mcscf                           cpu_time=     0.182 walltime=     0.201
