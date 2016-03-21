============================
file: LakeSuperiorN.lim

Solve the model in R with:

underdetermined:

require(LIM)
web.lim <- Setup("GreatLakesN3.lim.R")
pars<-Ldei(web.lim)
webranges<-Xranges(web.lim, ispos=TRUE, tol = 1e-3)
data.frame(webranges,parsimonious=pars$X)

xlim<-range(webranges)
dotchart(x=pars$X,labels=rownames(webranges),xlim=xlim,main="N-flux",pch=16)
cc<-1:nrow(webranges)
segments(x0=webranges[,1],y0=cc,x1=webranges[,2],y1=cc)

xs<-Xsample(web.lim)
pairs(xs,upper.panel = NULL, diag.panel = panel.hist, pch=".", cex=1, main="practice figures")

overdetermined:

web.lim <- Setup("GreatLakesN3.lim.R")
pars<-Lsei(web.lim)
data.frame(parsimonious=pars$X)

=============================

## EXTERNAL
	zSup_NH4
	zSup_NO3	
	zSup_orgN
	zHur_NH4
	zHur_NO3
	zHur_orgN
	zErie_NH4
	zErie_NO3
	zErie_orgN
	ySup_sedN
	yHur_sedN
	yErie_NH4
	yErie_NO3
	yErie_orgN
	yErie_sedN
## END EXTERNAL

## COMPONENT
	Sup_NH4	
	Sup_NO3	
	Sup_orgN	
	Sup_sedN
	Hur_NH4	
	Hur_NO3	
	Hur_orgN	
	Hur_sedN
	Erie_NH4	
	Erie_NO3	
	Erie_orgN	
	Erie_sedN
## END COMPONENT

## FLOWS  
	a1 : zSup_NH4->Sup_NH4
	a2 : zSup_NO3->Sup_NO3
	a3 : zSup_orgN->Sup_orgN
	a4 : zHur_NH4->Hur_NH4
	a5 : zHur_NO3->Hur_NO3
	a6 : zHur_orgN->Hur_orgN
	a7 : zErie_NH4->Erie_NH4
	a8 : zErie_NO3->Erie_NO3
	a9 : zErie_orgN->Erie_orgN
	a10 : Sup_NH4->Sup_NO3
	a11 : Sup_NH4->Sup_orgN
	a12 : Sup_NO3->Sup_orgN
	a13 : Sup_orgN->Sup_NH4
	a14 : Sup_orgN->Sup_sedN
	a15 : Sup_sedN->Sup_NH4
	a16 : Sup_sedN->Sup_NO3
	a17 : Sup_sedN->Sup_orgN
	a18 : Sup_sedN->ySup_sedN
	a19 : Sup_NH4->Hur_NH4
	a20 : Sup_NO3->Hur_NO3
	a21 : Sup_orgN->Hur_orgN
	a22 : Hur_NH4->Hur_NO3
	a23 : Hur_NH4->Hur_orgN
	a24 : Hur_NO3->Hur_orgN
	a25 : Hur_NO3->Hur_sedN
	a26 : Hur_orgN->Hur_NH4
	a27 : Hur_orgN->Hur_sedN
	a28 : Hur_sedN->Hur_NH4
	a29 : Hur_sedN->yHur_sedN
	a30 : Hur_NH4->Erie_NH4
	a31 : Hur_NO3->Erie_NO3
	a32 : Hur_orgN->Erie_orgN
	a33 : Erie_NH4->Erie_NO3
	a34 : Erie_NH4->Erie_orgN
	a35 : Erie_NO3->Erie_orgN
	a36 : Erie_NO3->Erie_sedN
	a37 : Erie_orgN->Erie_NH4
	a38 : Erie_orgN->Erie_sedN
	a39 : Erie_sedN->Erie_NH4
	a40 : Erie_sedN->Erie_orgN
	a41 : Erie_NH4->yErie_NH4
	a42 : Erie_NO3->yErie_NO3
	a43 : Erie_orgN->yErie_orgN
	a44 : Erie_sedN->yErie_sedN		
## END FLOWS

## PARAMETERS

## END PARAMETERS

## VARIABLES
	
## END VARIABLES

## EQUALITIES
	a1 + a13 + a15 - a10 - a11 - a19 = 0             !Sup_NH4
	a2 + a10 + a16 - a12 - a20 = 0                   !Sup_NO3
	a3 + a11 + a12 + a17 - a13 - a14 - a21 = 0       !Sup_orgN
	a14 - a15 - a16 - a17 - a18 = 0                  !Sup_sedN
	a4 + a19 + a26 + a28 - a22 - a23 - a30 = 0       !Hur_NH4
	a5 + a20 + a22 - a24 - a25 - a31 = 0             !Hur_NO3
	a6 + a21 + a23 + a24 - a26 - a27 - a32 = 0       !Hur_orgN
	a25 + a27 - a28 - a29 = 0                        !Hur_sedN
	a7 + a30 + a37 + a39 - a33 - a34 - a41 = 0       !Erie_NH4
	a8 + a31 + a33 - a35 - a36 - a42 = 0             !Erie_NO3
	a9 + a32 + a34 + a35 + a40 - a37 - a38 - a43 = 0 !Erie_orgN
	a36 + a38 - a39 - a40 - a44 = 0                  !Erie_sedN
## END EQUALITIES

## INEQUALITIES 
!Note all exponents are decreased by 10^3.  Does not like constraints >10^8
	a1 = [1.49e4, 2.23e4]   ! zSup_NH4->Sup_NH4
	a2 = [2.22e4, 3.33e4]   ! zSup_NO3->Sup_NO3
	a3 = [2.08e4, 3.12e4]   ! zSup_orgN->Sup_orgN
	a4 = [5.30e3, 2.12e4]   ! zHur_NH4->Hur_NH4
	a5 = [3.38e4, 1.35e5]   ! zHur_NO3->Hur_NO3
	a6 = [1.05e4, 4.21e4]   ! zHur_orgN->Hur_orgN
	a7 = [4.87e3, 1.95e4]   ! zErie_NH4->Erie_NH4
	a8 = [4.68e4, 1.87e5]   ! zErie_NO3->Erie_NO3
	a9 = [2.81e4, 1.12e5]   ! zErie_orgN->Erie_orgN
	a10 = [1.10e6, 1.66e6]  ! Sup_NH4->Sup_NO3
	a11 = [2.91e5,6.80e5]  ! Sup_NH4->Sup_orgN
	a12 = [8.31e5, 1.94e6] ! Sup_NO3->Sup_orgN
	!a13 >1e5	           ! Sup_orgN->Sup_NH4
	!a14 = [4.14e4, 9.66e4] ! Sup_orgN->Sup_sedN
	a15 = [4.20e1, 4.20e3]  ! Sup_sedN->Sup_NH4
	a16 = [5.20e3, 2.08e4] ! Sup_sedN->Sup_NO3 
	a17 = [1.49e4, 3.47e4] ! Sup_sedN->Sup_orgN
	a18 = [7.48e4, 2.99e5] ! Sup_sedN->ySup_sedN
	a19 = [1.06e2, 4.22e2]  ! Sup_NH4->Hur_NH4
	a20 = [9.66e3, 3.86e4]  ! Sup_NO3->Hur_NO3
	a21 = [1.84e3, 7.36e3]  ! Sup_orgN->Hur_orgN
	a22 = [1.33e5, 5.31e5]  ! Hur_NH4->Hur_NO3
	a23 = [7.94e5, 3.18e6]  ! Hur_NH4->Hur_orgN
	a24 = [8.63e4,3.45e5]  ! Hur_NO3->Hur_orgN
	a25 = [4.64e3,1.86e4]  ! Hur_NO3->Hur_sedN
	!a26 = [1e2,1e12]       ! Hur_orgN->Hur_NH4
	!a27 = [1.57e2,3.67e6]  ! Hur_orgN->Hur_sedN
	a28 = [2.68e3,1.07e4]  ! Hur_sedN->Hur_NH4
	a29 = [5.02e4, 2.37e5]  ! Hur_sedN->yHur_sedN
	a30 = [1.52e3, 6.10e3]  ! Hur_NH4->Erie_NH4
	a31 = [2.18e4, 8.70e4]  ! Hur_NO3->Erie_NO3
	a32 = [8.60e3, 3.44e4]  ! Hur_orgN->Erie_orgN
	a33 = [3.85e5,1.54e6]  ! Erie_NH4->Erie_NO3
	a34 = [2.83e5,1.13e6]  ! Erie_NH4->Erie_orgN
	!a35 = [3.53e4,8.25e4]  ! Erie_NO3->Erie_orgN
	a36 = [9.24e3,3.70e4]  ! Erie_NO3->Erie_sedN
	!a37 = [1e2,1e8]       ! Erie_orgN->Erie_NH4
	!a38 = [1.21e2,2.82e7]  ! Erie_orgN->Erie_sedN
	a39 = [6.36e3, 2.54e4]  ! Erie_sedN->Erie_NH4
	a40 = [1.97e4, 7.84e4]  ! Erie_sedN->Erie_orgN
	a41 = [1.94e3, 7.76e3]  ! Erie_NH4->yErie_NH4
	a42 = [9.36e3, 3.74e4]  ! Erie_NO3->yErie_NO3
	a43 = [1.27e4, 5.07e4]  ! Erie_orgN->yErie_orgN
	a44 = [1.96e4, 7.84e4]  ! Erie_sedN->yErie_sedN
## END INEQUALITIES

