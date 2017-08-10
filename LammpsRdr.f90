program	LammpsMdRdr
	character dumString*198,logFile*88,outFile*88,firstWord*4,CallIntraVdw*3,dumString2*20
	character(len=15) timeStepStr,equil_percentStr
        integer date_time(8)
	character*10 d(3)
	double precision :: AvgIntraVdw

     	call date_and_time(d(1), d(2), d(3), date_time)

	call getarg(1,logFile)
	call getarg(2,outFile)
	call getarg(3,equil_percentStr)

    	read (equil_percentStr, *) equil_percent

	call IntraVdwCalc(AvgIntraVdw,nMolec)

	open(61,file=outFile)
	write(61,*)TRIM(outFile),'	Date=',d(1)
	write(6,'(A9,1x,A9,1x,3(A9,1x),6(A15,1x),A6)')'T(K)','rho(gcc)','P(atm)','Z','Zstd','ePot(kcal/mol)','eMol(kcal/mol)','&
				eVdw(kcal/mol)','IVdw(kcal/mol)','RunTime(ns)','EqTime(ns)','nMolec'
	write(61,'(A9,1x,A9,1x,3(A9,1x),6(A15,1x),A6)')'T(K)','rho(gcc)','P(atm)','Z','Zstd','ePot(kcal/mol)','eMol(kcal/mol)','&
				eVdw(kcal/mol)','IVdw(kcal/mol)','RunTime(ns)','EqTime(ns)','nMolec'
	close(61)

	open(51,file=logFile,ioStat=ioErr)
	if(ioErr==1)write(*,*)'Error opening logFile'

10	continue
	noNewData=1
	do while(noNewData.ne.0)
		read(51,'(a)',ioStat=ioErr,END=861)dumString
		read(dumString,*,ioStat=ioErr)firstWord
		if(ioErr.ne.0)cycle
		if(firstWord=='Step')noNewData=0
		if(trim(dumString) .eq. "run		${RunPR}" .OR. trim(dumString) .eq. "run		${Run}") then
			read(51,*)dumString2,nRunSteps
		endif
		if(trim(dumString) .eq. "timestep	${timestep}") then
			read(51,*)dumString2,timeStep
		endif
	enddo

	nData=0
	tAvg=0
	pAvg=0
	zAvg=0
	rhoAvg=0
	ePotAvg=0
	eBondAvg=0
	eVdwAvg=0
	ioErr=0
	do while(ioErr==0)
		read(51,'(a188)',ioStat=ioErr,END=861)dumString
		!write(*,*)dumString

		read(dumString,*,ioStat=ioErr)Step,Elapsed,Time,CPULeft,Temp,Press,PotEng,rKinEng,TotEng,E_vdwl,E_bond,E_coul,Volume,Density !, zFactor 
		if(ioErr.NE.0)cycle
		if(Temp < 1)cycle
		if(Elapsed < equil_percent*timeStep/100.d0)cycle
		nData=nData+1
		tAvg=tAvg+Temp
		pAvg=pAvg+Press
		zFactor=Press*Volume*0.0101325/(1.3806488*Temp*Nmolec)
		zAvg=zAvg+zFactor
		zSqAvg=zSqAvg+zFactor*zFactor
		rhoAvg=rhoAvg+Density
		ePotAvg=ePotAvg+PotEng
		eBondAvg=eBondAvg+E_bond
		eVdwAvg=eVdwAvg+E_vdwl
		!write(*,*)'iData,zFactor',nData,zFactor
	enddo

	if(nData < 11)goto 10
	tAvg=tAvg/nData
	pAvg=pAvg/nData
	zAvg=zAvg/nData
	zSqAvg=zSqAvg/nData
	stDevZ=SQRT( nData*(zSqAvg-zAvg*zAvg)/(nData-1) )
	rhoAvg=rhoAvg/nData
	ePotAvg=ePotAvg/nData
	eBondAvg=eBondAvg/nData
	eVdwAvg=eVdwAvg/nData
	simTimeNs=timeStep*nRunSteps/1d6

	open(61,file=outFile,ACCESS='APPEND')
	write(61,'(f9.2,1x,f9.5,1x,3(f9.3,1x),4(f15.3,1x),2(f15.1,1x),I6)')tAvg,rhoAvg,pAvg,zAvg,stDevZ,ePotAvg,eBondAvg, &
							eVdwAvg,AvgIntraVdw,simTimeNs,equil_percent/100.d0*simTimeNs,nMolec
	write(6 ,'(f9.2,1x,f9.5,1x,3(f9.3,1x),4(f15.3,1x),2(f15.1,1x),I6)')tAvg,rhoAvg,pAvg,zAvg,stDevZ,ePotAvg,eBondAvg, &
							eVdwAvg,AvgIntraVdw,simTimeNs,equil_percent/100.d0*simTimeNs,nMolec
	close(61)

	goto 10

861	continue

	write(*,*)'Success! Your data are stored in ',TRIM(outFile)
	stop
	end program LammpsMdRdr

subroutine IntraVdwCalc(AvgIntraVdw,nMolec)
	implicit double precision(A-H,O-Z)
	double precision :: x(20),y(20),z(20),dum,xL,yL,zL
	character :: dumpFile*88,IntraVdwOutFile*88,parFile*88,dumString*50
	integer :: timestep(100000),iType(20),nSitesPerMolec,iSite,iStep
	double precision,intent(out) :: AvgIntraVdw
	integer,intent(out) :: nMolec
	common parFile
	common dumpFile,IntraVdwOutFile

	dumpFile="dump.txt"
	IntraVdwOutFile="IntraVdw.res"
	parFile="TG.inp"

	open(1,file=dumpFile)
		read(1,*) 
		read(1,*) 
		read(1,*)
		read(1,*) nAtoms
		read(1,*)
		read(1,*)
		read(1,*)
		read(1,*)
		read(1,*)
		do i=1,nAtoms
			read(1,*) i1,i2
			if(i2.eq.2) then
				nAtomsPerMolec=i1-1
				goto 20
			endif
		enddo
	20	continue
	close(1)
	nMolec=nAtoms/nAtomsPerMolec
	nMaxSnap=5000000
	open(4,file=IntraVdwOutFile)
		write(4,*) nAtoms,nMolec,trim(dumpFile)," ",trim(IntraVdwOutFile)," ",trim(parFile)
	open(3,file=dumpFile)
	nn=0
	do iStep=1,nMaxSnap
		read(3,*,iostat=iEnd) dumString
		if(iEnd<0) exit
		read(3,*) timestep(iStep)
		read(3,*) dumString
		read(3,*) nAtoms
		read(3,*) dumString
		read(3,*) dum,xL
		read(3,*) dum,yL
		read(3,*) dum,zL
		read(3,*) dumString
		do iMolec=1,nMolec
			do iSite=1, nAtomsPerMolec
				read(3,*) dum,dum,iType(iSite),x(iSite),y(iSite),z(iSite)
			enddo
			do i=1,nAtomsPerMolec-1
				do j=i+1,nAtomsPerMolec

				if ((x(i)-x(j))>xL/2) then
					rxijMIC = x(i)-x(j)-xL
				elseif((x(i)-x(j))<(-1*xL/2)) then
					rxijMIC = x(i)-x(j)+xL
				else
					rxijMIC = x(i)-x(j)
				endif

				if ((y(i)-y(j))>yL/2) then
					ryijMIC = y(i)-y(j)-yL
				elseif((y(i)-y(j))<(-1*xL/2)) then
					ryijMIC = y(i)-y(j)+yL
				else
					ryijMIC = y(i)-y(j)
				endif

				if ((z(i)-z(j))>zL/2) then
					rzijMIC = z(i)-z(j)-xL
				elseif((z(i)-z(j))<(-1*zL/2)) then
					rzijMIC = z(i)-z(j)+zL
				else
					rzijMIC = z(i)-z(j)
				endif

				call IntraStatus(i,j,istatus)
				!write(*,*)i,j,istatus
				if(istatus==0) then 
					r=sqrt(rxijMIC**2+ryijMIC**2+rzijMIC**2)
					call pot(iType(i),iType(j),r,Energy)
					eIntraNonBond=eIntraNonBond+Energy
				endif
			
				enddo
			enddo
		enddo
		write(4,*)timestep(iStep),eIntraNonBond
		if(iStep.ne.1) then		
		AvgIntraVdw=AvgIntraVdw+eIntraNonBond
		nn=nn+1
		endif
		eIntraNonBond=0.0
	enddo
	AvgIntraVdw=AvgIntraVdw/dble(nn)
	close(3)
	close(4)
end subroutine IntraVdwCalc

subroutine pot(ii,jj,r,Energy)
	implicit double precision(A-H,O-Z)
	character :: parFile*88
	character :: potentialName*20,keyword*20,dumString*50
	character (len=15), dimension(15) :: SiteTypeNames
	double precision :: sigma(15),epsilon(15),AorN(15)
	double precision, intent(out) :: Energy
	double precision, intent(in) :: r
	integer, intent(in) :: ii,jj
	integer :: nSiteTypes,i,j,k,l
	common parFile
	open(2,file=parFile)
		read(2,*)dumString,potentialName
		read(2,*)dumString,nSiteTypes
		read(2,*)dumString,(SiteTypeNames(j),j=1,nSiteTypes)
		read(2,*)dumString,(sigma(j),j=1,nSiteTypes)
		read(2,*)dumString,(epsilon(j),j=1,nSiteTypes)
		read(2,*)dumString,(AorN(j),j=1,nSiteTypes)
		read(2,*)dumString,rIncrement
		read(2,*)dumString,rInit
	close(2)

	zero = 0.0d0
	Conv=0.001987204d0
	rCutMax = 3.0*MAXVAL(sigma)
	sig=(sigma(ii)+sigma(jj))/2.d0
	eps=sqrt(epsilon(ii)*epsilon(jj))
	A=sqrt(1+AorN(ii)*AorN(jj)+AorN(ii)+AorN(jj))-1
	rCut=3.d0*sig
	xmin = -0.00019*A**4 + 0.0105*A**3 - 0.0245*A**2 + 0.0337*A  + 0.793700526
	Cij  = -1.0/(xmin**6 +  A*xmin**7 -  A*xmin**4 - xmin**3)
	foffset=0.0
	offset=0.0
	foffset= Cij*eps*(12*(sig/rCut)**12 + 14*A*(sig/rCut)**14 - 8*A*(sig/rCut)**8 -6*(sig/rCut)**6)
	offset= Cij*eps*((sig/rCut)**12 +A*(sig/rCut)**14 - A*(sig/rCut)**8 -(sig/rCut)**6)
	Energy=0
	Force=0
	if(r<=rCut)then
		Energy=(Cij*eps*((sig/r)**12 +A*(sig/r)**14 - A*(sig/r)**8 -(sig/r)**6) + (r/(rCut)-1)*foffset - offset)*Conv
		Force=(Cij*eps*(12.d0*(sig/r)**12 +14.d0*A*(sig/r)**14 - 8.d0*A*(sig/r)**8 -6.d0*(sig/r)**6) - (r/(rCut))*foffset)*Conv/r
	else
		Energy=0
		Force=0
	endif

end subroutine pot


subroutine IntraStatusDihedral(ii,jj,istatus)
	!METHOD:This subroutine reads the dihedrals from Connectivity.txt 
	!	(One of the input files used by DataFileGenerator executable file)
	!	and deteremines if both i and j belong to one dihedral. 
	!	If they do, then they will not be counted in intra VDW energy.

	implicit double precision(A-H,O-Z)
	integer, intent(out) :: istatus
	integer, intent(in) :: ii,jj
	INTEGER, DIMENSION(50, 4) :: DihArray
	INTEGER :: iiArray(50),jjArray(50)
	character :: dumStr*88

	iiArray=0
	jjArray=0
	istatus=2
	if(ii<=0 .OR. jj<=0) goto 40
	dumStr="SOMEWORD"
	open(111,file="Connectivity.txt")
		do while(dumStr.ne."Dihedrals")	
			read(111,*) dumStr,nDihedrals
			if(nDihedrals.eq.0) then
				istatus=2 
				goto 40
			endif
				
		end do
		do i=1,nDihedrals
			read(111,*)dum,dum,DihArray(i,1),DihArray(i,2),DihArray(i,3),DihArray(i,4)
		enddo


	nCount=0
	do i=1,nDihedrals
		do j=1,4
			if(DihArray(i,j)==ii) then
				nCount=nCount+1
				iiArray(nCount)=i	
			endif
		enddo	
	enddo
	nCount=0
	do i=1,nDihedrals
		do j=1,4
			if(DihArray(i,j)==jj) then
				nCount=nCount+1
				jjArray(nCount)=i	
			endif
		enddo	
	enddo

	iCount=0
	do i=1,50
		do j=1,50
			if(iiArray(i).eq.jjArray(j) .AND. iiArray(i).ne.0 .AND. jjArray(j).ne.0) then
				iCount=iCount+1
			endif
		enddo	
	enddo
	if(iCount>0) istatus=1
	if(iCount==0) istatus=0
	40 continue
	close(111)
end subroutine IntraStatusDihedral

subroutine IntraStatusAngle(ii,jj,istatus)
	!METHOD:This subroutine reads the angles from Connectivity.txt 
	!	(One of the input files used by DataFileGenerator executable file)
	!	and deteremines if both i and j belong to one angle. 
	!	If they do, then they will not be counted in intra VDW energy.

	implicit double precision(A-H,O-Z)
	integer, intent(out) :: istatus
	integer, intent(in) :: ii,jj
	INTEGER, DIMENSION(50, 4) :: AngArray
	INTEGER :: iiArray(50),jjArray(50)
	character :: dumStr*88

	iiArray=0
	jjArray=0
	istatus=2
	if(ii<=0 .OR. jj<=0) goto 40
	dumStr="SOMEWORD"
	open(111,file="Connectivity.txt")
		do while(dumStr.ne."Angles")	
			read(111,*) dumStr,nAngles
			if(nAngles.eq.0) then
				istatus=2 
				goto 40
			endif
				
		end do
		do i=1,nAngles
			read(111,*)dum,dum,AngArray(i,1),AngArray(i,2),AngArray(i,3)	!,AngArray(i,4)
		enddo


	nCount=0
	do i=1,nAngles
		do j=1,3
			if(AngArray(i,j)==ii) then
				nCount=nCount+1
				iiArray(nCount)=i	
			endif
		enddo	
	enddo
	nCount=0
	do i=1,nAngles
		do j=1,3
			if(AngArray(i,j)==jj) then
				nCount=nCount+1
				jjArray(nCount)=i	
			endif
		enddo	
	enddo

	iCount=0
	do i=1,50
		do j=1,50
			if(iiArray(i).eq.jjArray(j) .AND. iiArray(i).ne.0 .AND. jjArray(j).ne.0) then
				iCount=iCount+1
			endif
		enddo	
	enddo
	if(iCount>0) istatus=1
	if(iCount==0) istatus=0
	40 continue
	close(111)
end subroutine IntraStatusAngle

subroutine IntraStatusBond(ii,jj,istatus)
	!METHOD:This subroutine reads the bonds from Connectivity.txt 
	!	(One of the input files used by DataFileGenerator executable file)
	!	and deteremines if both i and j belong to one bond. 
	!	If they do, then they will not be counted in intra VDW energy.

	implicit double precision(A-H,O-Z)
	integer, intent(out) :: istatus
	integer, intent(in) :: ii,jj
	INTEGER, DIMENSION(50, 4) :: AngArray
	INTEGER :: iiArray(50),jjArray(50)
	character :: dumStr*88

	iiArray=0
	jjArray=0
	istatus=2
	if(ii<=0 .OR. jj<=0) goto 40
	dumStr="SOMEWORD"
	open(111,file="Connectivity.txt")
		do while(dumStr.ne."Bonds")	
			read(111,*) dumStr,nBonds
			if(nBonds.eq.0) then
				istatus=2 
				goto 40
			endif
				
		end do
		do i=1,nBonds
			read(111,*)dum,dum,AngArray(i,1),AngArray(i,2)	!,AngArray(i,3)	!,AngArray(i,4)
		enddo


	nCount=0
	do i=1,nBonds
		do j=1,2
			if(AngArray(i,j)==ii) then
				nCount=nCount+1
				iiArray(nCount)=i	
			endif
		enddo	
	enddo
	nCount=0
	do i=1,nBonds
		do j=1,2
			if(AngArray(i,j)==jj) then
				nCount=nCount+1
				jjArray(nCount)=i	
			endif
		enddo	
	enddo

	iCount=0
	do i=1,50
		do j=1,50
			if(iiArray(i).eq.jjArray(j) .AND. iiArray(i).ne.0 .AND. jjArray(j).ne.0) then
				iCount=iCount+1
			endif
		enddo	
	enddo
	if(iCount>0) istatus=1
	if(iCount==0) istatus=0
	40 continue
	close(111)
end subroutine IntraStatusBond

subroutine IntraStatus(ii,jj,istatus)
	!DETERMINES WHETHER SITE II AND SITE JJ ARE AT LEAST FOUR BONDS APART
	!IF ISTATUS=0 THEN I AND J ARE AT LEAST FOUR BONDS APART
	!IF ISTATUS=1 THEN I AND J ARE NOT AT LEAST FOUR BONDS APART

	!METHOD:This subroutine calls these subroutines, 
	! 	IntraStatusBond(ii,jj,i1)
	! 	IntraStatusAngle(ii,jj,i2)
	! 	IntraStatusDihedral(ii,jj,i3)
	!	and if i1,i2 and i3 are zero, returns 0
	!	otherwise returns 1

	implicit double precision(A-H,O-Z)
	integer, intent(out) :: istatus
	integer, intent(in) :: ii,jj

	call IntraStatusBond(ii,jj,i1)
	call IntraStatusAngle(ii,jj,i2)
	call IntraStatusDihedral(ii,jj,i3)
	
	if(i1==0.AND.i2==0.AND.i3==0) then 
		istatus=0
	else
		istatus=1
	endif

end subroutine IntraStatus


