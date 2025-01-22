! Program to find the 'chemical' equitable partition and corresponding quotient
! graph of a graph, i.e. the partition based on the chemical elements and 
! valences of the atoms or groups of atoms represented by graph nodes.
! 
! Copyright (c) 2022-2024 Vasilios Raptis <v.raptis@external.euc.ac.cy>
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! 
!-------------------------------------------------------------------------------
! 
program cep
    use strings
    implicit none
    integer, parameter :: INPUT=10,QGDOT=11,QGDAT=12
    integer, parameter :: ntrymax=1000
    type atom_t
        integer      cell
        integer      chemtype
        integer      idx
        integer      valence
        character(2) element
    end type atom_t
    integer itry
    logical check,exists
    integer, allocatable :: oldprof(:,:),profile(:,:)
    integer, allocatable :: adjacency(:,:),qgadj(:,:)
    type (atom_t), allocatable :: atom(:)
    character(RECLEN) inpfil
    ! Provide input file as commandline argument 
    CALL get_command_argument(1,inpfil)
    if(strlen(inpfil)==0) then
        write(*,'("No input filename provided!")')
        STOP
    endif
    inquire(file=inpfil(1:strlen(inpfil)),exist=exists)
    if(.NOT. exists) then
        write(*,'("Requested input file does not exist!")')
        STOP
    endif
    ! Step 1 
    ! DOT file parser 
    CALL parse_input(inpfil)
    ! Step 2
    ! Now, initialise the chemical cells partition 
    CALL initialise()
    ! The party is about to begin 
    itry=0
    check=.TRUE.
    do while(check)
        itry=itry+1
        if(itry>ntrymax) then
            write(*,'("Failed to converge after",I4," trials. Exiting.")')itry
            STOP
        endif
        ! Step 3
        ! Save current partition to compare with next one
        CALL backup()
        ! Step 4
        ! Now, reconstruct partition based on current atom profiles 
        ! First, find new cells 
        CALL find_new_cells()
        ! Step 5
        ! Then, allocate atoms in new cells 
        CALL partition()
        ! Step 6 
        ! Check convergence; it's a two-step procedure
        ! 6-1. If profile and oldprof have different sizes, not converged yet
        ! 6-2. If profile and oldprof have same size, compare elements 
        check=check_convergence()
    enddo
    write(*,'(/"*** Chemical Equitable Partition found! ***")')
    ! Done! Form adjacency matrix of CEP 
    CALL form_qgadj()
    ! write quotient graph to file qg.dot 
    CALL write_qg()
    ! write node connectivity profiles, partition cells, and qg-adjacency matrix
    CALL write_part()
contains
    subroutine parse_input(inpfil)    
        implicit none
        integer i,m,n,nchemtyp,nodes,p,q
        character(2) atlabel
        logical read_edges,read_nodes
        character(RECLEN) inpfil
        type (atom_t), allocatable :: tmpat(:) 
        !-------------------------------------
        read_edges=.FALSE.
        read_nodes=.FALSE.
        open (unit=INPUT,file=inpfil(1:strlen(inpfil)))
        nodes=0
        nchemtyp=0
        do while(get_record(INPUT))
            if(strstr("comment",record)>0) then
                if(strstr("nodes section",record)>0) then
                    read_edges=.FALSE.
                    read_nodes=.TRUE.
                    if(get_record(INPUT)) &
                        continue    ! read next line 
                else if(strstr("edges section",record)>0) then
                    read_edges=.TRUE.
                    read_nodes=.FALSE.
                    if(get_record(INPUT)) &
                        continue    ! read next line 
                endif
            endif
            ! First, save the nodes 
            if(read_nodes) then
                !--- backup atoms 
                if(nodes>0 .AND. allocated(atom)) then
                    if(allocated(tmpat)) &
                        deallocate(tmpat)
                    allocate(tmpat(nodes))
                    tmpat=atom
                    deallocate(atom)
                endif
                !--- read info
                read(record,*,err=10)i
                p=strstr("label",record)    ! assumed format: compulsory for each node record to have a label
                if(p>0) then
                    p=p+7   ! length of `label="`
                    read(record(p:),*)atlabel
                endif
                nodes=nodes+1
                !--- extend atoms array and restore atoms backed up
                if(allocated(atom)) &
                    deallocate(atom)
                allocate(atom(nodes))
                n=nodes-1
                if(n>0) then
                    atom(1:n)=tmpat
                    deallocate(tmpat)
                endif
                !--- store/initialise info
                atom(nodes)%idx=i
                atom(nodes)%element=atlabel
                atom(nodes)%valence=0
                atom(nodes)%chemtype=0
typloop:        do i=1,n
                    if(atom(nodes)%element==atom(i)%element) then
                        atom(nodes)%chemtype=atom(i)%chemtype
                        exit typloop
                    endif
                enddo typloop
                if(atom(nodes)%chemtype==0) then
                    nchemtyp=nchemtyp+1
                    atom(nodes)%chemtype=nchemtyp
                endif
10              continue
            endif
            ! Then, save the edges 
            if(read_edges) then
                !--- read info
                ! node-node connection is either an arc or an edge; we can read  
                ! any of them by looking for the corresponding string
                p=strstr("--",record)
                p=max(p,strstr("->",record)) ! it can't be both so one of them will be zero
                q=p+2   ! for number after edge/arc symbol 
                p=p-1   ! for number before edge/arc symbol
                read(record(1:p),*,end=11)m
                read(record(q: ),*,end=11)n
                p=0
                q=0
                do i=1,size(atom)
                    if(atom(i)%idx==m) &
                        p=p+1
                    if(atom(i)%idx==n) &
                        q=q+1
                enddo
                i=0
                if(p==0) &
                    i=i+1
                if(q==0) &
                    i=i+1
                !--- backup atoms 
                if(allocated(atom)) then
                    if(allocated(tmpat)) &
                        deallocate(tmpat)
                    allocate(tmpat(nodes))
                    tmpat=atom
                    deallocate(atom)
                endif
                !--- extend atoms array and restore atoms backed up
                if(allocated(atom)) &
                    deallocate(atom)
                if(allocated(tmpat)) then
                    allocate(atom(size(tmpat)+i))
                    i=size(tmpat)
                else
                    i=0
                endif
                if(i>0) then
                    atom(1:i)=tmpat
                    deallocate(tmpat)
                endif
                !--- store/initialise info
                if(p==0) then
                    i=i+1
                    atom(i)%idx=m
                endif
                if(q==0) then
                    i=i+1
                    atom(i)%idx=n
                endif
            endif
        enddo
11      continue    
        ! 
        !--- 
        ! Now that all nodes in the system have been registered 
        ! it is safe to identify adjacent pairs and update 
        ! adjacency matrix, accordingly 
        rewind(unit=INPUT)
        read_nodes=.FALSE.
        read_edges=.FALSE.
        n=size(atom)
        allocate(adjacency(n,n))
        adjacency=0
        do while(get_record(INPUT))
            if(strstr("comment",record)>0) then
                if(strstr("edges section",record)>0) then
                    read_edges=.TRUE.
                    if(get_record(INPUT)) &
                        continue    ! go to next line 
                endif
            endif
            if(read_edges) then
                !--- read info
                ! node-node connection is either an arc or an edge; we can read  
                ! any of them by looking for the corresponding string
                p=strstr("--",record) ! can't be both so one of them will be zero
                p=max(p,strstr("->",record))
                q=p+2   ! for number after edge symbol 
                p=p-1   ! for number before edge symbol
                read(record(1:p),*,end=21)m
                read(record(q: ),*,end=21)n
                adjacency(m,n)=adjacency(m,n)+1
                adjacency(n,m)=adjacency(n,m)+1
                atom(m)%valence=atom(m)%valence+1
                atom(n)%valence=atom(n)%valence+1
            endif
        enddo
21      close(unit=INPUT)
        return
    end subroutine parse_input
!---
    subroutine initialise()
        implicit none
        integer i,nat
        nat=size(atom)
        allocate(profile(2,nat))
        do i=1,nat
            atom(i)%cell=0
            profile(:,i)=(/atom(i)%valence,atom(i)%chemtype/)
        enddo
        return
    end subroutine initialise 
!---
    subroutine partition()
        implicit none
        integer i,j,k,nat,ncell
        nat=size(atom)
        ncell=maxval((/(atom(i)%cell,i=1,nat)/))
        write(*,'(1X,"ncell=",I2)')ncell
        if(allocated(profile)) &   
            deallocate(profile) 
        allocate(profile(1+ncell,nat))          !   take into account both chemical element and valence 
        profile=0
        do i=1,nat
            profile(1+ncell,i)=atom(i)%chemtype !   take into account both chemical element and valence
            do j=i+1,nat
                if(adjacency(i,j)==1) then
                    k=atom(j)%cell
                    profile(k,i)=profile(k,i)+1
                    k=atom(i)%cell
                    profile(k,j)=profile(k,j)+1
                endif
            enddo
        enddo
        return
    end subroutine partition
!---
    subroutine backup()
        implicit none
        integer i,nat,ncell
        nat=size(atom)
        ncell=maxval((/(atom(i)%cell,i=1,nat)/))
        if(allocated(oldprof)) &   
            deallocate(oldprof) 
        allocate(oldprof(ncell,nat))
        oldprof=profile
        return
    end subroutine backup
!---
    subroutine find_new_cells()
        implicit none
        integer i,im1,j,k,nat,Nc
        logical found
        nat=size(atom)
        atom(1)%cell=1
        Nc=1
        do i=2,nat
            im1=i-1
inner:      do j=1,im1
                found=.TRUE.     
                do k=1,size(profile(:,1))
                    found=found.AND.(profile(k,i)==profile(k,j))
                enddo
                if(found) &
                    exit inner
	        enddo inner
	        if(found) then
  	            atom(i)%cell=atom(j)%cell
  	        else
	        	Nc=Nc+1
          	    atom(i)%cell=Nc
  	        endif
        enddo
        return
    end subroutine find_new_cells
!---
    logical function check_convergence()
        implicit none
        integer i,j,nat
        logical :: fail=.TRUE.
        if(size(profile(:,1))==size(oldprof(:,1))) then
            fail=.FALSE.
            nat=size(Atom)
check:      do i=1,nat
                do j=1,size(profile(:,1))
                    if(profile(j,i)/=oldprof(j,i)) then
                        fail=.TRUE.
                        exit check
                    endif
                enddo
            enddo check
        endif
        check_convergence=fail
        return
    end function check_convergence    
!---
    subroutine form_qgadj()
        implicit none
        integer i,j,m,n,nat,ncell
        nat=size(atom)
        ncell=maxval((/(atom(i)%cell,i=1,nat)/))
        allocate(qgadj(ncell,ncell))
        qgadj=0
        do i=1,nat
            do j=1,nat
                m=atom(i)%cell
                n=atom(j)%cell
                qgadj(m,n)=qgadj(m,n)+adjacency(i,j)
            enddo
        enddo
        ! Now, it contains sum of edges from cell m to cell n 
        ! Divide by number of nodes in each cell to obtain CEP adj mat; 
        ! it works because we know it is CEP so a balance of the form 
        !   (# edges from m to n)*(# atoms in m) = (# edges from n to m)*(# atoms in n)
        ! should hold. 
        do m=1,ncell
            n=0
            do i=1,nat
                if(atom(i)%cell==m) &
                    n=n+1
            enddo
            qgadj(m,:)=qgadj(m,:)/n
        enddo
        return
    end subroutine form_qgadj    
!---
    subroutine write_qg()
        implicit none
        integer i,j,n,nat,ncell,nfmt
        character(RECLEN) fmtrec
        nat=size(atom)
        n=size(qgadj(:,1))
        ncell=maxval((/(atom(i)%cell,i=1,nat)/))
        open (unit=QGDOT,file="qg.dot")
        write(QGDOT,'("digraph G {"/"layout=""dot"""/"graph [shape=""square""];")')
        write(QGDOT,'("comment=""nodes section""")'   )
        write(QGDOT,'("node  [style=""filled""; shape=""circle""];")')
        do i=1,ncell
inner:      do j=1,nat
                if(atom(j)%cell==i) &
                    exit inner
            enddo inner
            write(QGDOT,'(I5," [fillcolor=""grey""   label=""",A2",",I2,"""]")')i,atom(j)%element,i
        enddo
        write(QGDOT,'("comment=""edges section""")'   )
        nfmt=max(n,maxval(qgadj))
        nfmt=ceiling(log10(real(nfmt)))
        write(fmtrec,'("I",I1)')min(9,nfmt+1)
        do i=1,ncell
            do j=1,ncell
                if(qgadj(i,j)/=0) then
                    if(i==j) then
                        write(QGDOT,'(4X,'//fmtrec//',1X,"->",'//fmtrec//',1X,"[ label=""",'//fmtrec//   &
                                     ',""" color=""red"" ];")')i,j,qgadj(i,j)
                    else
                        write(QGDOT,'(4X,'//fmtrec//',1X,"->",'//fmtrec//',1X,"[ label=""",'//fmtrec//   &
                                     ',"""];")')i,j,qgadj(i,j)
                    endif
                endif
            enddo
        enddo
        write(QGDOT,'("}")')
        close(unit=QGDOT)
        return
    end subroutine write_qg
!---
    subroutine write_part
        implicit none
        integer i,j,m,n,nat,ncell,nfmt,np
        real CR,infocont
        character(RECLEN) fmtrec
        real,allocatable :: prob(:)
        open (unit=QGDAT,file="qg.dat")
        write(QGDAT,'(/"Converged connectivity profiles:")')
        nat=size(atom)
        n=size(qgadj(:,1))
        nfmt=ceiling(log10(real(nat)))
        write(fmtrec,'("I",I1)')min(9,nfmt+1)
        do i=1,nat
            np=size(profile(:,1))
            write(QGDAT,'(1X,"Node: ",'//fmtrec//',"  Type: ",'//fmtrec//       &
                        ',"  Label: ",A2"  Cell: ",'//fmtrec//',"  Profile: ")',&
            advance='no') i,profile(np,i),atom(i)%element,atom(i)%cell
            do j=1,np-1
                write(QGDAT,'('//fmtrec//',1X)',advance='no')profile(j,i)
            enddo
            write(QGDAT,'("")')
        enddo
        write(QGDAT,'(/"Converged partition cells:")')
        ncell=maxval((/(atom(i)%cell,i=1,nat)/))
        nfmt=max(ncell,nat)
        nfmt=ceiling(log10(real(nat)))
        allocate(prob(ncell))
        prob=0.0
        write(fmtrec,'("I",I1)')min(9,nfmt+1)
        do i=1,ncell
            write(QGDAT,'(1X,"Cell: ",'//fmtrec//',3X,"Nodes: ")',advance='no')i
            do j=1,nat
                if(atom(j)%cell==i) then
                    write(QGDAT,'('//fmtrec//',1X)',advance='no')j
                    prob(i)=prob(i)+1.0
                endif
            enddo            
            write(QGDAT,'("")')
        enddo
        prob=prob/real(nat)
        ! QG matrix 
        write(QGDAT,'(/"Quotient-graph adjacency matrix:")')
        nfmt=max(ncell,maxval(qgadj))
        nfmt=ceiling(log10(real(nfmt)))
        write(fmtrec,'("I",I1)')min(9,nfmt+1)
        m=0
        do i=1,ncell
            do j=1,ncell
                write(QGDAT,'('//fmtrec//',1X)',advance="no")qgadj(i,j)
                if(qgadj(i,j)/=0) &
                    m=m+1
            enddo
            write(QGDAT,'("")')
        enddo
        ! Summary: number of cells, %CR, information content; first nr of cells
        nfmt=ceiling(log10(real(ncell)))
        write(fmtrec,'("I",I1)')min(9,nfmt+1)
        write(QGDAT,'(/  "Number of cells",17X,":",3X,'//fmtrec//')')ncell
        write(*    ,'(4X,"Number of cells",17X,":",3X,'//fmtrec//')')ncell
        ! Compression ratio 
        CR=real(ncell)/real(nat)        
        CR=CR*real(m)
        m=0
        do i=1,nat
            do j=i+1,nat
                if(adjacency(i,j)/=0) &
                    m=m+1
            enddo
        enddo
        CR=100.0*CR/real(m)
        write(QGDAT,'(   "% compression ratio",13X,":",F8.2)')CR
        write(*    ,'(4X,"% compression ratio",13X,":",F8.2)')CR
        ! Normalised information content 
        infocont=0.0
        do i=1,ncell
            infocont=infocont-prob(i)*log(prob(i))
        enddo
        infocont=infocont/log(2.0)  ! log base 2
        write(QGDAT,'(   "Information content",13X,":",F8.2)')infocont
        write(*    ,'(4X,"Information content",13X,":",F8.2)')infocont
        infocont=infocont/(log(real(nat))/log(2.0))   ! normalise 
        write(QGDAT,'(   "% normalised information content:",F8.2)')infocont*100.0
        write(*    ,'(4X,"% normalised information content:",F8.2)')infocont*100.0
        close(unit=QGDAT)
        return
    end subroutine write_part
end program cep

