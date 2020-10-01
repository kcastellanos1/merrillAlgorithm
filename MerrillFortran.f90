
!!!  Merrill subroutine

subroutine merrill( jout, ps, labels, iv, obj, start, v, verbose, refinement)

    ! This subroutine calculates the correct labels and replacement operations
    ! for Merrill's algorithm to solve a v-variable system of equations.
    ! Based off Merrill's algorithm which is a variant of Herb Scarf's 
    ! fixed point iteration algorithm.
    
    ! jout
    ! ps ...
    ! labels ...
    ! iv ...
    ! obj ...
    ! start = 1 if starting primitive simplex
    ! v ....

    integer v, start, jout, c, r, ind, p, i2, x1, major, refinement
    integer xt, jout1, jout2, i3, verbose
    double precision, dimension( v+1, v+1) :: ps
    integer, dimension(v+1) :: labels
    double precision, dimension(v) :: obj, iv 
    
    if(jout > v+1)then
        jout = 1
        write(*,*) "Err. in merrill subroutine:"
        write(*,*) "jout out of bounds - resetting to 1"
    end if
    
    ! Minor iteration
    
    ! Find duplicate label
    
    ! If the algorithm has exited and started again - i.e.
    ! the subroutine is called in the program then we know
    ! we are on the artificial simplex - use Scarf labeling
    ! Find first entry in objective vector that is positive
    ind = 0
    do r = 0, v-1
        if (ps(v+1-r,jout) == 0) then
            ind = v-r
            exit
        else if (obj(v-r)>0) then
            ind = v-r ! Index of positive entry
        end if
    end do
    
    labels(jout) = ind ! Set the label
    
    ! Write out simplex for troubleshooting
    if (verbose == 1) then
        write(*,*) "----------Prim. Simplex---------"
        do r = 1,v+1
            write(*,*) ps(r,:)
        end do
        write(*,*) "--------Labeling ", jout," as ", ind, "--------"
        write(*,*) "Labels: ", labels
    end if    
    
    ! Now we that we have the correct labels we can go
    ! back into the original algorithm. Remember the
    ! outside program is just generating the excess demands
    ! or objective function values. Once we have calculated
    ! excess demands we feed it back into the subroutine 'merrill'
    
    ! Exit flag -- This means we have a simplex that is either
    ! comletely labeled or
    ! in need of excess demands to be calculated
    major = 0
    xt = 0
    
    if (start == 1) then
        xt = 1
        major = 1
    end if
    
    do while (xt == 0)
    
    
        ! See if simplex is completely labeled
        x1 = 0
        if ( sum(ps(1,:)) == 1 ) then
            x1 = 1 ! First condition
        end if
        
        
        ! We have a ps with only one vertex (x1) on the artificial plane
        ! now we need to see if the simplex on the real plane is
        p = 0
        if (x1 == 1) then
            do i2 = 1, v
                do i3 = 1, v + 1
                    ! Second condition below here
                    if ( labels ( i3 ) == i2 .and. ps ( 1, i3 ) == 0 ) then
                        p = p + 1
                        exit
                    end if
                end do
            end do
        end if
        
        if (p == v .and. x1 == 1) then
            ! Simplex is completely labeled - move to major iteration
            !write(*,*) "--------------------------------------"
            !write(*,*) "            Major Iteration", sum(abs(obj))
            !write(*,*) "--------------------------------------"
            major = 1
            xt = 1
            do i2 = 2, v+1
                iv(i2-1) = ps(i2,jout)*refinement
                !write(*,*) "Vertex: ", iv(i2-1), "Value: ", obj(i2-1)
            end do
        end if
        
        
        ! If simplex is not completely labeled
        ! We need to replace the duplicate column
        do r = 1, v + 1
            if (labels(jout) == labels(r) .and. jout /= r) then
                jout = r ! Index of duplicate label
                exit
            end if
        end do
        
        if(verbose == 1)then
            write(*,*) "Replacing column ", jout
        end if
        
        ! Now we make the replacement
            jout1 = jout - 1
            jout2 = jout + 1
            if ( jout1 == 0 ) then
                jout1 = v + 1
            end if
            if ( jout2 > ( v + 1 ) ) then
                jout2 = 1
            end if
            
            ! replacement operation
            ps( :, jout ) = ps( :, jout2 ) + ps( :, jout1 ) - ps( :, jout )
            
            ! If we are on the artificial simplex calculate label
            if (ps(1,jout) == 1) then
                do r = 1, v
                    if (ps(r+1,jout)<iv(r)) then
                        labels(jout) = r
                        exit
                    end if
                end do
            else
                xt = 1
            end if
            
                            ! Write out simplex for troubleshooting
            if (verbose == 1) then
                write(*,*) "----------New Prim. Simplex---------"
                do r = 1,v+1
                    write(*,*) ps(r,:)
                end do
                write(*,*) "Labels: ", labels
            end if  
            
        
    end do
    
    
    
    ! ------------------------------------------------------------ Major iteration
    
    ! If this is the major iteration (or first time through) we
    ! need to setup the new primitive simplex from the initial
    ! vertex, as well as initialize the labeling vector.
    
    if (major == 1) then
    
        do c = 1, v
        
            labels(c+1) = c ! Setup Labels
            
            do r = 1, v
            
                ps( r+1, c+1) = iv(r) ! Setup primitive simplex
                
                if (c == r) then
                    ps( r+1, c+1) = iv(r) - 1 ! Diagonal entry shifted
                end if
            
            end do    
            
        end do
        
        
        do r = 1, v
            ps( r+1, 1) = iv(r) ! Set first column up as initial vertex
            ps( 1, r+1) = 1 ! Set first row as artifial simplex
            labels(r+1) = r
        end do
        
        ps(1,1) = 0
        jout = 1
    
    end if

    if(verbose == 1)then
        write(*,*) "-----------------"
        write(*,*) "Exiting Merrill"
        write(*,*) "Jout = ", jout
        write(*,*) "-----------------"
    end if
    
end subroutine merrill