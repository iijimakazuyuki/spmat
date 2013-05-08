module matcrs_mod
	implicit none
	type matcrs
		double precision, pointer :: e(:)
		integer :: n, m
		integer, pointer :: idx(:), col(:)
	end type
	
	!asñxNgÏ
	interface operator(*)
		module procedure spmatvec
	end interface
	
	!CRS`®©ç§sñÉÏ·
	interface assignment(=)
		module procedure crs2dense
	end interface
contains
	!asñxNgÏ
	function spmatvec(mat,vec) result(sp)
		double precision, allocatable :: sp(:)
		type(matcrs), intent(in) :: mat
		double precision, intent(in) :: vec(:)
		integer :: i, j
		allocate(sp(mat%n))
		sp = 0
		!$omp do private(j)
		do i=1, mat%n
			do j=mat%idx(i-1)+1, mat%idx(i)
				sp(i) = sp(i) + mat%e(j)*vec(mat%col(j))
			end do
		end do
		!$omp end do
	end function
	
	!CRS`®©ç§sñÉÏ·
	subroutine crs2dense(dense,crs)
		double precision, allocatable, intent(out) :: dense(:,:)
		type(matcrs), intent(in) :: crs
		integer :: i, j
		allocate(dense(crs%n,maxval(crs%col)))
		dense = 0
		do i=1, crs%n
			do j=crs%idx(i-1)+1, crs%idx(i)
				dense(i,crs%col(j)) = crs%e(j)
			end do
		end do
	end subroutine
	
	!§sñðoÍ
	subroutine print_matdense(mat)
		double precision, intent(in) :: mat(:,:)
		integer :: i
		do i=1, size(mat(:,1))
			print *, mat(i,:)
		end do
	end subroutine
	
	!CRS`®ÌasñðoÍ (Matrix Market coordinate format)
	subroutine print_matcrs(mat)
		type(matcrs), intent(in) :: mat
		integer :: i, j
		print *, mat%n, maxval(mat%col), mat%m
		do i=1, mat%n
			do j=mat%idx(i-1)+1, mat%idx(i)
				print *, i, mat%col(j), mat%e(j)
			end do
		end do
	end subroutine
	
	!CRS`®ÌasñðoÍ (original format)
	subroutine print_matcrs_array(mat)
		type(matcrs), intent(in) :: mat
		integer :: i, j
		print *, mat%n, mat%m
		print *, mat%e
		print *, mat%idx
		print *, mat%col
	end subroutine
	
	!CRS`®Ìasñð§sñ`®ÅoÍ
	subroutine print_matcrs_2d(mat)
		type(matcrs), intent(in) :: mat
		double precision, allocatable :: dense(:,:)
		dense = mat
		call print_matdense(dense)
	end subroutine
	
	!CRS`®ÌasñðÇÝÞ (original format)
	function read_matcrs_array(e, idx, col) result(mat)
		type(matcrs) :: mat
		double precision, target, allocatable, intent(inout) :: e(:)
		integer, target, allocatable, intent(inout) :: idx(:), col(:)
		read *, mat%n, mat%m
		call init_matcrs(mat, e, idx, col)
		read *, e, idx, col
	end function
	
	subroutine init_matcrs(mat, e, idx, col)
		type(matcrs) :: mat
		double precision, target, allocatable, intent(inout) :: e(:)
		integer, target, allocatable, intent(inout) :: idx(:), col(:)
		allocate(e(mat%m), idx(0:mat%n), col(mat%m))
		e = 0
		idx = 0
		col = 0
		mat%e => e
		mat%idx => idx
		mat%col => col
	end subroutine
	
end module