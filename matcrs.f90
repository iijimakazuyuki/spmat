module matcrs_mod
	implicit none
	type matcrs
		double precision, allocatable :: e(:)
		integer :: n, m
		integer, allocatable :: idx(:), col(:)
	end type
	
	!疎行列ベクトル積
	interface operator(*)
		module procedure spmatvec
	end interface
	
	!CRS形式から密行列に変換
	interface assignment(=)
		module procedure crs2dense
	end interface
contains
	!疎行列ベクトル積
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
	
	!CRS形式から密行列に変換
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
	
	!密行列を出力
	subroutine print_matdense(mat)
		double precision, intent(in) :: mat(:,:)
		integer :: i
		do i=1, size(mat(:,1))
			print *, mat(i,:)
		end do
	end subroutine
	
	!CRS形式の疎行列を出力 (Matrix Market coordinate format)
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
	
	!CRS形式の疎行列を出力 (original format)
	subroutine print_matcrs_array(mat)
		type(matcrs), intent(in) :: mat
		integer :: i, j
		print *, mat%n, mat%m
		print *, mat%e
		print *, mat%idx
		print *, mat%col
	end subroutine
	
	!CRS形式の疎行列を密行列形式で出力
	subroutine print_matcrs_2d(mat)
		type(matcrs), intent(in) :: mat
		double precision, allocatable :: dense(:,:)
		dense = mat
		call print_matdense(dense)
	end subroutine
	
	!CRS形式の疎行列を読み込む (original format)
	function read_matcrs_array() result(mat)
		type(matcrs) :: mat
		read *, mat%n, mat%m
		allocate(mat%e(mat%m), mat%idx(0:mat%n), mat%col(mat%m))
		read *, mat%e, mat%idx, mat%col
	end function
	
	subroutine init_matcrs(mat)
		type(matcrs), intent(in) :: mat
		allocate(mat%e(mat%m), mat%idx(0:mat%n), mat%col(mat%m))
	end subroutine
end module