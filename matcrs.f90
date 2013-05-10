module matcrs_mod
	implicit none
	type matcrs
		double precision, allocatable :: e(:)
		integer :: n, m
		integer, allocatable :: idx(:), col(:)
	end type
	
	type matcrs_part
		type(matcrs) :: mat
		integer :: inn, neib
		integer, allocatable :: map(:)
	end type
	
	!疎行列ベクトル積
	interface operator(*)
		module procedure spmatvec
	end interface
	
	!CRS形式から密行列に変換
	interface assignment(=)
		module procedure crs2dense
	end interface
	
	private
	public matcrs, matcrs_part
	public spmatvec, crs2dense, read_matcrs_array, init_matcrs, part_matcrs
	public print_matdense, print_matcrs, print_matcrs_array, print_matcrs_2d, print_matcrs_part_array
	public write_matcrs_part_array
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
		call init_matcrs(mat)
		read *, mat%e, mat%idx, mat%col
	end function
	
	subroutine init_matcrs(mat)
		type(matcrs), intent(inout) :: mat
		allocate(mat%e(mat%m), mat%idx(0:mat%n), mat%col(mat%m))
		mat%idx(0) = 0
	end subroutine
	
	!CRS形式の疎行列を分割
	function part_matcrs(mat, part, npp) result(smat)
		type(matcrs_part) :: smat
		type(matcrs), intent(in) :: mat
		integer, intent(in) :: part(npp), npp
		
		integer :: i, j, inv(mat%n)
		logical :: interior(mat%n), exterior(mat%n)
		
		smat%mat%n = npp
		smat%mat%m = 0
		do i=1, npp
			smat%mat%m = smat%mat%m + mat%idx(part(i)) - mat%idx(part(i)-1)
		end do
		call init_matcrs(smat%mat)
		
		do i=1, smat%mat%n
			smat%mat%idx(i) = mat%idx(part(i)) - mat%idx(part(i)-1) + smat%mat%idx(i-1)
			smat%mat%e(smat%mat%idx(i-1)+1:smat%mat%idx(i)) = mat%e(mat%idx(part(i)-1)+1:mat%idx(part(i)))
			smat%mat%col(smat%mat%idx(i-1)+1:smat%mat%idx(i)) = mat%col(mat%idx(part(i)-1)+1:mat%idx(part(i)))
		end do
		
		interior = .false.
		exterior = .false.
		do i=1, smat%mat%n
			do j=smat%mat%idx(i-1)+1, smat%mat%idx(i)
				if(has(part, smat%mat%col(j))) then
					interior(smat%mat%col(j)) = .true.
				else
					exterior(smat%mat%col(j)) = .true.
				end if
			end do
		end do
		
		smat%inn = count(interior)
		smat%neib = count(exterior)
		allocate(smat%map(smat%inn + smat%neib))
		
		inv = 0
		do i=1, smat%inn
			j = indexof(interior, .true.)
			smat%map(i) = j
			inv(j) = i
			interior(j) = .false.
		end do
		
		do i=1, smat%neib
			j = indexof(exterior, .true.)
			smat%map(i+smat%inn) = j
			inv(j) = i+smat%inn
			exterior(j) = .false.
		end do
		
		do i=1, smat%mat%m
			smat%mat%col(i) = inv(smat%mat%col(i))
		end do
		
	end function
	
	!
	subroutine print_matcrs_part_array(mat)
		type(matcrs_part) :: mat
		call print_matcrs_array(mat%mat)
		print *, mat%inn, mat%neib
		print *, mat%map
	end subroutine
	
	subroutine write_matcrs_part_array(n, mat)
		integer :: n
		type(matcrs_part) :: mat
		write(n,*) mat%mat%n, mat%mat%m
		write(n,*) mat%mat%e
		write(n,*) mat%mat%idx
		write(n,*) mat%mat%col
		write(n,*) mat%inn, mat%neib
		write(n,*) mat%map
	end subroutine
	
	logical function has(a, x)
		integer, intent(in) :: a(:), x
		integer :: i
		has = .false.
		do i=1, size(a)
			if(a(i) == x) then
				has = .true.
				exit
			end if
		end do
	end function
	
	integer function indexof(a, x)
		logical, intent(in) :: a(:), x
		integer :: i
		indexof = -1
		do i=1, size(a)
			if(a(i) .eqv. x) then
				indexof = i
				exit
			end if
		end do
	end function
	
end module