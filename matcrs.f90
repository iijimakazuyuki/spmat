module matcrs_mod
	implicit none
	
	!CRS
	type matcrs
		double precision, allocatable :: e(:)
		integer :: n, m, k
		integer, allocatable :: row(:), col(:)
	end type
	
	!
	type matcrs_part
		type(matcrs) :: mat
		integer :: inn, ext
		integer, allocatable :: map(:)
	end type
	
	!�a�s��x�N�g����
	interface operator(*)
		module procedure spmatvec
	end interface
	
	!CRS�`�����疧�s��ɕϊ�
	interface assignment(=)
		module procedure crs2dense
	end interface
	
	interface indexof
		module procedure indexof_b, indexof_i
	end interface
	
	private
	public matcrs, matcrs_part, operator(*), assignment(=)
	public spmatvec, crs2dense, read_matcrs_array, init_matcrs, part_matcrs
	public print_matdense, print_matcrs, print_matcrs_array, print_matcrs_2d, print_matcrs_part_array
	public write_matcrs_part_array, read_file_matcrs_part_array, order_matcrs
contains
	!�a�s��x�N�g����
	function spmatvec(mat,vec) result(sp)
		double precision, allocatable :: sp(:)
		type(matcrs), intent(in) :: mat
		double precision, intent(in) :: vec(:)
		integer :: i, j
		allocate(sp(mat%n))
		sp = 0
		!$omp do private(j)
		do i=1, mat%n
			do j=mat%row(i-1)+1, mat%row(i)
				sp(i) = sp(i) + mat%e(j)*vec(mat%col(j))
			end do
		end do
		!$omp end do
	end function
	
	!CRS�`�����疧�s��ɕϊ�
	subroutine crs2dense(dense,crs)
		double precision, allocatable, intent(out) :: dense(:,:)
		type(matcrs), intent(in) :: crs
		integer :: i, j
		allocate(dense(crs%n,crs%m))
		dense = 0
		do i=1, crs%n
			do j=crs%row(i-1)+1, crs%row(i)
				dense(i,crs%col(j)) = crs%e(j)
			end do
		end do
	end subroutine
	
	!���s����o��
	subroutine print_matdense(mat)
		double precision, intent(in) :: mat(:,:)
		integer :: i
		do i=1, size(mat(:,1))
			print *, mat(i,:)
		end do
	end subroutine
	
	!CRS�`���̑a�s����o�� (Matrix Market coordinate format)
	subroutine print_matcrs(mat)
		type(matcrs), intent(in) :: mat
		integer :: i, j
		print *, mat%n, mat%m, mat%k
		do i=1, mat%n
			do j=mat%row(i-1)+1, mat%row(i)
				print *, i, mat%col(j), mat%e(j)
			end do
		end do
	end subroutine
	
	!CRS�`���̑a�s����o�� (original format)
	subroutine print_matcrs_array(mat)
		type(matcrs), intent(in) :: mat
		integer :: i, j
		print *, mat%n, mat%m, mat%k
		print *, mat%e
		print *, mat%row
		print *, mat%col
	end subroutine
	
	!CRS�`���̑a�s��𖧍s��`���ŏo��
	subroutine print_matcrs_2d(mat)
		type(matcrs), intent(in) :: mat
		double precision, allocatable :: dense(:,:)
		dense = mat
		call print_matdense(dense)
	end subroutine
	
	!CRS�`���̑a�s���ǂݍ��� (original format)
	function read_matcrs_array() result(mat)
		type(matcrs) :: mat
		read *, mat%n, mat%m, mat%k
		call init_matcrs(mat)
		read *, mat%e, mat%row, mat%col
	end function
	
	function read_file_matcrs_part_array(n) result(mat)
		type(matcrs_part) :: mat
		integer, intent(in) :: n
		read(n,*) mat%mat%n, mat%mat%m, mat%mat%k
		call init_matcrs(mat%mat)
		read(n,*) mat%mat%e, mat%mat%row, mat%mat%col
		read(n,*) mat%inn, mat%ext
		allocate(mat%map(mat%inn + mat%ext))
		read(n,*) mat%map
	end function
	
	subroutine init_matcrs(mat)
		type(matcrs), intent(inout) :: mat
		allocate(mat%e(mat%k), mat%row(0:mat%n), mat%col(mat%k))
		mat%row(0) = 0
	end subroutine
	
	!CRS�`���̑a�s��𕪊�
	function part_matcrs(mat, part, npp) result(smat)
		type(matcrs_part) :: smat
		type(matcrs), intent(in) :: mat
		integer, intent(in) :: part(npp), npp
		
		integer :: i, j, inv(mat%n)
		logical :: interior(mat%n), exterior(mat%n)
		
		smat%mat%n = npp
		smat%mat%k = 0
		do i=1, npp
			smat%mat%k = smat%mat%k + mat%row(part(i)) - mat%row(part(i)-1)
		end do
		call init_matcrs(smat%mat)
		
		do i=1, smat%mat%n
			smat%mat%row(i) = mat%row(part(i)) - mat%row(part(i)-1) + smat%mat%row(i-1)
			smat%mat%e(smat%mat%row(i-1)+1:smat%mat%row(i)) = mat%e(mat%row(part(i)-1)+1:mat%row(part(i)))
			smat%mat%col(smat%mat%row(i-1)+1:smat%mat%row(i)) = mat%col(mat%row(part(i)-1)+1:mat%row(part(i)))
		end do
		smat%mat%m = maxval(smat%mat%col)
		
		interior = .false.
		exterior = .false.
		do i=1, smat%mat%n
			do j=smat%mat%row(i-1)+1, smat%mat%row(i)
				if(has(part, smat%mat%col(j))) then
					interior(smat%mat%col(j)) = .true.
				else
					exterior(smat%mat%col(j)) = .true.
				end if
			end do
		end do
		
		smat%inn = count(interior)
		smat%ext = count(exterior)
		allocate(smat%map(smat%inn + smat%ext))
		
		inv = 0
		do i=1, smat%inn
			j = indexof(interior, .true.)
			smat%map(i) = j
			inv(j) = i
			interior(j) = .false.
		end do
		
		do i=1, smat%ext
			j = indexof(exterior, .true.)
			smat%map(i+smat%inn) = j
			inv(j) = i+smat%inn
			exterior(j) = .false.
		end do
		
		do i=1, smat%mat%k
			smat%mat%col(i) = inv(smat%mat%col(i))
		end do
		
	end function
	
	!
	subroutine print_matcrs_part_array(mat)
		type(matcrs_part) :: mat
		call print_matcrs_array(mat%mat)
		print *, mat%inn, mat%ext
		print *, mat%map
	end subroutine
	
	subroutine write_matcrs_part_array(n, mat)
		integer :: n
		type(matcrs_part) :: mat
		write(n,*) mat%mat%n, mat%mat%m, mat%mat%k
		write(n,*) mat%mat%e
		write(n,*) mat%mat%row
		write(n,*) mat%mat%col
		write(n,*) mat%inn, mat%ext
		write(n,*) mat%map
	end subroutine
	
	function order_matcrs(mat, order) result(omat)
		type(matcrs) :: omat
		type(matcrs), intent(in) :: mat
		integer, intent(in) :: order(mat%n)
		integer :: i, j, osi, oei, si, ei, inv(mat%n)
		
		do i=1, mat%n
			inv(i) = indexof(order, i)
		end do
		omat%n = mat%n
		omat%m = mat%m
		omat%k = mat%k
		call init_matcrs(omat)
		do i=1, mat%n
			si = mat%row(order(i)-1) + 1
			ei = mat%row(order(i))
			osi = omat%row(i-1) + 1
			omat%row(i) = ei - si + osi
			oei = omat%row(i)
			omat%e(osi:oei) = mat%e(si:ei)
			do j=osi, oei
				omat%col(j) = inv(mat%col(si+j-osi))
			end do
		end do
	end function
	
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
	
	integer function indexof_b(a, x) result(idx)
		logical, intent(in) :: a(:), x
		integer :: i
		idx = -1
		do i=1, size(a)
			if(a(i) .eqv. x) then
				idx = i
				exit
			end if
		end do
	end function

	integer function indexof_i(a, x) result(idx)
		integer, intent(in) :: a(:), x
		integer :: i
		idx = -1
		do i=1, size(a)
			if(a(i) == x) then
				idx = i
				exit
			end if
		end do
	end function
end module

