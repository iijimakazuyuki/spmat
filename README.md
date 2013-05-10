spmat
=====

Fortran sparse matrix module

types
-----

###matcrs

Sparse matrix in CRS format

####n

Num of raws

####m

Num of elements

####e

Elements

####idx

Indices of each raw

####col

Columns of each element

###matcrs_part

Part of sparse matrix

####mat

Type of matcrs

####inn

Num of internal points

####neib

Num of external points

####map

Mapping of part indices to general indices

arithmetic
----------

####spmatvec

Sparse matrix-vector multiplication

####crs2dense

Convert CRS format to dense matrix (two-dimensional array)

io
--

####read_matcrs_array

Read matcrs data (original format: n m e idx col)

####print_matdense

Print a two-dimensional array as a dense matrix

####print_matcrs

Print matcrs data (Matrix Market coordinate format)

####print_matcrs_array

Print matcrs data (cf. read_matcrs_array)
Can be read using read_matcrs_array

####print_matcrs_2d

Print matcrs data as a dense matrix

####print_matcrs_part_array

Print matcrs_part data (original format: n m e idx col inn neib map)

####write_matcrs_part_array

Write matcrs_part data (cf. print_matcrs_part)

####read_file_matcrs_part_array

Read matcrs_part data from a file

operation
---------

####init_matcrs

Utility subroutine
Allocate dynamic arrays of matcrs

####part_matcrs

Part CRS matrix

####order_matcrs

Order CRS matrix
