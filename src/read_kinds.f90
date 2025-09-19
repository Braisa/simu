program read_kinds
    use get_kinds

    implicit none

    integer (int) :: i
    real (double) :: d

    print *, 'integer: ', int
    print *, 'double: ', double

    print *, 'integer bits: ', STORAGE_SIZE(i, int)
    print *, 'double bits: ', STORAGE_SIZE(d, double)
end program read_kinds