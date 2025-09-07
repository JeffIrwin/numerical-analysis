
!> Unit tests for chapter 6: eigenvalue problems
module numa__chapter_6

	use numa
	use numa__exercises

	implicit none

contains

!===============================================================================

integer function chapter_6_basic_qr() result(nfail)

	character(len = *), parameter :: label = "chapter_6_basic_qr"

	double precision :: diff
	double precision, allocatable :: a(:,:), d(:,:), s(:,:), eigvals(:), &
		expect(:)

	integer :: n, irep

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Fuzz test

	call rand_seed_determ()

	! There are fewer tests here and with smaller matrices compared to
	! qr_factor() because eig_basic_qr() is an expensive iterative algorithm,
	! involving a full QR decomposition and matmul at each iteration
	do n = 4, 15, 3

		allocate(s (n, n))

		!if (mod(n, 10) == 0) then
		!	print *, "Testing basic QR algorithm with n = " // to_str(n) // " ..."
		!end if

		do irep = 1, 1

			! Construct a random matrix `a` with known real eigenvalues

			! Known eigenvalues
			expect = zeros(n)
			call random_number(expect)

			d = diag(expect)
			call sort(expect)
			!print *, "expect  = ", expect
			print "(a,*(es15.5))", " expect  = ...", expect(n-3: n)

			!call random_number(s)  ! random matrix
			s = rand_f64(n, n)

			a = matmul(matmul(s, d), inv(s))
			!print *, "a = "
			!print "(4es15.5)", a

			eigvals = eig_basic_qr(a, iters = 5 * n)
			print "(a,*(es15.5))", " eigvals = ...", eigvals(n-3: n)

			! There is a large tolerance here because the basic QR algorithm
			! converges slowly
			diff = norm2(eigvals - expect)
			call test(diff, 0.d0, 1.d-2 * n, nfail, "eig_basic_qr 1")
			print *, "diff = ", diff
			print *, ""

		end do

		deallocate(s)
	end do

	!********
	!print *, ""

end function chapter_6_basic_qr

!===============================================================================

integer function chapter_6_hessenberg_qr() result(nfail)

	character(len = *), parameter :: label = "chapter_6_hessenberg_qr"

	double precision :: diff, t, t0, t_kr, t_pq
	double precision, allocatable :: a(:,:), d(:,:), s(:,:), eigvals(:), &
		expect(:), eigvecs(:,:), a0(:,:)

	integer :: n, irep, iters

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	!********
	! Fuzz test

	call rand_seed_determ()

	t_pq = 0.d0
	t_kr = 0.d0

	do n = 5, 25, 7

		allocate(s (n, n))

		print *, "Testing Hessenberg with n = " // to_str(n) // " ..."

		do irep = 1, 1

			! Construct a random matrix `a` with known real eigenvalues

			! Known eigenvalues
			expect = zeros(n)
			call random_number(expect)

			!expect(3) = expect(2)

			d = diag(expect)
			call sort(expect)
			!print *, "expect  = ", expect
			print "(a,*(es18.8))", " expect  = ...", expect(n-3: n)

			!call random_number(s)  ! random matrix
			s = rand_f64(n, n)
			!call print_mat(s, "s = ")

			a = matmul(matmul(s, d), inv(s))
			!print *, "a = "
			!print "(4es18.8)", a

			a0 = a

			! Benchmark results:  null space (kernel) algorithm is nearly twice
			! as fast as the Q product algorithm for iters = 5 * n**2:
			!
			!     Q product  fuzz time = 2.6070019999999996 s
			!     Null space fuzz time = 1.5295850000000000 s
			!
			! Q product costs more per iteration, because it has to update the
			! product.  Null space is faster per iteration, but it has a more
			! expensive step at the end to calculate the eigenvectors.  Null
			! space has to perform n QR decompositions to get the eigenvectors,
			! while the Q product method only has to do one LU solve on the full
			! n x n matrix (but many more LU solves on smaller matrices of size
			! n-1, n-2, n-3, ...)
			!
			! If you only do a very small number of iterations, Q product can
			! come out faster, but it doesn't converge very well.  For practical
			! numbers of iterations, the null space algo seems better
			!
			! Also, the `diff vec` printed below has a 10x smaller residual for
			! the null space algo
			!
			! Shifting algorithms should converge faster, so it's good to keep
			! the Q product algorithm for now.  It may be the better performer
			! with shifting

			iters = 5 * n**2
			!iters = n**2
			!iters = n
			!iters = 10*n

			!********
			print *, "Hessenberg QR with Q product:"

			! Hessenberg QR doesn't converge in any fewer iterations than basic
			! QR, but each iteration is only O(n**2) instead of O(n**3), so it's
			! cheap to just do a bunch of iterations
			call cpu_time(t0)
			eigvals = eig_hess_qr(a, iters, eigvecs)
			call cpu_time(t)
			t_pq = t_pq + t - t0

			!print *, "eigvecs = "
			!print "(5es15.5)", eigvecs

			! Check the eigenvalues
			diff = norm2(sorted(eigvals) - expect)
			call test(diff, 0.d0, 1.d-4 * n, nfail, "eig_hess_qr val 1")
			print *, "diff val = ", diff

			!print *, "a * eigvecs / eigvecs = "
			!print "(5es15.5)", matmul(a0, eigvecs) / eigvecs
			!!print *, "spread = "
			!!print "(4es15.5)", spread(eigvals, 1, n)

			! Check the eigenvectors: A * eigvecs = eigvals * eigvecs
			diff = norm2(matmul(a0, eigvecs) / eigvecs - spread(eigvals, 1, n))
			call test(diff, 0.d0, 1.d-3 * n, nfail, "eig_hess_qr vec 2")
			print *, "diff vec = ", diff

			print *, ""
			!********
			print *, "Hessenberg QR with kernel:"

			a = a0
			call cpu_time(t0)
			eigvals = eig_hess_qr_kernel(a, iters, eigvecs)
			call cpu_time(t)
			t_kr = t_kr + t - t0

			! The magnitudes of eigvecs may differ from above, but their
			! directions are the same, thus the spread test passes in all cases
			!print *, "eigvecs = "
			!print "(5es15.5)", eigvecs

			! Check the eigenvalues
			diff = norm2(sorted(eigvals) - expect)
			call test(diff, 0.d0, 1.d-4 * n, nfail, "eig_hess_qr_kernel val 1")
			print *, "diff val = ", diff

			!print *, "a * eigvecs / eigvecs = "
			!print "(5es15.5)", matmul(a0, eigvecs) / eigvecs
			!!print *, "spread = "
			!!print "(4es15.5)", spread(eigvals, 1, n)

			! Check the eigenvectors: A * eigvecs = eigvals * eigvecs
			diff = norm2(matmul(a0, eigvecs) / eigvecs - spread(eigvals, 1, n))
			call test(diff, 0.d0, 1.d-3 * n, nfail, "eig_hess_qr_kernel vec 2")
			print *, "diff vec = ", diff

			print *, ""

			!********
			!stop

		end do

		deallocate(s)
	end do

	write(*,*) "Q product  fuzz time = ", to_str(t_pq), " s"
	write(*,*) "Null space fuzz time = ", to_str(t_kr), " s"

	!********
	print *, ""

end function chapter_6_hessenberg_qr

!===============================================================================

integer function chapter_6_francis_qr() result(nfail)

	character(len = *), parameter :: label = "chapter_6_francis_qr"

	double precision :: diff
	double precision, allocatable :: a(:,:), d(:,:), s(:,:), &
		expect(:), a0(:,:)!, ar(:,:), ai(:,:)
	double complex, allocatable :: eigvals(:), eigvals2(:), eigvecs(:,:)

	integer :: n, irep, p0

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	! Matrix from literature with complex eigenvalues

	n = 6
	allocate(a(n,n))
	a(:,1) = [ 7,  3,  4, -11, -9, -2]
	a(:,2) = [-6,  4, -5,   7,  1, 12]
	a(:,3) = [-1, -9,  2,   2,  9,  1]
	a(:,4) = [-8,  0, -1,   5,  0,  8]
	a(:,5) = [-4,  3, -5,   7,  2, 10]
	a(:,6) = [ 6,  1,  4, -11, -7, -1]
	print *, "a = "
	print "("//to_str(n)//"es19.9)", a
	a0 = a
	!    --> spec(a')
	!     ans  =
	!
	!       5. + 6.i
	!       5. - 6.i
	!       1. + 2.i
	!       1. - 2.i
	!       3. + 0.i
	!       4. + 0.i

	eigvals = eig_francis_qr(a, eigvecs)
	print *, "eigvals = [real, imag]"
	print "(2es15.5)", sorted(eigvals)

	! Eigenvectors are optional to save work.  Beware, rounding and ordering may be different
	a = a0
	eigvals2 = eig_francis_qr(a)
	diff = diff_complex_vecs(eigvals, eigvals2)
	print *, "diff eigvals2 = ", diff
	call test(diff, 0.d0, 1.d-8, nfail, "eig_francis_qr no vecs")

	!print *, "eigvals2 = [real, imag]"
	!print "(2es15.5)", sorted(eigvals2)

	!print *, "eigvecs = "
	!print "("//to_str(2*n)//"es15.5)", eigvecs

	!print *, "a * eigvecs / eigvecs = "
	!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs) / eigvecs

	!print *, "a * eigvecs = "
	!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs)

	!print *, "eigvals * eigvecs = "
	!print "("//to_str(2*n)//"es15.5)", spread(eigvals, 1, n) * eigvecs

	!print *, "a * eigvecs - eigvals * eigvecs = "
	!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs

	diff = norm2(abs( &
		matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
	print *, "diff = ", diff
	call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_francis_qr vec 1")

	!! Some components of some eigenvectors are 0, so dividing like this is
	!! ill-conditioned.  Need to multiply instead
	!diff = norm2(abs(matmul(a0, eigvecs) / eigvecs - spread(eigvals, 1, n)))
	!call test(diff, 0.d0, 1.d-3 * n, nfail, "eig_francis_qr vec 2")

	!********
	! Fuzz test

	call rand_seed_determ()

	p0 = -1
	do n = 5, 35, 3

		allocate(s (n, n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing real eig_francis_qr() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 2

			! Construct a random matrix `a` with known real eigenvalues

			! Known eigenvalues
			expect = zeros(n)
			call random_number(expect)

			!expect(3) = expect(2)

			d = diag(expect)
			call sort(expect)
			!print "(a,*(es18.8))", " expect  = ...", expect(n-3: n)

			call random_number(s)  ! random matrix
			a = matmul(matmul(s, d), inv(s))
			!print *, "a = "
			!print "("//to_str(n)//"es19.9)", a

			a0 = a

			!********

			eigvals = eig_francis_qr(a, eigvecs)

			!print *, "eigvecs = "
			!print "("//to_str(2*n)//"es15.5)", eigvecs

			! Check the eigenvalues

			diff = norm2(sorted(dble(eigvals)) - expect)

			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_francis_qr val 1")
			!print *, "diff val = ", diff
			!print *, "expect = ", expect
			!print *, "actual = ", sorted(dble(eigvals))

			!********

			!print *, "a * eigvecs / eigvecs = "
			!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs) / eigvecs
			!!print *, "spread = "
			!!print "("//to_str(n)//"4es15.5)", spread(eigvals, 1, n)

			! Check the eigenvectors: A * eigvecs = eigvals * eigvecs
			!
			! idk if norm2(abs( [complex_matrix] )) is a good norm but i don't feel like making my
			! own matrix norm rn.  need to figure out what built-in norm2() is
			! doing on matrices, whether it just reshapes them to vec or does
			! some more involved matrix norm

			!print *, "eigvals = ", eigvals

			diff = norm2(abs( &
				matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_francis_qr vec 3")
			!print *, "diff vec = ", diff
			!print *, ""

		end do

		deallocate(s)
	end do
	deallocate(a)

	p0 = -1
	do n = 5, 65, 7
		! eig_francis_qr() can't converge over about n == 80

		allocate(a(n, n))

		if (n/10 > p0) then
		!if (n/1 > p0) then
			p0 = n/10
			print *, "Testing complex eig_francis_qr() with n = " // to_str(n) // " ..."
		end if

		!do irep = 1, 10
		do irep = 1, 1
			!print *, "irep = ", irep

			! Construct a random matrix `a`.  It can have complex eigenvalues,
			! most likely a mixture of real and complex eigenvalues
			call random_number(a)
			!print *, "a = "
			!print "("//to_str(n)//"es19.9)", a
			a0 = a

			!********

			eigvals = eig_francis_qr(a, eigvecs)

			!print *, "eigvals = [real, imag]"
			!print "(2es15.5)", eigvals

			!print *, "eigvecs = "
			!print "("//to_str(2*n)//"es15.5)", eigvecs

			! Eigenvalues cannot be verified by themselves because we don't know
			! what they should be a priori
			!
			! We can still verify that `a0 * eigvecs == eigvals * eigvecs`

			!print *, "a * eigvecs / eigvecs = "
			!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs) / eigvecs
			!!print *, "spread = "
			!!print "("//to_str(n)//"4es15.5)", spread(eigvals, 1, n)

			!print *, "a * eigvecs - eigvals * eigvecs = "
			!print "("//to_str(2*n)//"es15.5)", matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs

			!print *, "checking diff"
			!diff = norm2(abs(matmul(a0, eigvecs) / eigvecs - spread(eigvals, 1, n)))
			diff = norm2(abs( &
				matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_francis_qr vec 4")
			!print *, "diff vec = ", diff

			!print *, ""

		end do

		deallocate(a)
	end do
	!********
	print *, ""

end function chapter_6_francis_qr

!===============================================================================

integer function chapter_6_eig_lapack() result(nfail)

	character(len = *), parameter :: label = "chapter_6_eig_lapack"

	double precision :: diff
	double precision, allocatable :: a(:,:), d(:,:), s(:,:), &
		expect(:), a0(:,:)!, ar(:,:), ai(:,:)
	double complex, allocatable :: eigvals(:), eigvals2(:), eigvecs(:,:)

	integer :: n, irep, p0

	write(*,*) CYAN // "Starting " // label // "()" // COLOR_RESET

	nfail = 0

	! Matrix from literature with complex eigenvalues

	n = 6
	allocate(a(n,n))
	a(:,1) = [ 7,  3,  4, -11, -9, -2]
	a(:,2) = [-6,  4, -5,   7,  1, 12]
	a(:,3) = [-1, -9,  2,   2,  9,  1]
	a(:,4) = [-8,  0, -1,   5,  0,  8]
	a(:,5) = [-4,  3, -5,   7,  2, 10]
	a(:,6) = [ 6,  1,  4, -11, -7, -1]
	print *, "a = "
	print "("//to_str(n)//"es19.9)", a
	a0 = a
	!    --> spec(a')
	!     ans  =
	!
	!       5. + 6.i
	!       5. - 6.i
	!       1. + 2.i
	!       1. - 2.i
	!       3. + 0.i
	!       4. + 0.i

	eigvals = eig_lapack(a, eigvecs)
	print *, "eigvals = [real, imag]"
	print "(2es15.5)", sorted(eigvals)

	! Eigenvectors are optional to save work.  Beware, rounding and ordering may be different
	a = a0
	eigvals2 = eig_lapack(a)
	diff = diff_complex_vecs(eigvals, eigvals2)
	print *, "diff eigvals2 = ", diff
	call test(diff, 0.d0, 1.d-8, nfail, "eig_lapack no vecs")

	diff = norm2(abs( &
		matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
	print *, "diff = ", diff
	call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_lapack vec 1")

	!stop

	!! Some components of some eigenvectors are 0, so dividing like this is
	!! ill-conditioned.  Need to multiply instead
	!diff = norm2(abs(matmul(a0, eigvecs) / eigvecs - spread(eigvals, 1, n)))
	!call test(diff, 0.d0, 1.d-3 * n, nfail, "eig_lapack vec 2")

	!********
	! Fuzz test

	call rand_seed_determ()

	p0 = -1
	do n = 5, 25, 5

		allocate(s (n, n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing real eig_lapack() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 2

			! Construct a random matrix `a` with known real eigenvalues

			! Known eigenvalues
			expect = zeros(n)
			call random_number(expect)

			d = diag(expect)
			call sort(expect)

			call random_number(s)  ! random matrix
			a = matmul(matmul(s, d), inv(s))
			!print *, "a = "
			!print "("//to_str(n)//"es19.9)", a
			a0 = a

			!********

			eigvals = eig_lapack(a, eigvecs)

			! Check the eigenvalues
			diff = norm2(sorted(dble(eigvals)) - expect)
			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_lapack val 1")
			!********

			! Check the eigenvectors: A * eigvecs = eigvals * eigvecs
			!
			! idk if norm2(abs( [complex_matrix] )) is a good norm but i don't feel like making my
			! own matrix norm rn.  need to figure out what built-in norm2() is
			! doing on matrices, whether it just reshapes them to vec or does
			! some more involved matrix norm
			diff = norm2(abs( &
				matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_lapack vec 3")

		end do

		deallocate(s)
	end do
	deallocate(a)

	p0 = -1
	!do n = 5, 175, 5  ! eig_lapack() is robust even for larger matrices
	do n = 5, 95, 11

		allocate(a(n, n))

		if (n/10 > p0) then
			p0 = n/10
			print *, "Testing complex eig_lapack() with n = " // to_str(n) // " ..."
		end if

		do irep = 1, 1
		!do irep = 1, 2

			!print *, "irep = ", irep

			! Construct a random matrix `a`.  It can have complex eigenvalues,
			! most likely a mixture of real and complex eigenvalues
			call random_number(a)
			!print *, "a = "
			!print "("//to_str(n)//"es19.9)", a
			a0 = a
			!********

			eigvals = eig_lapack(a, eigvecs)

			! Eigenvalues cannot be verified by themselves because we don't know
			! what they should be a priori
			!
			! We can still verify that `a0 * eigvecs == eigvals * eigvecs`
			diff = norm2(abs( &
				matmul(a0, eigvecs) - spread(eigvals, 1, n) * eigvecs))
			call test(diff, 0.d0, 1.d-6 * n, nfail, "eig_lapack vec 4")

		end do

		deallocate(a)
	end do
	!********
	print *, ""

end function chapter_6_eig_lapack

!===============================================================================

end module numa__chapter_6

