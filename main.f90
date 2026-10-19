program read_spectra_v2
    
    implicit none

    integer, parameter :: n1=525, n2=573, n3=645, n4=741, n5=861
    real(kind=8) :: S1(n1), S2(n2), S3(n3), S4(n4), S5(n5)
    real(kind=8) :: w1(n1), w2(n2), w3(n3), w4(n4), w5(n5)
    real(kind=8) :: x1(n1), x2(n2), x3(n3), x4(n4), x5(n5)
    real(kind=8) :: S1_theor(n1), S2_theor(n2), S3_theor(n3), S4_theor(n4), S5_theor(n5)
    real(kind=8) :: weight1(n1)=1, weight2(n2)=1, weight3(n3)=1, weight4(n4)=1, weight5(n5)=1
    integer :: unit, ios
    real(kind=8) :: PI = 3.1415926535897

    call get_data_from_file(S1, 1, n1)
    call get_data_from_file(S2, 2, n2)
    call get_data_from_file(S3, 3, n3)
    call get_data_from_file(S4, 4, n4)
    call get_data_from_file(S5, 5, n5)

    ! get spectrum fit for each spectrum and save to text file

    ! SPECTRUM 1
    ! for weight = 1
    call compute_spectrum_fit(S1_theor, x1, w1, S1, weight1, n1, 1)
    call save_data(x1, w1, S1_theor, S1, n1, "S1.txt")

    ! for W=F(exp)**2
    weight1 = S1**2
    call compute_spectrum_fit(S1_theor, x1, w1, S1, weight1, n1, 1)
    call save_data(x1, w1, S1_theor, S1, n1, "S1_weighted.txt")

    ! SPECTRUM 2
    weight2 = S2**2
    call compute_spectrum_fit(S2_theor, x2, w2, S2, weight2, n2, 2)
    call save_data(x2, w2, S2_theor, S2, n2, "S2_weighted.txt")

    ! SPECTRUM 3
    weight3 = S3**2
    call compute_spectrum_fit(S3_theor, x3, w3, S3, weight3, n3, 3)
    call save_data(x3, w3, S3_theor, S3, n3, "S3_weighted.txt")

    ! SPECTRUM 4
    weight4 = S4**2
    call compute_spectrum_fit(S4_theor, x4, w4, S4, weight4, n4, 4)
    call save_data(x4, w4, S4_theor, S4, n4, "S4_weighted.txt")

    ! SPECTRUM 5
    weight5 = S5**2
    call compute_spectrum_fit(S5_theor, x5, w5, S5, weight5, n5, 5)
    call save_data(x5, w5, S5_theor, S5, n5, "S5_weighted.txt")
    

contains

    subroutine save_data(x, w, spectrum_theor, spectrum_exp, num_points, save_filename)
        ! save data to text file

        integer, intent(in) :: num_points
        character(len=*), intent(in) :: save_filename
        real(kind=8), intent(in) :: x(num_points), w(num_points), spectrum_theor(num_points), spectrum_exp(num_points)
        integer :: k

        ! create file and save values to text file
        open(unit=10, file=save_filename, status="replace")
        do k = 1, num_points
            write(10, *) x(k), w(k), spectrum_theor(k), spectrum_exp(k)
        end do
        close(10)

    end subroutine save_data

    subroutine get_data_from_file(spectrum, file_index, num_points)
        ! get data from text file 

        integer, intent(in) :: file_index, num_points
        real(kind=8), intent(out) :: spectrum(num_points)
        character(len=30) :: filename
        integer :: j
        real :: val

        j=1

        ! make file name
        write(filename, '(A,I1,A)') 'spectre_', file_index, '.txt'

        ! open file
        open(newunit=unit, file=filename, iostat=ios)

        ! check file opened with no error

        if (ios /= 0) then
            print *, 'cannot open ', trim(filename)
            stop
        end if

        ! read all lines of files and populate spectrum with data
        do
            read(unit, *, iostat=ios) val
            if (ios < 0) exit ! end of file
            
            spectrum(j) = val
            j = j+1
        end do

        close(unit)
        print *, trim(filename), ': ', j-1, ' values read'

    end subroutine get_data_from_file

    subroutine compute_spectrum_fit(spectrum_theor, x, w, spectrum_exp, weights, num_points, index)
        ! compute theoretical spectrum uing least squares fitting

        integer, intent(in) :: num_points, index
        real(kind=8), intent(in) :: weights(num_points), spectrum_exp(num_points)
        real(kind=8), intent(out) :: x(num_points), w(num_points), spectrum_theor(num_points)
        real(kind=8) :: sigma, w_mean, a0, a1, a2, spectrum_inverse(num_points), wm, gamma, big_S
        integer j

        print *, "Spectrum fit for spectra ", index

        w = 0

        ! make list of w values
        do j=1, num_points
            w(j) = 2280.0 + (j-1)*0.01
        end do

        ! calculate mean and rms of w
        call calc_mean_and_rms(sigma, w_mean, w, num_points)
        print *, "sigma", sigma, "w mean", w_mean

        ! make lists of x values
        x = (w - w_mean) / sigma
        
        ! take inverse of F
        spectrum_inverse = 1/spectrum_exp

        ! calculate coefficients
        call calc_coefficients(a0, a1, a2, x, spectrum_inverse, num_points, weights)

        ! calculate spectra parameters
        wm = w_mean - sigma*a1/(2*a2)
        gamma = sigma*SQRT( a0/a2 - (a1**2)/ (4*(a2**2)) )
        big_S = PI*sigma / SQRT( a0*a2 - (a1**2)/4 )

        ! display parameters
        print *, "wm", wm, "gamma", gamma, "s", big_S

        ! make theoretical data set for S
        spectrum_theor = a2*x**2 + a1*x + a0
        spectrum_theor = 1/spectrum_theor

    end subroutine compute_spectrum_fit

    subroutine calc_mean_and_rms(rms, mean, values, num)
        ! calculate mean and rms of list of values

        real(kind=8), intent(out) :: rms, mean
        integer, intent(in) :: num 
        real(kind=8), intent(in) :: values(num)
        real(kind=8) :: mean_sqr

        ! mean: mean of values
        ! mean_sqr: mean of square of values

        mean = SUM(values) / num
        mean_sqr = SUM(values**2) / num

        rms = SQRT(mean_sqr - mean**2)

    end subroutine calc_mean_and_rms

    subroutine calc_coefficients(a0, a1, a2, x, Y, n, weight)
        ! calculate quadratic coefficients for a set of weights
        
        real(kind=8), intent(out) :: a0,a1,a2
        integer, intent(in) :: n
        real(kind=8), intent(in) :: x(n), Y(n), weight(n)
        real(kind=8) :: avg_1=0, avg_x=0, avg_x2=0, avg_x3=0, avg_x4=0, avg_xY=0, avg_x2Y=0, avg_Y=0
    
        ! calculate weighted means
        avg_1 = SUM(weight)/n
        avg_x = SUM(weight*x)/n
        avg_x2 = SUM(weight*x**2)/n
        avg_x3 = SUM(weight*x**3)/n
        avg_x4 = SUM(weight*x**4)/n
        avg_xY = SUM(weight*x*Y)/n
        avg_x2Y = SUM(weight*x**2*Y)/n
        avg_Y = SUM(weight*Y)/n

        ! print *, "avg_1", avg_1, "avg_x", avg_x, "avg_x2", avg_x2, "avg_x3", avg_x3, "avg_x4", avg_x4, "avg_xy", avg_xy, "avg_x2y", avg_x2y, "avg_y", avg_y

        ! calculate coefficients using determinant equations
        a0 = ( avg_y*avg_x2*avg_x4 + avg_xY*avg_x3*avg_x2 + avg_x*avg_x3*avg_x2Y - avg_x2Y*avg_x2**2 - avg_y*avg_x3**2 - avg_x*avg_xy*avg_x4) / &
        ( avg_1*avg_x2*avg_x4 + 2*avg_x*avg_x2*avg_x3 - avg_x2**3 - avg_1*avg_x3**2 - avg_x**2*avg_x4 )

        a1 = ( avg_1*avg_xy*avg_x4 + avg_x2*avg_x3*avg_y + avg_x*avg_x2*avg_x2y - avg_x2**2*avg_xy - avg_x*avg_x4*avg_y - avg_1*avg_x3*avg_x2y ) / &
        ( avg_1*avg_x2*avg_x4 + 2*avg_x*avg_x2*avg_x3 - avg_x2**3 - avg_1*avg_x3**2 - avg_x**2*avg_x4 )

        a2 = ( avg_1*avg_x2*avg_x2y + avg_x*avg_x3*avg_y + avg_x*avg_x2*avg_xy - avg_x2**2*avg_y - avg_x**2*avg_x2y - avg_1*avg_x3*avg_xy ) / &
        ( avg_1*avg_x2*avg_x4 + 2*avg_x*avg_x2*avg_x3 - avg_x2**3 - avg_1*avg_x3**2 - avg_x**2*avg_x4 )

        ! display coefficients
        print *, "a0:", a0, "a1:", a1, "a2:", a2

    end subroutine calc_coefficients

end program read_spectra_v2