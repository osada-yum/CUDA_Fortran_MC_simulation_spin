module output_utilities_m
  use, intrinsic :: iso_fortran_env
  use variance_kahan_m
  use variance_covariance_kahan_m
  implicit none
contains
  subroutine output_parameters_from_disorder(nall, mcs, order_parameter, order_parameter_y, autocorrelation)
    integer(int64), intent(in) :: nall
    integer(int32), intent(in) :: mcs
    type(variance_kahan), intent(in) :: autocorrelation(mcs)
    type(variance_covariance_kahan), intent(in) :: order_parameter(mcs), order_parameter_y(mcs)
    integer(int32) :: i
    write(output_unit, '(a)') "# N, Nsample, time, <m>, <e>, <m^2>, <e^2>, N*Var[mx], N*Var[e], N*Cov[mx,e], <A>, <A^2>, N*Var[A], <m_y>"
    do i = 1, mcs
       write(output_unit, '(*(g0, 1x))') nall, order_parameter(i)%num_sample(), i, &
            & order_parameter(i)%mean1(), order_parameter(i)%mean2(), &
            & order_parameter(i)%square_mean1(), order_parameter(i)%square_mean2(), &
            & nall * order_parameter(i)%var1(), nall * order_parameter(i)%var2(), &
            & nall * order_parameter(i)%cov(), &
            & autocorrelation(i)%mean(), &
            & autocorrelation(i)%square_mean(), &
            & nall * autocorrelation(i)%var(), &
            & order_parameter_y(i)%mean1()
    end do
  end subroutine output_parameters_from_disorder

  subroutine output_abs_parameters_from_disorder(nall, mcs, order_parameter_abs, order_parameter_x, order_parameter_y, autocorrelation)
    integer(int64), intent(in) :: nall
    integer(int32), intent(in) :: mcs
    type(variance_covariance_kahan), intent(in) :: order_parameter_abs(mcs)
    type(variance_kahan), intent(in) :: order_parameter_x(mcs), order_parameter_y(mcs), autocorrelation(mcs)
    integer(int32) :: i
    write(output_unit, '(a)') "# N, Nsample, time, <|m|>, <e>, <m^2>, <e^2>, N*Var[|m|], N*Var[e], N*Cov[|x|,e], <A>, <A^2>, N*Var[A]"&
         & //", <m_x>, <m_x^2>, N*Var[mx], <m_y>, <m_y^2>, N*Var[my]"
    do i = 1, mcs
       write(output_unit, '(*(g0, 1x))') nall, order_parameter_abs(i)%num_sample(), i, &
            & order_parameter_abs(i)%mean1(), order_parameter_abs(i)%mean2(), &
            & order_parameter_abs(i)%square_mean1(), order_parameter_abs(i)%square_mean2(), &
            & nall * order_parameter_abs(i)%var1(), nall * order_parameter_abs(i)%var2(), &
            & nall * order_parameter_abs(i)%cov(), &
            & autocorrelation(i)%mean(), autocorrelation(i)%square_mean(), &
            & nall * autocorrelation(i)%var(), &
            & order_parameter_x(i)%mean(), order_parameter_x(i)%square_mean(), &
            & nall * order_parameter_x(i)%var(), &
            & order_parameter_y(i)%mean(), order_parameter_y(i)%square_mean(), &
            & nall * order_parameter_y(i)%var()
    end do
  end subroutine output_abs_parameters_from_disorder
end module output_utilities_m
