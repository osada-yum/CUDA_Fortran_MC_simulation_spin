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

  subroutine output_abs_parameters_from_disorder(nall, mcs, order_parameter_abs, order_parameter_xy, autocorrelation)
    integer(int64), intent(in) :: nall
    integer(int32), intent(in) :: mcs
    type(variance_covariance_kahan), intent(in) :: order_parameter_abs(mcs), order_parameter_xy(mcs)
    type(variance_kahan), intent(in) :: autocorrelation(mcs)
    real(real64) :: chi
    integer(int32) :: i
    write(output_unit, '(a)') "# N, Nsample, time" // &
         & ", <|m|>, <e>" // &
         & ", <m^2>, <e^2>" //&
         & ", <|m|e>, (<m^2> - (<mx>^2 + <my>^2))" // &
         & ", <A>, <A^2>" // &
         & ", <mx>, <my>" // &
         & ", <mx^2>, <my^2>, <mx*my>"
    do i = 1, mcs
       chi = order_parameter_abs(i)%square_mean1() - (order_parameter_xy(i)%mean1() ** 2 + order_parameter_xy(i)%mean2() ** 2)
       write(output_unit, '(*(g0, 1x))') nall, order_parameter_abs(i)%num_sample(), i, &
            & order_parameter_abs(i)%mean1(), order_parameter_abs(i)%mean2(), &
            & order_parameter_abs(i)%square_mean1(), order_parameter_abs(i)%square_mean2(), &
            & order_parameter_abs(i)%mean_v1v2(), chi, &
            & autocorrelation(i)%mean(), autocorrelation(i)%square_mean(), &
            & order_parameter_xy(i)%mean1(), order_parameter_xy(i)%mean2(), &
            & order_parameter_xy(i)%square_mean1(), order_parameter_xy(i)%square_mean2(), order_parameter_xy(i)%mean_v1v2()
    end do
  end subroutine output_abs_parameters_from_disorder
end module output_utilities_m
