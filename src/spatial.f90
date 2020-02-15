module spatial

implicit none

type derivative_type
  integer :: order
  integer :: scheme
  integer :: stencil_width
  integer, pointer :: stencil_weights(:)
  integer, pointer :: stencil_index(:)
end type derivative_type

contains

subroutine spatial_init( deriv, order, scheme )

type(derivative_type), pointer :: pderiv 

if ( associated(pderiv) ) &
     call error('spatial_init$','Derivative is already allocated$')

allocate( pderiv )
pderiv%order = order
pderiv%scheme = scheme
pderiv%stencil_width = schemes%width(scheme)
allocate( pderiv%stencil_weights(pderiv%stencil_width) )
pderiv%stencil_weights = schemes%weights(1:schemes%width(scheme))
pderiv%stencil_index(



end subroutine spatial_init

subroutine spatial_derivative( deriv, f, df )

type(derivative_type) :: deriv


integer :: nweights
real :: weights(:)
real :: shift(:)

df = zero
do iw = 1, nweights
  df = df + weights(iw) * f(:,index+shift(iw))
end do

end subroutine spatial_derivative



end module spatial
