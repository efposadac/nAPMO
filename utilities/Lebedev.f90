module lebedev
    use lebedev_1
    use lebedev_2
    use lebedev_3
    use lebedev_4
    use lebedev_5
    implicit none
    !>
    !! @brief Handles Lebedev routines
    !! @author E. F. Posada (eposada@sissa.it)
    !! @param x, y, z: Cartesian coordinates
    !! @param w: weights
    !! @see V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59, No. 3, 477 (1999)

contains 

    subroutine Lebedev_compute(nlebedev, lorder, t, p, w)
        !References:
        !
        !1. V.I. Lebedev, and D.N. Laikov, **A quadrature formula for the sphere of the 131st algebraic order of accuracy**, Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
        !2. V.I. Lebedev, **A quadrature formula for the sphere of 59th algebraic order of accuracy**,  Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286.
        !3. V.I. Lebedev, and A.L. Skorokhodov, **Quadrature formulas of orders 41, 47, and 53 for the sphere**, Russian Acad.Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592.
        !4. V.I. Lebedev, **Spherical quadrature formulas exact to orders 25-29**, Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.
        !5. V.I. Lebedev, **Quadratures on a sphere**, Computational Mathematics and Mathematical Physics, Vol. 16, 1976, pp. 10-24.
        !6. V.I. Lebedev, **Values of the nodes and weights of ninth to seventeenth order Gauss-Markov quadrature formula invariant under the octahedron group with inversion**, Computational Mathematics and Mathematical Physics, Vol. 15, 1975, pp. 44-51.
        !
        !This subroutine takes the integers ``nlebedev`` and ``lorder`` as input and returns ``lorder`` number of 
        !abscissas and weights of the Lebedev quadrature in the array LQ. The output is as follows:
        !   
        !``lq[:,0] = phi``
        !
        !``lq[:,1] = theta``
        !
        !``lq[:,2] = weight``
        !
        !Here the point, :math:`(x,y,x)`, on the surface of the sphere is given by:
        !
        !:math:`x = sin \theta * cos \phi`
        !
        !:math:`y = sin \theta * sin \phi`
        !
        !:math:`z = cos \theta`
        !
        implicit none
        integer :: lorder !Number of angular points
        integer :: nlebedev !Order or Lebedev quadrature.
        real(8), dimension(lorder) :: t, p, w !Spherical grid

        if ( nlebedev .ge. 1  .and. nlebedev .le. 5   ) call lebedev1(nlebedev,lorder, t, p, w)
        if ( nlebedev .ge. 6  .and. nlebedev .le. 10  ) call lebedev2(nlebedev,lorder, t, p, w)
        if ( nlebedev .ge. 11 .and. nlebedev .le. 15  ) call lebedev3(nlebedev,lorder, t, p, w)
        if ( nlebedev .ge. 16 .and. nlebedev .le. 20  ) call lebedev4(nlebedev,lorder, t, p, w)
        if ( nlebedev .ge. 21 .and. nlebedev .le. 24  ) call lebedev5(nlebedev,lorder, t, p, w)

    end subroutine Lebedev_compute

end module lebedev
