#include "math.h"

void Math::Simpson_Integral
(
    const int mesh,
    const double *func,
    const double *rab,
    double &asum
)
{
    /*     simpson's rule integration. On input:
    !      mesh = mhe number of grid points (should be odd)
    !      func(i)= function to be integrated
    !      rab(i) = r(i) * dr(i)/di * di
    !      For the logarithmic grid not including r=0 :
    !      r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
    !      For the logarithmic grid including r=0 :
    !      r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
    !      Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
    !      where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
    */
    //  simpson's rule integrator for function stored on the
    //  radial logarithmic mesh
    //	routine assumes that mesh is an odd number so run check
    if (mesh % 2 == 0)
    {
        cout << "\n error in subroutine simpson ";
        cout << "\n routine assumes mesh is odd but mesh = "
             << mesh << endl;
        return;
    }

    asum = 0.00;
    const double r12 = 1.00 / 12.00;
    double f3 = func [0] * rab [0] * r12;
    for (int i = 1;i < mesh;i += 2)
    {
        const double f1 = f3;
        const double f2 = func [i] * rab [i] * r12;
        f3 = func [i + 1] * rab [i + 1] * r12;
        asum += 4.00 * f1 + 16.00 * f2 + 4.00 * f3;
    }
    return;
}// end subroutine simpson
