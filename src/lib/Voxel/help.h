#include <gprim.hpp>


namespace netgen
{
  
double Det (const Mat<3,3> & m) 
{
  return 
    m(0,0) * m(1,1) * m(2,2)
    + m(1,0) * m(2,1) * m(0,2)
    + m(2,0) * m(0,1) * m(1,2)
    - m(0,0) * m(2,1) * m(1,2)
    - m(1,0) * m(0,1) * m(2,2)
    - m(2,0) * m(1,1) * m(0,2);
}

void CalcInverse (const Mat<3,3> & m, Mat<3,3> & inv)
{
  double det = Det (m);
  if (det == 0) 
    {
      inv = 0;
      return;
    }

  double idet = 1.0 / det;
  inv(0,0) =  idet * (m(1,1) * m(2,2) - m(1,2) * m(2,1));
  inv(1,0) = -idet * (m(1,0) * m(2,2) - m(1,2) * m(2,0));
  inv(2,0) =  idet * (m(1,0) * m(2,1) - m(1,1) * m(2,0));

  inv(0,1) = -idet * (m(0,1) * m(2,2) - m(0,2) * m(2,1));
  inv(1,1) =  idet * (m(0,0) * m(2,2) - m(0,2) * m(2,0));
  inv(2,1) = -idet * (m(0,0) * m(2,1) - m(0,1) * m(2,0));

  inv(0,2) =  idet * (m(0,1) * m(1,2) - m(0,2) * m(1,1));
  inv(1,2) = -idet * (m(0,0) * m(1,2) - m(0,2) * m(1,0));
  inv(2,2) =  idet * (m(0,0) * m(1,1) - m(0,1) * m(1,0));
}
}

