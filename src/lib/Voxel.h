#ifndef VOXEL_H
#define VOXEL_H

// c++ -O2 -fopenmp voxel.cpp help.cpp adtree.cpp optmem.cpp -I../../libsrc/include
// c++ -O2 -fopenmp voxel.cpp help.cpp adtree.cpp optmem.cpp -I/home/parallels/ngsolve/sources/netgen-5.1/libsrc/include -o voxel
// tar -czf voxel.tar.gz *.cpp *.pqr


#include <gprim.hpp>
#include <lib/Voxel/adtree.h>
#include <lib/Voxel/help.h>
#include <lib/Voxel/optmem.h>

using namespace netgen;


class Voxel {

public:

	Voxel(){
		rball = 1.4;
	};
	~Voxel(){};

int calcSurface(mMesh &mSurface, vector<Atom> &atomList, INI& ini, string fileName)
{
    int gridSize = atoi(ini.get<string>("model.grid_resolution").c_str());
	rball = atof(ini.get<string>("experiment.probe_radius").c_str());

  int nx = gridSize;   // number of intervals
  int ny = gridSize;
  int nz = gridSize;

  Array<Point<3> > gridpoints(long(nx+1)*(ny+1)*(nz+1));
  Array<float> values(long(nx+1)*(ny+1)*(nz+1));
  Array<Vec<3> > dvalues(long(nx+1)*(ny+1)*(nz+1));

  clock_t t1 = clock();
  rmax = 0;

  for (unsigned int i = 0; i < atomList.size(); i++)
    {
	  Point3<float> coord = atomList.at(i).getCoord();
	  double r = atomList.at(i).getRadius();
	  pnts.Append (Point<3>(coord.X(),coord.Y(),coord.Z()));
	  rad.Append (r);
	  if (r > rmax) rmax = r;
    }

  cout << "rmax = " << rmax << endl;
  int increase = ceil(rmax)+1;

  clock_t t2 = clock();

  netgen::Box<3> box;
  box.Set (pnts[0]);
  for (int i = 1; i < pnts.Size(); i++) box.Add (pnts[i]);
  box.Increase (increase);
  Point<3> pmin = box.PMin();
  Point<3> pmax = box.PMax();


  searchtree = new Point3dTree (pmin, pmax);
  for (int i = 0; i < pnts.Size(); i++)
    searchtree -> Insert (pnts[i], i);


  for (int ix = 0; ix <= nx; ix++)
    for (int iy = 0; iy <= ny; iy++)
      for (int iz = 0; iz <= nz; iz++)
        {
          int ind = ix + (nx+1) * (iy + (ny+1)*iz);
          gridpoints[ind] = Point<3> (pmin(0) + (pmax(0)-pmin(0))*ix/nx,
                                      pmin(1) + (pmax(1)-pmin(1))*iy/ny,
                                      pmin(2) + (pmax(2)-pmin(2))*iz/nz);
        }


  SetValues (pmin, pmax, gridpoints, values, dvalues);

  clock_t t3 = clock();

  double h = Dist(pmin, pmax)/nx;
  for (int i = 0; i < values.Size(); i++)
    if (fabs (values[i]) < 0.02 * dvalues[i].Length() * h)
      {
        gridpoints[i] -= (values[i]/dvalues[i].Length2()) * dvalues[i];
        values[i] = 0;
      }

  stlout.open (fileName.c_str());
  stlout << "solid" << endl;
  stlout.precision(12);
  MakeSTL (pmin, pmax, nx, ny, nz, gridpoints, values, dvalues);
  stlout << "endsolid" << endl;

  clock_t t4 = clock();

  cout << "load: " << double(t2-t1) / CLOCKS_PER_SEC << " secs" << endl;
  cout << "find f: " << double(t3-t2) / CLOCKS_PER_SEC << " secs" << endl;
  cout << "write stl (" << fileName <<"): " << double(t4-t3) / CLOCKS_PER_SEC << " secs" << endl;

  stlout.close();
  return 0;
}

private:

	double rball;
	double rmax;
	ofstream stlout;

	Array<Point<3> > pnts;
	Array<double> rad;
	Point3dTree * searchtree;

void SetValues (Point<3> pmin, Point<3> pmax,
                // int nx, int ny, int nz,
                Array<Point<3> > & gridpoints,
                Array<float> & values,
                Array<Vec<3> >  & dvalues)
{
  int cnt = 0;

  int cnt3a = 0, cnt3b = 0, cnt3c = 0;

  cout << "Setting up level-set function" << endl;

#pragma omp parallel
  {
    Array<int> indices1, indices;


#pragma omp for
    for (int i = 0; i < gridpoints.Size(); i++)
      {

#pragma omp atomic
        cnt++;


        if (cnt % 1000 == 0)
          {
#pragma omp critical(print)
            {
              cout << "\r" << int (100 * double(cnt)/gridpoints.Size()) << " %" << flush;
            }
          }
        /*
          int ind = ix + (nx+1) * (iy + (ny+1)*iz);
          Point<3> p(pmin(0) + (pmax(0)-pmin(0))*ix/nx,
          pmin(1) + (pmax(1)-pmin(1))*iy/ny,
          pmin(2) + (pmax(2)-pmin(2))*iz/nz);
        */
        int ind = i;
        int extend = 6;
        if (rmax > 5)
        	extend = 2*rmax;

        Point<3> p = gridpoints[ind];

        searchtree -> GetIntersecting (p - Vec<3>(extend,extend,extend), p+Vec<3>(extend,extend,extend), indices1);
        double minf = 1e10;
        Vec<3> df;


        indices.SetSize(0);
        for (int j = 0; j < indices1.Size(); j++)
          if (Dist2 (pnts[indices1[j]], p) < sqr(rad[indices1[j]]+2*rball))
            indices.Append (indices1[j]);

        for (int jj = 0; jj < indices.Size(); jj++)
          {
            int j = indices[jj];
            // double fj = rad[j] * (Dist2 (p, pnts[j]) / sqr(rad[j]) - 1);
            double fj = Dist (p, pnts[j]) - rad[j];
            if (fj < minf)
              {
                minf = fj;
                df = 1.0/Dist(p,pnts[j]) * (p-pnts[j]);
              }
          }

        for (int jj = 0; jj < indices.Size(); jj++)
          for (int kk = 0; kk < jj; kk++)
            {
              int j = indices[jj];
              int k = indices[kk];

              double r1 = rad[j];
              double r2 = rad[k];

              Point<3> c1 = pnts[j];
              Point<3> c2 = pnts[k];

              Vec<3> ex = c2-c1;
              double dist2 = ex.Length2();
              if (dist2 > sqr(r1+r2+2*rball)) continue;
              double dist = sqrt(dist2);
              ex /= dist;

              double px = ex * (p - c1);
              if (px <= -rball || px > dist+rball) continue;

              double py = sqrt ( Dist2(p,c1) - px*px );

              double tx = 1.0 / (2*dist) * (sqr(dist)+(r1-r2)*(2*rball+r1+r2));
              double ty = sqrt (sqr(r1+rball) - tx*tx);

              // if ( (py / px < ty / tx) &&  (py / (dist-px) < ty / (dist-tx)))
              if ( (py * tx < ty * px) && ( (dist-tx)*py < ty * (dist-px) ) )
                {
                  // double fj = rball * (1-(sqr(px-tx)+sqr(py-ty)) / sqr(rball));
                  double fj = rball - sqrt(sqr(px-tx)+sqr(py-ty));

                  double x0 = tx - ty/(ty-py) * (tx-px);
                  Point<3> p0 = c1+x0*ex;
                  Vec<3> dir = p - p0;
                  dir.Normalize();
                  if (fj < minf)
                    {
                      minf = fj;
                      df = dir;
                    }
                }
            }

        for (int jj = 0; jj < indices.Size(); jj++)
          for (int kk = 0; kk < jj; kk++)
            for (int ll = 0; ll < kk; ll++)
              {
#pragma omp atomic
                cnt3a++;

                int j = indices[jj];
                int k = indices[kk];
                int l = indices[ll];

                double r1 = rad[j];
                double r2 = rad[k];
                double r3 = rad[l];

                Point<3> c1 = pnts[j];
                Point<3> c2 = pnts[k];
                Point<3> c3 = pnts[l];

                if (Dist2(c1, p) > sqr(r1+2*rball)) continue;
                if (Dist2(c2, p) > sqr(r2+2*rball)) continue;
                if (Dist2(c3, p) > sqr(r3+2*rball)) continue;


                Vec<3> t1 = c2-c1;
                Vec<3> t2 = c3-c1;
                Vec<3> n = Cross(t1, t2);
                n.Normalize();

                Mat<2> mat;
                Vec<2> rhs, sol;
                mat(0,0) = t1*t1;
                mat(0,1) = mat(1,0) = t1*t2;
                mat(1,1) = t2*t2;

                rhs(0) = 0.5 * (sqr(r1+rball) - sqr(r2+rball) + Dist2(c1,c2));
                rhs(1) = 0.5 * (sqr(r1+rball) - sqr(r3+rball) + Dist2(c1,c3));

                mat.Solve (rhs, sol);

                Point<3> cp = c1 + sol(0) * t1 + sol(1) * t2;
                if ( sqr(r1+rball) <= Dist2(c1,cp) ) continue;

                double lamn = sqrt ( sqr(r1+rball) - Dist2(c1,cp) );
                // c = cp +/- lamn n

                Vec<2> rhs2, sol2;
                rhs2(0) = t1 * (p-c1);
                rhs2(1) = t2 * (p-c1);
                mat.Solve (rhs2, sol2);
                double lamn2 = n * (p-c1);

#pragma omp atomic
                    cnt3b++;


                if ( sol2(0) > sol(0) * fabs(lamn2)/lamn &&
                     sol2(1) > sol(1) * fabs(lamn2)/lamn &&
                     (1-sol2(0)-sol2(1)) > (1-sol(0)-sol(1)) * fabs(lamn2)/lamn)

                  {
#pragma omp atomic
                    cnt3c++;


                    Point<3> sp1 = cp + lamn * n;
                    Point<3> sp2 = cp - lamn * n;


                    /*
                    double fj = max (rball - Dist(sp1, p), rball - Dist(sp2, p));
                    if (fj < minf)
                      {
                        minf = fj;

                        Vec<3> dir;
                        if (lamn2 > 0)
                          dir = sp1-p;
                        else
                          dir = sp2-p;
                        dir.Normalize();
                        df = dir;
                      }
                    */
                    if (lamn2 > 0)
                      {
                        Vec<3> dir = sp1-p;
                        double len = dir.Length();
                        double fj = rball - len;
                        if (fj < minf)
                          {
                            dir.Normalize();
                            minf = fj;
                            df = dir;
                          }
                      }
                    else
                      {
                        Vec<3> dir = sp2-p;
                        double len = dir.Length();
                        double fj = rball - len;
                        if (fj < minf)
                          {
                            dir.Normalize();
                            minf = fj;
                            df = dir;
                          }
                      }

                  }
              }





        values[ind] = minf;
        dvalues[ind] = df;
      }
  }

  cout << "\r" << "100 %" << endl;
  // cout << "cnt trip-balls: " << cnt3a << " " << cnt3b << " " << cnt3c << endl;
}



void WriteTrig (Point<3> p1, Point<3> p2, Point<3> p3, Vec<3> n)
{

  if (Dist (p1, p2) < 1e-10) return;
  if (Dist (p1, p3) < 1e-10) return;
  if (Dist (p2, p3) < 1e-10) return;

  Vec<3> nt = Cross(p2-p1,p3-p1);
  if (nt.Length() < 1e-5 * (p2-p1).Length() * (p3-p1).Length())
    {
      if (nt.Length() < 1e-12 * (p2-p1).Length() * (p3-p1).Length()) return;
     // cout << "flat trig: " << nt.Length() / ((p2-p1).Length() * (p3-p1).Length()) << endl;
    }
  if (nt * n < 0) Swap (p2, p3);

  stlout << "facet normal " << n(0) << " " << n(1) << " " << n(2) << "\n";
  stlout << "  outer loop" << "\n";
  stlout << "    vertex " << p1(0) << " " << p1(1) << " " << p1(2) << "\n";
  stlout << "    vertex " << p2(0) << " " << p2(1) << " " << p2(2) << "\n";
  stlout << "    vertex " << p3(0) << " " << p3(1) << " " << p3(2) << "\n";
  stlout << "  endloop" << "\n";
  stlout << "endfacet" << "\n";
}

double CutEdge (double f1, double f2)
{
  if (fabs (f1) < 1e-10) return 0;
  return f1/(f1-f2);
}

void MakeTetSTL (Point<3> pnts[], double valtet[4], Vec<3> dvaltet[4])
{
  int npos = 0;
  for (int j = 0; j < 4; j++)
    if (valtet[j] > 0) npos++;

  if (npos == 0 || npos == 4) return;


  const int edges[6][2] =
    { { 0, 1 }, { 0, 2 }, { 0, 3 },
      { 2, 3 }, { 1, 3 }, { 1, 2 } };

//  double lame[6] = { -1 };
  Point<3> cutp[4];
  int cutpi = 0;
  for (int j = 0; j < 6; j++)
    {
      int pi1 = edges[j][0];
      int pi2 = edges[j][1];
      if ((valtet[pi1] > 0) != (valtet[pi2] > 0))
        {
//          Vec<3> ve = pnts[pi2]-pnts[pi1];
          double lam = CutEdge (valtet[pi1], valtet[pi2]);
          cutp[cutpi] = pnts[pi1] + lam * (pnts[pi2]-pnts[pi1]);
          cutpi++;
        }
    }


  // calc grad f (normal vector)
  Mat<3,3> mat;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      mat(i,j) = pnts[i](j) - pnts[3](j);
  Vec<3> rhs, sol;
  for (int i = 0; i < 3; i++)
    rhs(i) = valtet[i] - valtet[3];
  mat.Solve (rhs, sol);

  if (npos == 1 || npos == 3)   // cut 3 edges
    WriteTrig (cutp[0], cutp[1], cutp[2], sol);

  if (npos == 2)   // cut 4 edges
    {
      WriteTrig (cutp[0], cutp[1], cutp[2], sol);
      WriteTrig (cutp[0], cutp[2], cutp[3], sol);
    }
}


void MakeCubeSTL (//Point<3> p1, Point<3> p2,
                  Point<3> pointcube[8],
                  double valcube[8], Vec<3> dvalcube[8])
{
  const int tets[6][4] =
    { { 0, 1, 2, 4 },
      { 1, 2, 4, 5 },
      { 2, 4, 5, 6 },
      { 1, 3, 2, 5 },
      { 3, 2, 5, 6 },
      { 3, 6, 5, 7 } };

  double valtet[4];
  Vec<3> dvaltet[4];
  Point<3> pnts[4];

  for (int j = 0; j < 6; j++)
    {
      for (int k = 0; k < 4; k++)
        {
          valtet[k] = valcube[tets[j][k]];
          dvaltet[k] = dvalcube[tets[j][k]];
          pnts[k] = pointcube[tets[j][k]];
        }

      MakeTetSTL (pnts, valtet, dvaltet);
    }
}



void MakeSTL (Point<3> pmin, Point<3> pmax,
              int nx, int ny, int nz,
              Array<Point<3> > & gridpoints,
              Array<float> & values,
              Array<Vec<3> > & dvalues)
{
  double valcube[8];
  Vec<3> dvalcube[8];
  Point<3> pointcube[8];

  cout << "Writing stl file" << endl;
  for (int ix = 0; ix < nx; ix++)
    {
      cout << "\r" << ix << "/" << nx << flush;
      for (int iy = 0; iy < ny; iy++)
        for (int iz = 0; iz < nz; iz++)
          {
            for (int jx = 0, jj = 0; jx < 2; jx++)
              for (int jy = 0; jy < 2; jy++)
                for (int jz = 0; jz < 2; jz++, jj++)
                  {
                    int ind = ix+jx + (nx+1) * (iy+jy + (ny+1)*(iz+jz));
                    valcube[jj] = values[ind];
                    dvalcube[jj] = dvalues[ind];
                    pointcube[jj] = gridpoints[ind];
                  }

            int npos = 0;
            for (int j = 0; j < 8; j++)
              if (valcube[j] > 0) npos++;
            if (npos == 0 || npos == 8) continue;

            MakeCubeSTL (pointcube, valcube, dvalcube);
          }
    }
  cout << "\r" << nx << "/" << nx << endl;
}


};



#endif
