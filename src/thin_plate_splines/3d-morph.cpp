#include <math.h>
#include "linalg3d.h"
#include "3d-morph.h"

int main(int argc, char**argv)
{
  std::cout << "Hello World!" << std::endl;

  std::vector<Coord_Diff> *p = new std::vector<Coord_Diff>;

  std::auto_ptr< std::vector<Coord_Diff> > samples(p);

  for (int i = 0; i < 2; ++i)
  {
	  for (int j = 0; j < 2; ++j)
	  {
		  for (int k = 0; k < 2; ++k)
		  {
			  Coord_Diff sample;
			  sample.x = 0.1 * i;
			  sample.y = 0.1 * j;
			  sample.z = 0.1 * k;
			  sample.dx = 0.1;
			  sample.dy = 0.1;
			  sample.dz = 0.1;
			  samples->push_back(sample);
		  }
    }
  }

  std::cout << "Start!" << std::endl;

  TPS_Morpher morpher(samples, 1.0);
  std::vector<Point> pts;
  Point pt;
  pt.x = 0.05;
  pt.y = 0.05;
  pt.z = 0.05;
  pts.push_back(pt);
  morpher.morph(pts);
  printf("%f %f %f\n", pts[0].x, pts[0].y, pts[0].z);

  return 0;
}