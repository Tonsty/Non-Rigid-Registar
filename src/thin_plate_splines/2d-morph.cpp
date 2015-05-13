#include <math.h>
#include "linalg3d.h"
#include "2d-morph.h"

int main(int argc, char**argv)
{
  std::cout << "Hello World!" << std::endl;

  std::vector<Coord_Diff> *p = new std::vector<Coord_Diff>;

  std::auto_ptr< std::vector<Coord_Diff> > samples(p);

  for (int i = 0; i < 10; ++i)
  {
	  Coord_Diff sample;
	  sample.x = 0.5*cos(1.0f * i/5 * PI);
	  sample.y = 0.5*sin(1.0f * i/5 * PI);
	  sample.dx = 0.2;
	  sample.dy = 0.2;
	  samples->push_back(sample);
  }

  std::cout << "Start!" << std::endl;

  TPS_Morpher morpher(samples, 0.5);
  std::vector<Point> pts;
  Point pt;
  pt.x = 0.3;
  pt.y = 0.3;
  pts.push_back(pt);
  morpher.morph(pts);
  printf("%f %f\n", pts[0].x, pts[0].y);

  return 0;
}