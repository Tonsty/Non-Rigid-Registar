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
    for (int j = 0; j < 10; ++j)
    {
      Coord_Diff sample;
      sample.x = cos(i/5 * PI);
      sample.y = sin(j/5 * PI);
      sample.dx = 10;
      sample.dy = 10;
      samples->push_back(sample);
    }
  }

  std::cout << "Start!" << std::endl;

  TPS_Morpher morpher(samples, 1.0);
  std::vector<Point> pts;
  Point pt;
  pt.x = 1.0;
  pt.y = 2.0;
  pts.push_back(pt);
  morpher.morph(pts);
  printf("%f %f\n", pts[0].x, pts[0].y);

  return 0;
}