#ifndef __L2_DISTANCE_H__
#define __L2_DISTANCE_H__

#include <utilities/types.h>

#include <vector>
#include <fstream>
#include <functional>
#include <cmath>
#include <string>
#include <sstream>


typedef     std::vector<double>                                     Point;
typedef     std::vector<Point>                                      PointContainer;

struct WeightedL2Distance:
    public std::binary_function<const Point&, const Point&, double>
{
    result_type     operator()(const Point& p1, const Point& p2) const
    {
        AssertMsg(p1.size() == p2.size(), "Points must be in the same dimension (in L2Distance): dim1=%d, dim2=%d", p1.size()-1, p2.size()-1);

        /* the distance of a point to itself is the radius at which the "power distance" ball around it appears:
           d(p,p) := sqrt(-w_p) */
        if (p1 == p2)
            return sqrt(-p1[p1.size()-1]); 

        /* otherwise, the distance is the radius at which the power balls around p, q intersect:
           d(p,q) := sqrt( ( (w_p-w_q)^2 / 4||p-q||^2 ) - (w_p+w_q) / 2 + ||p-q||^2 / 4 ) */
        result_type sq_l2dist = 0;
        result_type wp = p1[p1.size()-1];
        result_type wq = p2[p2.size()-1];
        for (size_t i = 0; i < p1.size()-1; ++i)
            sq_l2dist += (p1[i] - p2[i])*(p1[i] - p2[i]);
        return sqrt( (wp-wq)*(wp-wq) / (4 * sq_l2dist) - (wp+wq) / 2 + sq_l2dist / 4 );
    }
};

void    read_weighted_points(const std::string& infilename, PointContainer& points)
{
    std::ifstream in(infilename.c_str());
    std::string   line;
    while(std::getline(in, line))
    {
        if (line[0] == '#') continue;               // comment line in the file
        std::stringstream linestream(line);
        double x;
        points.push_back(Point());
        while (linestream >> x)
            points.back().push_back(x);
    }
}

#endif // __L2_DISTANCE_H__
