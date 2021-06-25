#ifndef __L2_DISTANCE_H__
#define __L2_DISTANCE_H__

#include <utilities/types.h>
#include <utilities/log.h>

#include <vector>
#include <fstream>
#include <functional>
#include <cmath>
#include <string>
#include <sstream>


typedef     std::vector<double>                                     Point;
typedef     std::vector<Point>                                      PointContainer;

struct L2Distance:
    public std::binary_function<const Point&, const Point&, double>
{
    result_type     operator()(const Point& p1, const Point& p2) const
    {
        AssertMsg(p1.size() == p2.size(), "Points must be in the same dimension (in L2Distance): dim1=%d, dim2=%d", p1.size(), p2.size());
        result_type sum = 0;
        for (size_t i = 0; i < p1.size(); ++i)
            sum += (p1[i] - p2[i])*(p1[i] - p2[i]);

        return sqrt(sum);
    }
};

void    read_points(const std::string& infilename, PointContainer& points)
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
