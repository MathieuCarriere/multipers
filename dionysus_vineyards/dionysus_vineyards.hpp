#include <utilities/log.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <topology/lsvineyard.h>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <time.h>

namespace bl = boost::lambda;


typedef     double                                      VertexValue;
typedef     unsigned                                    Vertex;
typedef     std::vector<VertexValue>                    VertexVector;
struct SubscriptFunctor: public std::unary_function<Vertex, VertexValue>
{
                                SubscriptFunctor(const VertexVector& v): vec(&v)    {}
        float                   operator()(Vertex i) const                          { return (*vec)[i]; }
        SubscriptFunctor&       operator=(const SubscriptFunctor& other)            { vec = other.vec; return *this; }
        const VertexVector*     vec;
};
typedef     SubscriptFunctor                            VertexEvaluator;
typedef     std::vector<VertexVector>                   VertexVectorVector;
typedef     LSVineyard<Vertex, VertexEvaluator>         PLVineyard;
typedef     PLVineyard::Simplex                         Smplx;              // gotta start using namespaces

std::vector<std::vector<std::vector<double>>> vineyards(const std::vector<std::vector<double> >& vertices_values, const std::string& complex_fn, const int& discard_inf, const int& trajectories){

  clock_t start, end;

  // Read in the complex
  PLVineyard::LSFiltration simplices;
  std::ifstream   in(complex_fn.c_str());
  std::string     line;
  while (std::getline(in, line)){
    std::istringstream  strin(line);
    simplices.push_back(Smplx(std::istream_iterator<Vertex>(strin), std::istream_iterator<Vertex>()));
  }
  //std::cout << "Simplices read:" << std::endl;
  //std::copy(simplices.begin(), simplices.end(), std::ostream_iterator<Smplx>(std::cout, "\n"));

  // Read in vertex values
  VertexVectorVector vertices;
  int n = vertices_values.size();
  for (int i = 0; i < n; i++){
    vertices.push_back(VertexVector(vertices_values[i].begin(), vertices_values[i].end()));
  }
  //std::cout << "Vertex values read:" << std::endl;
  //for (size_t i = 0; i < vertices.size(); ++i)
  //{
  //      std::copy(vertices[i].begin(), vertices[i].end(), std::ostream_iterator<VertexValue>(std::cout, " "));
  //      std::cout << std::endl;
  //}

  // Setup the vineyard
  start = clock();
  VertexEvaluator veval(vertices[0]);
  PLVineyard::VertexComparison vcmp(veval);
  PLVineyard::SimplexComparison scmp(vcmp);
  simplices.sort(scmp);
  PLVineyard v(boost::counting_iterator<Vertex>(0), boost::counting_iterator<Vertex>(vertices[0].size()), simplices, veval);
  end = clock(); 
  std::cout << double(end-start)/CLOCKS_PER_SEC << std::endl;

  // Compute vineyard
  for (size_t i = 1; i < vertices.size(); ++i){
    veval = VertexEvaluator(vertices[i]);
    v.compute_vineyard(veval);
  }

  // Retrieve vineyard
  std::vector<std::vector<std::vector<double>>> V;
  if (trajectories)  V = v.vineyard().get_vines(discard_inf);
  else  V = v.vineyard().get_dgms(discard_inf, n);

  return V;

}
