#ifndef POISSON_VTK_HH
#define POISSON_VTK_HH
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <vector>

namespace PoissonVtk {
// ASCII Float64 retains the exported coefficient values. The caller supplies
// simplex cells in global-DOF numbering; only rank zero writes the final mesh.
template<class VisitCells>
void write(std::filesystem::path const& file, int dimension,
           std::vector<std::array<double,3>> const& points,
           std::vector<double> const& values, std::size_t count, VisitCells visitCells) {
  if ((dimension!=2 && dimension!=3) || points.empty() || points.size()!=values.size() || !count)
    throw std::invalid_argument("Invalid Poisson VTK geometry");
  for(std::size_t i=0;i<points.size();++i) {
    if(!std::isfinite(values[i])) throw std::invalid_argument("Nonfinite VTK solution");
    for(double x:points[i]) if(!std::isfinite(x)) throw std::invalid_argument("Nonfinite VTK coordinate");
  }
  if(count>std::size_t(std::numeric_limits<std::int64_t>::max())/(dimension+1))
    throw std::invalid_argument("VTK connectivity exceeds Int64");
  auto temporary=file; temporary += ".tmp";
  if(!file.parent_path().empty()) std::filesystem::create_directories(file.parent_path());
  std::ofstream out;out.exceptions(std::ios::failbit|std::ios::badbit);out.open(temporary);
  out<<std::setprecision(17);
  out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n";
  out<<"<Piece NumberOfPoints=\""<<points.size()<<"\" NumberOfCells=\""<<count<<"\">\n";
  out<<"<PointData Scalars=\"u\">\n<DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">\n";
  for(double v:values) out<<v<<'\n';
  out<<"</DataArray>\n<DataArray type=\"Int64\" Name=\"global_dof\" format=\"ascii\">\n";
  for(std::size_t i=0;i<points.size();++i) out<<i<<'\n';
  out<<"</DataArray>\n</PointData>\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(auto const& p:points) out<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';
  out<<"</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  std::size_t seen=0;
  visitCells([&](auto const& ids,int) {
    if(ids.size()!=std::size_t(dimension+1)) throw std::invalid_argument("VTK cell must be a simplex");
    for(auto id:ids) {
      if(std::size_t(id)>=points.size()) throw std::invalid_argument("VTK vertex index out of bounds");
      out<<id<<' ';
    }
    out<<'\n';++seen;
  });
  if(seen!=count) throw std::invalid_argument("VTK cell count mismatch");
  out<<"</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for(std::size_t i=1;i<=count;++i) out<<i*(dimension+1)<<'\n';
  out<<"</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for(std::size_t i=0;i<count;++i) out<<(dimension==3 ? 10 : 5)<<'\n';
  out<<"</DataArray>\n</Cells>\n<CellData Scalars=\"subdomain\">\n<DataArray type=\"Int32\" Name=\"subdomain\" format=\"ascii\">\n";
  visitCells([&](auto const&,int subdomain) {out<<subdomain<<'\n';});
  out<<"</DataArray>\n</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  out.close();std::filesystem::rename(temporary,file);
}
}
#endif
