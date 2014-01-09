/*=========================================================================

 Program:   ParaView
 Module:    $RCSfile XRAGEAnalysisAdaptor.h,v $

 Copyright (c) Kitware, Inc.
 All rights reserved.
 See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

    This software is distributed WITHOUT ANY WARRANTY; without even
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
    PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef MPASAdaptor_h
#define MPASAdaptor_h

#include <map>
#include <string>
#include <vector>

using namespace std;

//////////////////////////////////////////////////////////////////////////
//
// Namespace data consists of pointers to data within Fortran simulator
// Names of variables were kept the same
// Used a namespace so that all methods could refer to it
//
// Primal and dual mesh are opposite each other so vertices in one are
// the centers in the other so the same methods can be used to build
// just by passing different parameters
//
// MPASAdaptor supports three kinds of meshes and two kinds of data
// 2D and 3D Cartesian really are (x,y) grid with depth in the z dimension
// Longitude/Latitude Cartesian is (lon,lat) grid with depth in the z dimension
// Both of these contain 3D cells
// Spherical is 3D (x,y,z) with a single depth and contains 2D cells
//
//////////////////////////////////////////////////////////////////////////

namespace MPAS
{
  extern int rank;
  extern int totalRank;

  // Above parameters assigned for primal mesh
  extern int nPrimalCells;
  extern int nPrimalVerts;
  extern int nPrimalVertsPerCell;
  extern int nPrimalGhosts;
  extern int* primalGhostCell;
  extern int* primalGhostHalo;
  extern int totalPrimalCells;
  extern int* nEdgesOnCell;            // Primal cell number of sides
  extern int* verticesOnCell;          // Point indices for dual cells
  extern int* vertexMask;              // Valid values

  // Above parameters assigned for dual mesh
  extern int nDualCells;
  extern int nDualVerts;
  extern int nDualVertsPerCell;
  extern int nDualGhosts;
  extern int* dualGhostCell;           // Ghost cell indices
  extern int* dualGhostHalo;           // Ghost cell levels
  extern int totalDualCells;
  extern int* cellsOnVertex;           // Point indices for primal cells
  extern int* cellMask;                // Valid values

  extern int nVertLevels;              // Number of vertex depth levels

  extern const int PRIMAL;
  extern const int DUAL;

  extern string rankName[];
  extern string maskName[];

  // Boundary cells are not complete enough to draw and must be omitted
  extern map<string, vector<bool> > usedCells;
}

//////////////////////////////////////////////////////////////////////////
//
// Cartesian mesh creation and data load methods called from MPASAdaptor
//
//////////////////////////////////////////////////////////////////////////

void create_xy3D_grids(
                 double* xCell,
                 double* yCell,
                 double* xVert,
                 double* yVert,
                 float zFactor);

void create_xy3D_mesh(
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex,
                 double* xCenter, double* yCenter,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell,
                 double* offset,
                 double* poffset,
                 float zFactor);

//////////////////////////////////////////////////////////////////////////
//
// Spherical mesh creation and data load methods called from MPASAdaptor
//
//////////////////////////////////////////////////////////////////////////

void create_xyz2D_grids(
                 double* xCell,
                 double* yCell,
                 double* zCell,
                 double* xVert,
                 double* yVert,
                 double* zVert);

void create_xyz2D_mesh(
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex, double* zVertex,
                 double* xCenter, double* yCenter, double* zCenter,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell);

//////////////////////////////////////////////////////////////////////////
//
// Spherical mesh creation and data load methods called from MPASAdaptor
//
//////////////////////////////////////////////////////////////////////////

void create_xyz3D_grids(
                 double* xCell,
                 double* yCell,
                 double* zCell,
                 double* xVert,
                 double* yVert,
                 double* zVert,
                 float zFactor);

void create_xyz3D_mesh(
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex, double* zVertex,
                 double* xCenter, double* yCenter, double* zCenter,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell);

//////////////////////////////////////////////////////////////////////////
//
// Lon/lat mesh creation and data load methods called from MPASAdaptor
//
//////////////////////////////////////////////////////////////////////////

void create_lonlat2D_grids(
                 double* xCell,
                 double* yCell,
                 double* xVert,
                 double* yVert);

void create_lonlat2D_mesh(
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell);

//////////////////////////////////////////////////////////////////////////
//
// Lon/lat mesh creation and data load methods called from MPASAdaptor
//
//////////////////////////////////////////////////////////////////////////

void create_lonlat3D_grids(
                 double* xCell,
                 double* yCell,
                 double* xVert,
                 double* yVert,
                 float zFactor);

void create_lonlat3D_mesh(
                 int meshType,
                 int numberOfCells,
                 int numberOfVertices,
                 int verticesPerCell,
                 int numberOfGhosts,
                 int* ghostCell, int* ghostLevel,
                 double* xVertex, double* yVertex,
                 int* nEdgesOnCell,
                 int* vertices,
                 vector<bool> const& makeCell,
                 float zFactor);

#endif
