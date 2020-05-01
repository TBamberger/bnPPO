#pragma once

// @param [in] initType Initial point pattern is generated according to this (0: random, 1: darts, 2: jittered grid, 3: regular grid)
// @param [out] outMatrix Resulting point pattern will be written here (MATLAB Nx2 matrix). Memory has to be preallocated.
int optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, int initType, double* outMatrix, double aspectRatio);

// @param [in] inMatrix Input point pattern for which blue noise properties will be restored. Order of point coordinates has to be according to MATLAB Nx2 matrix.
// @param [out] outMatrix Resulting point pattern will be written here (MATLAB Nx2 matrix). Memory has to be preallocated.
int optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, double* inMatrix, double* outMatrix, double aspectRatio);

// @param [in] inMatrix Input point pattern for which blue noise properties will be restored. Order of point coordinates has to be according to MATLAB Nx2 matrix.
// @param [out] outMatrix Resulting point pattern will be written here (MATLAB Nx2 matrix). Memory has to be preallocated.
int optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, double* inMatrix, double* inMatrixFixedTile, double* outMatrix, double aspectRatio);
