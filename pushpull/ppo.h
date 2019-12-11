#pragma once

// todo: second signature that takes an existing pattern for initialization

// Points will be written to outMatrix
void optimizePattern(double dMin, double rC, double areaDeltaMax, unsigned int nPoints, int initType, double * outMatrix);