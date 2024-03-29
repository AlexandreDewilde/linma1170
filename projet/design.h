#pragma once
#ifndef DESIGN_H // le header guard
#define DESIGN_H 

void designTuningFork(double const r1, double const r2, double const e, double const l, double const meshSizeFactor, char const* filename);

void designTuningForkSymmetric(double const r1, double const r2, double const e, double const l, double const meshSizeFactor, char const* filename);

void designTuningForkNLayer(double const lh, double const wh, double const* d, double const* dec, double const* h, double const* l, size_t const n, double meshSizeFactor, char const* filename);

void designTuningForkSymmetricNLayer(double const lh, double const wh, double const* d, double const* dec, double const* h, double const* l, size_t const n, double const meshSizeFactor, char const* filename);

#endif