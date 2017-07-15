#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include <stdlib.h>

void interpolate(double* weights_x, double* weights_y, size_t N,
                 double* caveman, double* predictions, size_t Npred)
{
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* spline_steffen = gsl_spline_alloc(gsl_interp_steffen, N);

  gsl_spline_init(spline_steffen, weights_x, weights_y, N);

  int i;

  for (i = 0; i < Npred; i++)
  {
    predictions[i] = gsl_spline_eval(spline_steffen, caveman[i], acc);
  }

  gsl_spline_free(spline_steffen);
  gsl_interp_accel_free(acc);
}

double interpolateSingleValue(double* weights_x, double* weights_y, size_t N,
                              double interval)
{
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* spline_steffen = gsl_spline_alloc(gsl_interp_steffen, N);

  gsl_spline_init(spline_steffen, weights_x, weights_y, N);

  double threshold = gsl_spline_eval(spline_steffen, interval, acc);

  gsl_spline_free(spline_steffen);
  gsl_interp_accel_free(acc);

  return threshold;
}
