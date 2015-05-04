/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 *   Contributing authors:
 *   Seyed Mehdi Vaez Allaei (University of Tehran) smvaez at gmail.com
 *   Morteza Jalalvand (University of Tehran) jalalvand.m at gmail.com
 * ---------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(addcapillaryforce,FixAddCapillaryForce)

#else

#ifndef LMP_FIX_ADDCAPILLARYFORCE_H
#define LMP_FIX_ADDCAPILLARYFORCE_H

#include "fix.h"
#include <list>

namespace LAMMPS_NS {

class FixAddCapillaryForce : public Fix {
 public:
  FixAddCapillaryForce(class LAMMPS *, int, char **);
  virtual ~FixAddCapillaryForce();
  int setmask();
  void init();
  void setup(int);
  virtual void pre_force(int);
  virtual void pre_force_respa(int, int, int);

  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 protected:
   int startvar, stopvar, slopevar, interceptvar;
   char *startstr, *stopstr, *slopestr, *interceptstr;
   int startstyle, stopstyle, slopestyle, interceptstyle;
   int varflag;

  int nlevels_respa,current_nmax,nlocal_neigh,maxbridge;
  double start_trigger, stop_trigger;
  double slope, intercept;

  double **trigger;

  std::list <int> *bridgelist;

  void variable_update();
  void bridge_update();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for fix addcapillayforce does not exist

Self-explanatory.

E: Variable for fix addcapillayforce is invalid style

Self-explanatory.

*/
