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

#include "math.h"
#include "string.h"
#include "fix_addcapillaryforce.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include <list>

using namespace LAMMPS_NS;
using namespace FixConst;

enum {CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixAddCapillaryForce::FixAddCapillaryForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR, "Illegal fix addcapillayforce command");

  restart_global = 1;
  restart_peratom = 1;
  create_attribute = 1;
  current_nmax = 0;

  startstr = stopstr = slopestr = interceptstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    startstr = new char[n];
    strcpy(startstr,&arg[3][2]);
  } else {
    start_trigger = force->numeric(FLERR, arg[3]);
    startstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    stopstr = new char[n];
    strcpy(stopstr,&arg[4][2]);
  } else {
    stop_trigger = force->numeric(FLERR, arg[4]);
    stopstyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    slopestr = new char[n];
    strcpy(slopestr,&arg[5][2]);
  } else {
    slope = force->numeric(FLERR, arg[5]);
    slopestyle = CONSTANT;
  }
  if (strstr(arg[6],"v_") == arg[6]) {
    int n = strlen(&arg[6][2]) + 1;
    interceptstr = new char[n];
    strcpy(interceptstr,&arg[6][2]);
  } else {
    intercept = force->numeric(FLERR, arg[6]);
    interceptstyle = CONSTANT;
  }

  atom->add_callback(0);
  atom->add_callback(1);

  if (atom->map_style == 0) {
    atom->map_init();
    atom->map_set();
  }
}

/* ---------------------------------------------------------------------- */

FixAddCapillaryForce::~FixAddCapillaryForce()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  delete [] bridgelist;
  delete [] startstr;
  delete [] stopstr;
  delete [] slopestr;
  delete [] interceptstr;
  memory->destroy(trigger);
}

/* ---------------------------------------------------------------------- */

int FixAddCapillaryForce::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddCapillaryForce::init()
{
  // check variables

  if (startstr) {
    startvar = input->variable->find(startstr);
    if (startvar < 0) error->all(FLERR, "Variable name for fix addcapillayforce"
                                       " does not exist");
    if (input->variable->equalstyle(startvar)) startstyle = EQUAL;
    else if (input->variable->atomstyle(startvar)) startstyle = ATOM;
    else error->all(FLERR, "Variable for fix addcapillayforce "
                          "is invalid style");
  }
  if (stopstr) {
    stopvar = input->variable->find(stopstr);
    if (stopvar < 0) error->all(FLERR, "Variable name for fix addcapillayforce"
                                      " does not exist");
    if (input->variable->equalstyle(stopvar)) stopstyle = EQUAL;
    else if (input->variable->atomstyle(stopvar)) stopstyle = ATOM;
    else error->all(FLERR, "Variable for fix addcapillayforce "
                          "is invalid style");
  }
  if (slopestr) {
    slopevar = input->variable->find(slopestr);
    if (slopevar < 0) error->all(FLERR, "Variable name for fix addcapillayforce"
                                       " does not exist");
    if (input->variable->equalstyle(slopevar)) slopestyle = EQUAL;
    else error->all(FLERR, "Variable for fix addcapillayforce "
                          "is invalid style");
  }
  if (interceptstr) {
    interceptvar = input->variable->find(interceptstr);
    if (interceptvar < 0) error->all(FLERR, "Variable name for fix "
                                           "addcapillayforce does not exist");
    if (input->variable->equalstyle(interceptvar)) interceptstyle = EQUAL;
    else error->all(FLERR, "Variable for fix addcapillayforce"
                          " is invalid style");
  }

  if (startstyle == ATOM || stopstyle == ATOM) varflag = ATOM;
  else if (startstyle == EQUAL || stopstyle == EQUAL || slopestyle == EQUAL ||
           interceptstyle == EQUAL) varflag = EQUAL;
  else varflag = CONSTANT;

  grow_arrays(atom->nmax);

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixAddCapillaryForce::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    pre_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    pre_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddCapillaryForce::pre_force(int vflag)
{
  int i,j,ii,inum,*mask,*ilist,*numneigh,**firstneigh;
  std::list <int>::iterator it;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,rinv;
  double capillaryforce,fx,fy,fz;

  double **x = atom->x;
  double **f = atom->f;
  double *radius = atom->radius;

  int *tag = atom->tag;
  NeighList *list = force->pair->list;
  int nlocal = atom->nlocal;
  inum = list->inum;
  ilist = list->ilist;
  mask = atom->mask;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (varflag != CONSTANT) variable_update();

  bridge_update();

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    for (it = bridgelist[i].begin(); it != bridgelist[i].end(); it++) {
      j = atom->map(*it);
      if (!(mask[j] & groupbit) || j<i) continue;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      domain->minimum_image(delx, dely, delz);
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;
      r = sqrt(rsq);
      rinv = 1/r;

      capillaryforce = (slope*(r-radsum) + intercept)*rinv;
      fx = delx*capillaryforce;
      fy = dely*capillaryforce;
      fz = delz*capillaryforce;
      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;
      if (j < nlocal) {
        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAddCapillaryForce::pre_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) pre_force(vflag);
}

/* ----------------------------------------------------------------------
 *   updates all variables
 * ------------------------------------------------------------------------- */

void FixAddCapillaryForce::variable_update()
{
  modify->clearstep_compute();

  if (startstyle == EQUAL)
    start_trigger = input->variable->compute_equal(startvar);
  if (stopstyle == EQUAL)
    stop_trigger = input->variable->compute_equal(stopvar);
  if (slopestyle == EQUAL) slope = input->variable->compute_equal(slopevar);
  if (interceptstyle == EQUAL)
    intercept = input->variable->compute_equal(interceptvar);
  if (varflag == ATOM) {
    if (startstyle == ATOM)
      input->variable->compute_atom(startvar,igroup,&trigger[0][1],2,0);
    else for (int i=0; i<current_nmax; i++) trigger[i][1] = start_trigger;
    if (stopstyle == ATOM)
      input->variable->compute_atom(stopvar,igroup,&trigger[0][0],2,0);
    else for (int i=0; i<current_nmax; i++) trigger[i][0] = stop_trigger;
  }

  modify->addstep_compute(update->ntimestep + 1);
}

/* ----------------------------------------------------------------------
 *   updates list of all capillary bridges
 * ------------------------------------------------------------------------- */

void FixAddCapillaryForce::bridge_update()
{
  int i,j,ii,jj,inum,jnum,nlocal;
  std::list <int>::iterator it;
  bool found;
  static bool previous = false;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq;
  double istart, istop, jstart, jstop;
  int *mask,*ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double *radius = atom->radius;

  int *tag = atom->tag;
  NeighList *list = force->pair->list;
  mask = atom->mask;
  nlocal = atom->nlocal;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  nlocal_neigh = (inum)?ilist[inum-1]:0;
  maxbridge = 0;

  if (varflag != ATOM) {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (!(mask[i] & groupbit)) continue;
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      radi = radius[i];

      // find and delete broken bridges of i

      for (it = bridgelist[i].begin(); it != bridgelist[i].end();) {
        j = atom->map(*it);
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        domain->minimum_image(delx, dely, delz);
        rsq = delx*delx + dely*dely + delz*delz;
        radj = radius[j];
        radsum = radi + radj + stop_trigger;

        if (rsq <= radsum*radsum) it++;
        else {
          it = bridgelist[i].erase(it);
          if (j < nlocal) bridgelist[j].remove(tag[i]);
        }
      }

      // loop over neighbors of i

      jlist = firstneigh[i];
      jnum = numneigh[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        if (!(mask[j] & groupbit)) continue;

        // skip j if it's already in the bridge list of i

        found = false;
        for (it = bridgelist[i].begin(); it != bridgelist[i].end(); it++)
          if (found = (tag[j] == *it)) break;
        if (found) continue;

        // form new bridges

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        domain->minimum_image(delx, dely, delz);
        rsq = delx*delx + dely*dely + delz*delz;
        radj = radius[j];
        radsum = radi + radj + start_trigger;

        if (rsq <= radsum*radsum) {
          bridgelist[i].push_back(tag[j]);
          if (j < nlocal) bridgelist[j].push_back(tag[i]);
        }
      }
      maxbridge = MAX (maxbridge, bridgelist[i].size());
    }
  } else {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (!(mask[i] & groupbit)) continue;
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      radi = radius[i];
      istop = trigger[i][0];
      istart = trigger[i][1];

      // find and delete broken bridges of i

      for (it = bridgelist[i].begin(); it != bridgelist[i].end();) {
        j = atom->map(*it);
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        domain->minimum_image(delx, dely, delz);
        rsq = delx*delx + dely*dely + delz*delz;
        radj = radius[j];
        jstop = trigger[j][0];
        radsum = radi + radj + istop + jstop;

        if (rsq <= radsum*radsum) it++;
        else {
          it = bridgelist[i].erase(it);
          if (j < nlocal) bridgelist[j].remove(tag[i]);
        }
      }

      // loop over neighbors of i

      jlist = firstneigh[i];
      jnum = numneigh[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        if (!(mask[j] & groupbit)) continue;

        // skip j if it's already in the bridge list of i

        found = false;
        for (it = bridgelist[i].begin(); it != bridgelist[i].end(); it++)
          if (found = (tag[j] == *it)) break;
          if (found) continue;

          // form new bridges

          delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        domain->minimum_image(delx, dely, delz);
        rsq = delx*delx + dely*dely + delz*delz;
        radj = radius[j];
        jstart = trigger[j][1];
        radsum = radi + radj + istart + jstart;

        if (rsq <= radsum*radsum) {
          bridgelist[i].push_back(tag[j]);
          if (j < nlocal) bridgelist[j].push_back(tag[i]);
        }
      }
      maxbridge = MAX (maxbridge, bridgelist[i].size());
    }
  }
}

/* ----------------------------------------------------------------------
 *   memory usage of local atom-based arrays
 * ------------------------------------------------------------------------- */

double FixAddCapillaryForce::memory_usage()
{
  double bytes = current_nmax * sizeof(std::list <int>);
  for (int i=0; i<current_nmax; i++)
    bytes += bridgelist[i].size()*(sizeof(int) + 2*sizeof(int*));
  if (varflag == ATOM) bytes += current_nmax*2*sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
 *   pack entire state of Fix into one write
 * ------------------------------------------------------------------------- */

void FixAddCapillaryForce::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    int size = 4 * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&start_trigger,sizeof(double),1,fp);
    fwrite(&stop_trigger,sizeof(double),1,fp);
    fwrite(&slope,sizeof(double),1,fp);
    fwrite(&intercept,sizeof(double),1,fp);
  }
}

/* ----------------------------------------------------------------------
 *   use state info from restart file to restart the Fix
 * ------------------------------------------------------------------------- */

void FixAddCapillaryForce::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  start_trigger = list[n++];
  stop_trigger = list[n++];
  slope = list[n++];
  intercept = list[n++];
}

/* ----------------------------------------------------------------------
 *   allocate local atom-based arrays
 * ------------------------------------------------------------------------- */

void FixAddCapillaryForce::grow_arrays(int nmax)
{
  if (current_nmax) {
    std::list <int> *tmp = new std::list <int> [nmax];
    for (int i=0; i<current_nmax; i++) tmp[i] = bridgelist[i];
    delete [] bridgelist;
    bridgelist = tmp;
  }
  else bridgelist = new std::list <int> [nmax];
  if (varflag == ATOM)
    memory->grow(trigger, nmax, 2, "addcapillayforce:trigger");
  current_nmax = nmax;
}

/* ----------------------------------------------------------------------
 *   copy values within local atom-based arrays
 * ------------------------------------------------------------------------- */

void FixAddCapillaryForce::copy_arrays(int i, int j, int delflag)
{
  std::list <int>::iterator it;
  bridgelist[j].clear();
  for (it = bridgelist[i].begin(); it != bridgelist[i].end(); it++)
    bridgelist[j].push_back(*it);
}

/* ----------------------------------------------------------------------
 *   pack values in local atom-based arrays for exchange with another proc
 * ------------------------------------------------------------------------- */

int FixAddCapillaryForce::pack_exchange(int i, double *buf)
{
  int m = 0;
  std::list <int>::iterator it;
  buf[m++] = bridgelist[i].size();
  for (it = bridgelist[i].begin(); it != bridgelist[i].end(); it++)
    buf[m++] = *it;
  bridgelist[i].clear();
  return m;
}

/* ----------------------------------------------------------------------
 *   unpack values into local atom-based arrays after exchange
 * ------------------------------------------------------------------------- */

int FixAddCapillaryForce::unpack_exchange(int i, double *buf)
{
  int j, nlocal = atom->nlocal, m = 0, npartner = static_cast <int> (buf[m++]);
  int tagj;
  maxbridge = MAX (maxbridge, npartner);
  bridgelist[i].clear();
  for (int n=0; n<npartner; n++) {
    tagj = static_cast <int> (buf[m++]);
    bridgelist[i].push_back(tagj);
  }
  return m;
}

/* ----------------------------------------------------------------------
 *   pack values in local atom-based arrays for restart file
 * ------------------------------------------------------------------------- */

int FixAddCapillaryForce::pack_restart(int i, double *buf)
{
  int m = 0;
  std::list <int>::iterator it;
  buf[m++] = bridgelist[i].size() + 1;
  for (it = bridgelist[i].begin(); it != bridgelist[i].end(); it++)
    buf[m++] = *it;
  return m;
}

/* ----------------------------------------------------------------------
 *   unpack values from atom->extra array to restart the fix
 * ------------------------------------------------------------------------- */

void FixAddCapillaryForce::unpack_restart(int i, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int n = 0; n < nth; n++) m += static_cast<int> (extra[i][m]);

  int npartner = static_cast <int> (extra[i][m++]) - 1;
  maxbridge = MAX (maxbridge, npartner);
  bridgelist[i].clear();
  for (int n=0; n<npartner; n++)
    bridgelist[i].push_back(static_cast <int> (extra[i][m++]));
}

/* ----------------------------------------------------------------------
 *   maxsize of any atom's restart data
 * ------------------------------------------------------------------------- */

int FixAddCapillaryForce::maxsize_restart()
{
  // maxbridge_all = max # of bridges across all procs

  int maxbridge_all;
  MPI_Allreduce (&maxbridge, &maxbridge_all, 1, MPI_INT, MPI_MAX, world);
  return maxbridge_all + 1;
}

/* ----------------------------------------------------------------------
 *   size of atom nlocal's restart data
 * ------------------------------------------------------------------------- */

int FixAddCapillaryForce::size_restart(int i)
{
  return bridgelist[i].size() + 1;
}
