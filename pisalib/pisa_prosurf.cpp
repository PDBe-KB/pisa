// $Id: pisa_prosurf.cpp $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  ProSurf <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::ProSurf
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2003-2014
//
// =================================================================
//

#include <string.h>
#include <math.h>

#include "pisa_prosurf.h"
#include "mmdb2/mmdb_math_fft.h"
#include "mmdb2/mmdb_tables.h"

namespace pisa  {

  // =================================================================

  mmdb::realtype ASP[nASPs] = {

    0.0,        //  #0  "Other" atoms

    0.016,      //  #1  Protein - Carbon
    0.0408388,  //  #2  Protein - Sulphur
   -0.0114138,  //  #3  Protein - neutral NO
   -0.0370456,  //  #4  Protein - charged N
   -0.0169549,  //  #5  Protein - charged O
   -0.0330343,  //  #6  DNA - Carbon     C_D_S
   -0.0330343,  //  #7  DNA - Carbon     C_D_B
   -0.0330343,  //  #8  DNA - Carbon     C_ARO
    0.014375,   //  #9  DNA - Carbon     C_ALI
   -0.00023691, //  #10 DNA - Nitrogen   N_ARO
   -0.019118,   //  #11 DNA - Nitrogen   N_AMI
    0.0532585,  //  #12 DNA - Oxygen     O_D_P
    0.0231976,  //  #13 DNA - Oxygen     O_D_B
    0.0231976,  //  #14 DNA - Oxygen     O_D_S
    0.0566945,  //  #15 DNA - Oxygen     O_BYL
   -0.0234746   //  #16 DNA - Phosphorus P_DNA

  };


  // =========================  ProSurf  ============================

  ProSurf::ProSurf() : mmdb::io::Stream()  {
    InitProSurf();
  }

  ProSurf::ProSurf ( mmdb::io::RPStream Object )
         : mmdb::io::Stream(Object)  {
    InitProSurf();
  }

  ProSurf::~ProSurf()  {
    FreeSCMemory();
  }


  void  ProSurf::ResetResults()  {

    selHndSurf1 = 0;   // sel-n handle for surface atoms of 1st structure
    selHndSurf2 = 0;   // sel-n handle for surface atoms of 2nd structure
    selHndInt1  = 0;   // sel handle for interface atoms of 1st structure
    selHndInt2  = 0;   // sel handle for interface atoms of 2nd structure

    surfArea1   = 0.0; // surface area of 1st structure
    surfArea2   = 0.0; // surface area of 2nd structure
    surfArea12  = 0.0; // surface area of both structures in contact
    intArea1    = 0.0; // area of the interface of 1st structure
    intArea2    = 0.0; // area of the interface of 2nd structure

  }


  void  ProSurf::InitProSurf()  {

    ResetResults();

    scx = NULL;
    scy = NULL;
    scz = NULL;
    sca = NULL;
    scn = 0;

    codeNo = 36;

    max_radius     = 4.0;  // a generous estimate for maximal atom radius

    probe_radius   = 1.4;  // prove (solvent) sphere radius, A
    excludeSolvent = true;

  }

  void ProSurf::FreeSCMemory()  {
    mmdb::FreeVectorMemory ( scx,0 );
    mmdb::FreeVectorMemory ( scy,0 );
    mmdb::FreeVectorMemory ( scz,0 );
    mmdb::FreeVectorMemory ( sca,0 );
    scn = 0;
  }

  mmdb::realtype ProSurf::getMaxGap()  {
    return 2.0*probe_radius;
  }

  void  ProSurf::setCodeNo ( int code_no )  {
    codeNo = code_no;
    calcSphericalCode();
  }

  void  ProSurf::setProbeRadius ( mmdb::realtype probeRadius )  {
    probe_radius = probeRadius;
  }

  void  ProSurf::setExcludeSolvent ( bool exclSolvent )  {
    excludeSolvent = exclSolvent;
  }

  #define minCodeNo  6

  void ProSurf::calcSphericalCode()  {
  mmdb::realtype dalpha,dbeta,beta,z,xy,a;
  int      i,j,k,nr;

    FreeSCMemory();
    if (codeNo<minCodeNo)  codeNo = minCodeNo;

    dalpha = mmdb::Pi/(codeNo-1);

    scn = 0;
    for (i=1;i<=codeNo;i++)
      if ((i==1) || (i==codeNo))
            scn++;
      else  scn  += mmdb::IMax(minCodeNo,mmdb::mround(codeNo*sin((i-1)*dalpha)));

    mmdb::GetVectorMemory ( scx,scn,0 );
    mmdb::GetVectorMemory ( scy,scn,0 );
    mmdb::GetVectorMemory ( scz,scn,0 );
    mmdb::GetVectorMemory ( sca,scn,0 );

    k = 0;
    for (i=1;i<=codeNo;i++)
      if (i==1)  {
        scx[k] = 0.0;
        scy[k] = 0.0;
        scz[k] = 1.0;
        sca[k] = 2.0*mmdb::Pi*(1.0-cos(dalpha/2.0));
        k++;
      } else if (i==codeNo)  {
        scx[k] = 0.0;
        scy[k] = 0.0;
        scz[k] = -1.0;
        sca[k] = 2.0*mmdb::Pi*(1.0-cos(dalpha/2.0));
        k++;
      } else  {
        nr     = mmdb::IMax(minCodeNo,mmdb::mround(codeNo*sin((i-1)*dalpha)));
        dbeta  = 2.0*mmdb::Pi/nr;
        z      = cos((i-1)*dalpha);
        xy     = sin((i-1)*dalpha);
        a      = 2.0*mmdb::Pi*(cos((i-1.5)*dalpha)-cos((i-0.5)*dalpha))/nr;
        for (j=0;j<nr;j++)  {
          beta   = j*dbeta;
          scx[k] = xy*cos(beta);
          scy[k] = xy*sin(beta);
          scz[k] = z;
          sca[k] = a;
          k++;
        }
      }

    a = 0.0;
    for (i=0;i<scn;i++)
      a += sca[i];

  }


  void ProSurf::GetAtomSet ( mmdb::PManager M, mmdb::cpstr CID,
                             int & selHnd,    mmdb::PPAtom & atom,
                             int & nAtoms )  {
  int i,j;

    selHnd = M->NewSelection();
    if (!CID)
      M->Select ( selHnd,mmdb::STYPE_ATOM,"/1",mmdb::SKEY_NEW );
    else if (!CID[0])
      M->Select ( selHnd,mmdb::STYPE_ATOM,"/1",mmdb::SKEY_NEW );
    else
      M->Select ( selHnd,mmdb::STYPE_ATOM,CID ,mmdb::SKEY_NEW );

    M->GetSelIndex ( selHnd,atom,nAtoms );
    if (excludeSolvent)  {
      j = 0;
      for (i=0;i<nAtoms;i++)
        if (atom[i]->isSolvent())  {
          j = 1;
          M->SelectAtom ( selHnd,atom[i],mmdb::SKEY_CLR,false );
        }
      if (j>0)  {
        M->MakeSelIndex ( selHnd );
        M->GetSelIndex  ( selHnd,atom,nAtoms );
      }
    }

  }


  PROSURF_RC  ProSurf::calcSurface ( mmdb::PManager M,
                                     PMolRefIndex   molRef,
                                     mmdb::cpstr    CID )  {
  mmdb::PPAtom  atom;
  int           selHnd,nAtoms;
  PROSURF_RC    rc;

    GetAtomSet ( M,CID,selHnd,atom,nAtoms );

    if (nAtoms<=0)  {
      ResetResults();
      rc = PROSURF_noAtoms;
    } else
      rc = calcSurface ( atom,nAtoms,molRef,NULL );

    M->DeleteSelection ( selHnd );

    return rc;

  }


  PROSURF_RC ProSurf::calcSurface ( mmdb::PPAtom  atom,
                                    int           nAtoms,
                                    PMolRefIndex  molRef,
                                    mmdb::rvector atomSAS )  {
  mmdb::PManager M;
  mmdb::PBrick   brick;
  mmdb::rvector  rA;
  mmdb::ovector  scb;
  mmdb::realtype ri,rj,rij, xi,yi,zi, xij,yij,zij, dx,dy,dz, d2;
  mmdb::realtype sarea;
  int            i,j,k, ix,iy,iz, bx,by,bz, nAcc;
  PROSURF_RC     rc;

    //  1. General initializations and checks

    ResetResults();

    if (nAtoms<=0)  return PROSURF_noAtoms;

    M = NULL;
    for (i=0;(i<nAtoms) && (!M);i++)
      if (atom[i])
        M = mmdb::PManager(atom[i]->GetCoordHierarchy());
    if (!M)  return PROSURF_noManager;

    if (!scx)  calcSphericalCode();

    rc = PROSURF_Ok;

    //  2. Count non-covered points of the spherical codes superimposed
    //     on all atoms. Use bricking algorithm for speed.

    M->MakeBricks ( atom,nAtoms,1.5*(max_radius+probe_radius),
                    2.0*(max_radius+probe_radius) );

    mmdb::GetVectorMemory ( rA,nAtoms,0 );
    for (i=0;i<nAtoms;i++)
      if (atom[i])
        rA[i] = molRef->getAtomRadius(atom[i]) + probe_radius;

    if (atomSAS)
      for (i=0;i<nAtoms;i++)
        atomSAS[i] = 0.0;

    mmdb::GetVectorMemory ( scb,scn,0 );
    selHndSurf1 = M->NewSelection();
    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        M->GetBrickCoor ( atom[i],bx,by,bz );
        for (k=0;k<scn;k++)
          scb[k] = true;  // initially mark all surface as accessible
        nAcc = scn;
        ri   = rA[i];
        xi   = atom[i]->x;
        yi   = atom[i]->y;
        zi   = atom[i]->z;
        for (ix=bx-1;(ix<=bx+1) && nAcc;ix++)
          for (iy=by-1;(iy<=by+1) && nAcc;iy++)
            for (iz=bz-1;(iz<=bz+1) && nAcc;iz++)  {
              brick = M->GetBrick ( ix,iy,iz );
              if (brick)  {
                for (j=0;(j<brick->nAtoms) && nAcc;j++)
                  if (brick->id[j]!=i)  {
                    rj  = rA[brick->id[j]];
                    xij = xi - brick->atom[j]->x;
                    yij = yi - brick->atom[j]->y;
                    zij = zi - brick->atom[j]->z;
                    d2  = xij*xij + yij*yij + zij*zij;
                    rij = ri + rj;
                    if (d2<=rij*rij)  {
                      rj = rj*rj + 0.00001;  // take care of round-offs
                      for (k=0;(k<scn) && nAcc;k++)
                        if (scb[k])  {
                          dx = xij + ri*scx[k];
                          dy = yij + ri*scy[k];
                          dz = zij + ri*scz[k];
                          d2 = dx*dx + dy*dy + dz*dz;
                          if (d2<=rj)  {  // kth SC point is inaccessible
                            scb[k] = false;
                            nAcc--;
                          }
                        }
                    }
                  }
              }
            }
        if (nAcc)  {
          // atom is accessible, select it
          M->SelectAtom ( selHndSurf1,atom[i],mmdb::SKEY_OR,false );
          // count the accessible area
          sarea = 0.0;
          for (k=0;k<scn;k++)
            if (scb[k])  sarea += sca[k];
          sarea *= ri*ri;
          surfArea1 += sarea;
          if (atomSAS)  atomSAS[i] = sarea;
        }
      }

    mmdb::FreeVectorMemory ( scb,0 );
    mmdb::FreeVectorMemory ( rA ,0 );

    //  3. Make index of accessible atoms

    M->MakeSelIndex ( selHndSurf1 );

    //  4. Release allocated bricks

    M->RemoveBricks();

    return rc;

  }


  PROSURF_RC ProSurf::calcInterface ( mmdb::PPAtom  atom1,
                                      int           nAtoms1,
                                      mmdb::PPAtom  atom2,
                                      int           nAtoms2,
                                      PMolRefIndex  molRef,
                                      SURFSEL_KEY   surfFlag,
                                      mmdb::rvector atomSAS1,
                                      mmdb::rvector intAtomSAS1,
                                      mmdb::rvector atomSAS2,
                                      mmdb::rvector intAtomSAS2 )  {
  mmdb::PManager M;
  mmdb::PPAtom   atom;
  mmdb::PBrick   brick;
  mmdb::rvector  rA;
  mmdb::ovector  scb1,scb2;
  mmdb::realtype ri,rj,rij, xi,yi,zi, xij,yij,zij, dx,dy,dz, d2;
  mmdb::realtype darea, sarea,iarea,carea;
  int            nAtoms, i,j,k, ix,iy,iz, bx,by,bz, atmid;
  int            nAcc1,nAcc2;
  bool           done,selSurface,checkOnSurface,doAtom;
  PROSURF_RC     rc;

    //  1.  General initializations and checks

    i = selHndSurf1;
    j = selHndSurf2;

    ResetResults();  // zeroes the surfaces and selection handles

    M = NULL;
    for (i=0;(i<nAtoms1) && (!M);i++)
      if (atom1[i])
        M = mmdb::PManager(atom1[i]->GetCoordHierarchy());
    if (!M)  return PROSURF_noManager;

    if (surfFlag==SURFF_Selected)  {
      selHndSurf1    = i;
      selHndSurf2    = j;
      selSurface     = (selHndSurf1<=0) || (selHndSurf2<=0);
      checkOnSurface = !selSurface;
    } else  {
      selSurface     = (surfFlag!=SURFF_dontSelect);
      checkOnSurface = false;
    }

    if (!scx)  calcSphericalCode();

    if ((nAtoms1<=0) || (nAtoms2<=0))  return PROSURF_noAtoms;

    if (atomSAS1)
      for (i=0;i<nAtoms1;i++)
        atomSAS1[i] = 0.0;

    if (intAtomSAS1)
      for (i=0;i<nAtoms1;i++)
        intAtomSAS1[i] = 0.0;

    if (atomSAS2)
      for (i=0;i<nAtoms2;i++)
        atomSAS2[i] = 0.0;

    if (intAtomSAS2)
      for (i=0;i<nAtoms2;i++)
        intAtomSAS2[i] = 0.0;

    rc = PROSURF_Ok;

    //  2. Build a common atom array

    nAtoms = nAtoms1 + nAtoms2;
    atom = new mmdb::PAtom[nAtoms];
    j = 0;
    for (i=0;i<nAtoms1;i++)
      atom[j++] = atom1[i];
    for (i=0;i<nAtoms2;i++)
      atom[j++] = atom2[i];

    //  3. Count non-covered points of the spherical codes superimposed
    //     on all atoms. Use bricking algorithm for speed.

    M->MakeBricks ( atom,nAtoms,1.5*(max_radius+probe_radius),
                    2.0*(max_radius+probe_radius) );

    mmdb::GetVectorMemory ( rA,nAtoms,0 );
    for (i=0;i<nAtoms;i++)
      if (atom[i])
        rA[i] = molRef->getAtomRadius(atom[i]) + probe_radius;

    mmdb::GetVectorMemory ( scb1,scn,0 );
    mmdb::GetVectorMemory ( scb2,scn,0 );

    if (selSurface)  {
      selHndSurf1 = M->NewSelection();
      selHndSurf2 = M->NewSelection();
    }
    selHndInt1 = M->NewSelection();
    selHndInt2 = M->NewSelection();

    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        if (checkOnSurface)  {
          // if surface selection is provided, then we look for interface
          // atoms only among those on the surface
          if (i<nAtoms1)  doAtom = atom[i]->isInSelection(selHndSurf1);
                    else  doAtom = atom[i]->isInSelection(selHndSurf2);
        } else  {
          // if surface selection is not given, we look through all atoms
          doAtom = true;
        }
        if (doAtom)  {

          xi = atom[i]->x;
          yi = atom[i]->y;
          zi = atom[i]->z;
          ri = rA[i];  // VdW radius of ith atom

          M->GetBrickCoor ( atom[i],bx,by,bz );

          //  scb1[k]/scb2[k] refer segments on the sphere of atom
          //  atom[i] overlaped by atoms belonging to 1st/2nd structures
          for (k=0;k<scn;k++)  {
            scb1[k] = true;  // initially mark all surface as accessible
            scb2[k] = true;
          }
          nAcc1 = scn;  // number of surface segments not covered by
          nAcc2 = scn;  // spheres belonging to atoms from 1st/2nd
                        // structures
          // search for overlaps in the neighbouring bricks
          done = false;
          for (ix=bx-1;(ix<=bx+1) && (!done);ix++)
            for (iy=by-1;(iy<=by+1) && (!done);iy++)
              for (iz=bz-1;(iz<=bz+1) && (!done);iz++)  {
                brick = M->GetBrick ( ix,iy,iz );
                if (brick)  {
                  for (j=0;(j<brick->nAtoms) && (!done);j++)
                    if (atom[i]!=brick->atom[j])  {
                      atmid = brick->id[j];
                      if (((atmid<nAtoms1) && (nAcc1>0)) ||
                          ((atmid>=nAtoms1) && (nAcc2>0)))  {
                        rj  = rA[atmid];
                        xij = xi - brick->atom[j]->x;
                        yij = yi - brick->atom[j]->y;
                        zij = zi - brick->atom[j]->z;
                        d2  = xij*xij + yij*yij + zij*zij;
                        rij = ri + rj;
                        if (d2<=rij*rij)  {
                          rj = rj*rj + 0.00001; // account of round-offs
                          if (atmid<nAtoms1)  {
                            // brick->atom[j] is an atom from 1st
                            // structure. Look for segments of
                            // atom[i]-th sphere that overlap with
                            // sphere of atom brick->atom[j]
                            for (k=0;(k<scn) && (nAcc1>0);k++)
                              if (scb1[k])  {
                                dx = xij + ri*scx[k];
                                dy = yij + ri*scy[k];
                                dz = zij + ri*scz[k];
                                d2 = dx*dx + dy*dy + dz*dz;
                                if (d2<=rj)  {
                                  // kth SC point is inaccessible
                                  scb1[k] = false;
                                  nAcc1--;
                                }
                              }
                          } else  {
                            //   brick->atom[j] is an atom from 2nd
                            // structure.
                            //   Look for segments of atom[i]-th sphere
                            // that overlap with sphere of atom
                            // brick->atom[j]
                            for (k=0;(k<scn) && (nAcc2>0);k++)
                              if (scb2[k])  {
                                dx = xij + ri*scx[k];
                                dy = yij + ri*scy[k];
                                dz = zij + ri*scz[k];
                                d2 = dx*dx + dy*dy + dz*dz;
                                if (d2<=rj)  {
                                  // kth SC point is inaccessible
                                  scb2[k] = false;
                                  nAcc2--;
                                }
                              }
                          }
                          done = (nAcc1<=0) && (nAcc2<=0);
                        }
                      }
                    }
                }
              }

          d2    = ri*ri; // VdW surface factor
          iarea = 0.0;   // interface area between the structures
          sarea = 0.0;   // surface area of 1st structure
          carea = 0.0;   // common surface area of both structures

          if ((i<nAtoms1) && (nAcc1>0))  {
            // the atom is found on the surface of 1st structure

            for (k=0;k<scn;k++)
              if (scb1[k])  {
                // kth segment is not covered by 1st structure atoms.
                // Add its area to the surface area of 1st structure
                darea  = sca[k]*d2;
                sarea += darea;
                if (scb2[k])  {
                  // the segment is not covered by 2nd structure atoms
                  // Add its area to the common surface area
                  carea += darea;
                } else  {
                  // the segment is covered by 2nd structure atoms.
                  // Add its area to the interface area of 1st structure
                  iarea += darea;
                }
              }
            if ((sarea>0.0) && selSurface)
              M->SelectAtom ( selHndSurf1,atom[i],mmdb::SKEY_OR,false );
            surfArea1  += sarea;
            surfArea12 += carea;
            if (iarea>0.0)  {
              M->SelectAtom ( selHndInt1,atom[i],mmdb::SKEY_OR,false );
              intArea1 += iarea;
            }

            if (atomSAS1   )  atomSAS1   [i] = sarea;
            if (intAtomSAS1)  intAtomSAS1[i] = iarea;

          } else if ((i>=nAtoms1) && (nAcc2>0))  {
            // the atom is found on the surface of 2nd structure

            for (k=0;k<scn;k++)
              if (scb2[k])  {
                // kth segment is not covered by 2nd structure atoms.
                // Add its area to the surface area of 1st structure
                darea  = sca[k]*d2;
                sarea += darea;
                if (scb1[k])  {
                  // the segment is not covered by 1st structure atoms
                  // Add its area to the common surface area
                  carea += darea;
                } else  {
                  // the segment is covered by 1st structure atoms.
                  // Add its area to the interface area of 1st structure
                  iarea += darea;
                }
              }
            if ((sarea>0.0) && selSurface)
              M->SelectAtom ( selHndSurf2,atom[i],mmdb::SKEY_OR,false );
            surfArea2  += sarea;
            surfArea12 += carea;
            if (iarea>0.0)  {
              M->SelectAtom ( selHndInt2,atom[i],mmdb::SKEY_OR,false );
              intArea2 += iarea;
            }

            k = i-nAtoms1;
            if (atomSAS2   )  atomSAS2   [k] = sarea;
            if (intAtomSAS2)  intAtomSAS2[k] = iarea;

          }

        }

      }

    mmdb::FreeVectorMemory ( scb1,0 );
    mmdb::FreeVectorMemory ( scb2,0 );
    mmdb::FreeVectorMemory ( rA  ,0 );
    delete[] atom;

    //  4. Make indexes of accessible atoms

    if (selSurface) {
      M->MakeSelIndex ( selHndSurf1 );
      M->MakeSelIndex ( selHndSurf2 );
      if (M->GetSelLength(selHndSurf1)<=0)  {
        M->DeleteSelection(selHndSurf1);
        selHndSurf1 = 0;
      }
      if (M->GetSelLength(selHndSurf2)<=0)  {
        M->DeleteSelection(selHndSurf2);
        selHndSurf2 = 0;
      }
    }
    M->MakeSelIndex ( selHndInt1  );
    M->MakeSelIndex ( selHndInt2  );
    if (M->GetSelLength(selHndInt1)<=0)  {
      M->DeleteSelection(selHndInt1);
      selHndInt1 = 0;
    }
    if (M->GetSelLength(selHndInt2)<=0)  {
      M->DeleteSelection(selHndInt2);
      selHndInt2 = 0;
    }

    //  5. Release allocated bricks

    M->RemoveBricks();

    return rc;

  }


  PROSURF_RC  ProSurf::calcInterface ( mmdb::PManager M,
                                       mmdb::cpstr    CID1,
                                       mmdb::cpstr    CID2,
                                       PMolRefIndex   molRef,
                                       SURFSEL_KEY    surfFlag )  {
  mmdb::PPAtom atom1,atom2;
  int          selHnd1,selHnd2, nAtoms1,nAtoms2;
  PROSURF_RC   rc;

    GetAtomSet ( M,CID1,selHnd1,atom1,nAtoms1 );
    GetAtomSet ( M,CID2,selHnd2,atom2,nAtoms2 );

    if ((nAtoms1<=0) || (nAtoms2<=0))  {
      ResetResults();
      rc = PROSURF_noAtoms;
    } else
      rc = calcInterface ( atom1,nAtoms1,atom2,nAtoms2,
                           molRef,surfFlag,NULL,NULL,NULL,NULL );

    M->DeleteSelection ( selHnd1 );
    M->DeleteSelection ( selHnd2 );

    return rc;

  }



  // ========================  SolvEnergy  ==========================


  SolvEnergy::SolvEnergy() : mmdb::io::Stream()  {
    InitSolvEnergy();
  }

  SolvEnergy::SolvEnergy ( mmdb::io::RPStream Object )
            : mmdb::io::Stream(Object)  {
    InitSolvEnergy();
  }

  SolvEnergy::~SolvEnergy()  {}


  void SolvEnergy::ResetResults()  {
  int i;

    solvEnergy = 0.0;  // solvation energy of folding:
    DeltaG     = 0.0;  // solvation energy gain of interface
    PValue     = 0.0;  // P-value of the interface solvation energy gain
    aveDeltaG  = 0.0;  // average Delta G for the interface
    probDeltaG = 0.0;  // probable Delta G for the interface
    for (i=0;i<nASPs;i++)
      aspArea[i] = 0.0; // partial atom surface areas

  }


  void SolvEnergy::InitSolvEnergy()  {
    ResetResults();
  }


  int aspType ( mmdb::PAtom a, mmdb::cpstr aname,
                int typeCharged, int typeNeutral )  {
  mmdb::PAtom a1;
    a1 = a->GetResidue()->GetAtom ( aname,NULL,NULL );
    if (a1)  {
      // su11 is temporarily assigned for atom SAS
      if (a->su11>a1->su11)  return typeCharged;
                                     // the most exposed one is charged
      return typeNeutral;
    } else
      return typeCharged;
  }

  int  aspType ( mmdb::PAtom a, mmdb::cpstr aname1, mmdb::cpstr aname2,
                 int typeCharged, int typeNeutral )  {
  mmdb::PAtom  a1,a2;
    a1 = a->GetResidue()->GetAtom ( aname1,NULL,NULL );
    a2 = a->GetResidue()->GetAtom ( aname2,NULL,NULL );
    if (a1 && a2)  {
      // su11 is temporarily assigned for atom SAS
      if (a->su11==mmdb::RMax(a->su11,mmdb::RMax(a1->su11,a2->su11)))
         return typeCharged;  // the most exposed one is charged
      return typeNeutral;
    } else if (a1)  {
      if (a->su11>a1->su11)  return typeCharged;
      return typeNeutral;
    } else if (a2)  {
      if (a->su11>a2->su11)  return typeCharged;
      return typeNeutral;
    } else
      return typeCharged;
  }

  void  writeASPs ( mmdb::io::RFile f )  {
  int        i;
  mmdb::byte Version=1;
    f.WriteByte ( &Version );
    for (i=0;i<nASPs;i++)
      f.WriteReal ( &(ASP[i]) );
  }


  void  readASPs ( mmdb::io::RFile f )  {
  int        i;
  mmdb::byte Version;
    f.ReadByte ( &Version );
    for (i=0;i<nASPs;i++)
      f.ReadReal ( &(ASP[i]) );
  }


  int specialAspType ( mmdb::PAtom a )  {
  mmdb::ResName name;

    if (!strcmp(a->element," N"))  {
      if (!strcmp(a->name," NE2"))  {
        // in non-His, this atom is uncharged
        if (strcmp(a->GetResName(),"HIS"))
          return ASP_NeutralNO;
        // else check for the presence and exposure of ND1 in His
        return aspType ( a," ND1",ASP_ChargedN,ASP_NeutralNO );
      }
      // ND1 appears only in His, check for the presence and exposure
      // of NE2
      if (!strcmp(a->name," ND1"))
        return aspType ( a," NE2",ASP_ChargedN,ASP_NeutralNO );
      // remaining nitrogens to check are NE, NH1 and NH2, which are
      // met only in Arg. Return the most exposed of them as the charged
      // one, others uncharged
      if (!strcmp(a->name," NE "))
        return aspType ( a," NH1"," NH2",ASP_ChargedN,ASP_NeutralNO );
      if (!strcmp(a->name," NH1"))
        return aspType ( a," NE "," NH2",ASP_ChargedN,ASP_NeutralNO );
      if (!strcmp(a->name," NH2"))
        return aspType ( a," NE "," NH1",ASP_ChargedN,ASP_NeutralNO );
      return ASP_ChargedN;
    }

    if (!strcmp(a->element," O"))  {
      if (a->name[3]=='1')  {
        strcpy ( name,a->name );
        name[3] = '2';
        return aspType ( a,name,ASP_ChargedO,ASP_NeutralNO );
      }
      if (a->name[3]=='2')  {
        strcpy ( name,a->name );
        name[3] = '1';
        return aspType ( a,name,ASP_ChargedO,ASP_NeutralNO );
      }
      return ASP_NeutralNO;
    }

    return ASP_Other;

  }


  mmdb::realtype SolvEnergy::calcSurfSolvEnergy (
                                  mmdb::PPAtom   atom,
                                  mmdb::rvector  atomSAS,
                                  int            nAtoms,
                                  mmdb::rvector  atomSE,
                                  PMolRefIndex   molRef )  {
  //   This function calculates the solvation energy of folding
  // for proteins, or other structures where accessible solvent areas
  // in reference (unfolded) state is found in MolRef. The accessible
  // surface areas for standard asp types is stored in aspArea[].
  // Other asp types, if found in MolRef, are taken into account,
  // however their areas are not returned.
  //   Areas of standard asp types (amino acids and nucleic acids) is
  // used in the asp fitting procedure. Other asp types are assigned
  // manually in MolRef.
  PAtmRefData     ARefData;
  mmdb::realtype  refASA;
  int             i,t;

    //  1.  Atom's accessible surface areas (calculated by ProSurf)
    //      are temporarily stored in atom[]->su11.
    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        refASA        = atomSAS[i];
        atomSAS[i]    = atom[i]->su11;
        atom[i]->su11 = refASA;
      }

    if (atomSE)
      for (i=0;i<nAtoms;i++)
        atomSE[i] = 0.0;

    for (i=0;i<nASPs;i++)
      aspArea[i] = 0.0;

    solvEnergy = 0.0;

    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        ARefData = molRef->getAtmRefData ( atom[i] );
        if (ARefData)  {
          refASA = ARefData->refASA;
          if (refASA>=-0.5)  {
            refASA = atom[i]->su11 - refASA;
            t = ARefData->aspId;
            if (t==ASP_Special)
              t = specialAspType ( atom[i] );
            if (t<0)  t = 0;
            aspArea[t] += refASA;
            if (t<=0)  {
              refASA      = ARefData->asp*atom[i]->su11;
              solvEnergy += refASA;
              if (atomSE)
                atomSE[i] = refASA;
            } else if (atomSE)
              atomSE[i] = ASP[t]*atom[i]->su11;
          }
        }
      }

    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        refASA        = atom[i]->su11;
        atom[i]->su11 = atomSAS[i];
        atomSAS[i]    = refASA;
      }

    for (i=0;i<nASPs;i++)
      solvEnergy += aspArea[i]*ASP[i];

    return solvEnergy;

  }


  mmdb::realtype SolvEnergy::calcIntSolvEnergy (
                               mmdb::PPAtom   atom,
                               mmdb::rvector  atomSAS,
                               mmdb::rvector  intSAS,
                               int            nAtoms,
                               mmdb::rvector  atomSE,
                               PMolRefIndex   molRef )  {
  //   The solvation energy gain at interface formation is calculated
  // as a difference between solvation energy of atoms atom[i],
  // i=0..nAtoms-1, with solvent-accessible surfaces given in
  // atomSAS[i], and their solvation energy when SAS is decreased by
  // areas given in intSAS[i], the letter being the SAS buried by
  // the interface.
  //   Although a previously-calculated value for solvation energy
  // without interface could be used, we make a full-range of
  // calculations here because a) reorientation of the molecule has
  // a numerical impact on SAS and b) interface formation often results
  // in reallocating charges between pairs of chargable oxygens and
  // nitrogens in aminoacids. Change in SAS at reorientation is
  // important to take account of, because the solvation energy gain
  // represents a difference effect. Calculating both the non-interfaced
  // and interfaced states at the same orientation drastically reduces
  // the numerical effect of changing SASs on the interface solvation
  // energy gain.
  PAtmRefData     ARefData;
  mmdb::realtype  sAspArea[nASPs];
  mmdb::realtype  refASA;
  int             i,t;

    //  1. Make a temporary storage for SASs in Atom structures, which
    //     will be used for calculating charges on oxygens and nitrogens.
    //     Vector atomSAS will keep original data from atom structures.

    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        refASA        = atomSAS[i];
        atomSAS[i]    = atom[i]->su11;
        atom[i]->su11 = refASA;
      }

    if (atomSE)
      for (i=0;i<nAtoms;i++)
        atomSE[i] = 0.0;

    //  2. Calculate solvation energy of folding in non-interfaced state

    // This is solvation energy of folding in non-interfaced state at
    // given orientation of the molecule (as reflected in atomSAS).
    // This value is retained by the class.
    solvEnergy = 0.0;

    for (i=0;i<nASPs;i++)
      sAspArea[i] = 0.0;

    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        ARefData = molRef->getAtmRefData ( atom[i] );
        if (ARefData)  {
          refASA = ARefData->refASA;
          if (refASA>=-0.5)  {
            refASA = atom[i]->su11 - refASA;
            t = ARefData->aspId;
            if (t==ASP_Special)
              t = specialAspType ( atom[i] );
            if (t<0)  t = 0;
            sAspArea[t] += refASA;
            if (t<=0) {
              refASA      = ARefData->asp*atom[i]->su11;
              solvEnergy += refASA;
              if (atomSE)
                atomSE[i] = refASA;
            } else if (atomSE)
              atomSE[i] = ASP[t]*refASA;
          }
        }
      }

    for (i=0;i<nASPs;i++)
      solvEnergy += sAspArea[i]*ASP[i];


    //  3. Calculate solvation energy of folding in the interfaced state

    //  Substract area, buried by the interface, from atom SAS. This will
    //  result in reallocating charges on some oxygens and nitrogens that
    //  border the interface, therefore we will repeat solvation energy
    //  calculations in full.
    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        atom[i]->su11 -= intSAS[i];
        // the following line corrects tiny round-offs
        if (atom[i]->su11<1.0e-8)  atom[i]->su11 = 0.0;
      }

    DeltaG = 0.0;

    for (i=0;i<nASPs;i++)
      aspArea[i] = 0.0;

    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        ARefData = molRef->getAtmRefData ( atom[i] );
        if (ARefData)  {
          refASA = ARefData->refASA;
          if (refASA>=-0.5)  {
            refASA = atom[i]->su11 - refASA;
            t = ARefData->aspId;
            if (t==ASP_Special)
              t = specialAspType ( atom[i] );
            if (t<0)  t = 0;
            aspArea[t] += refASA;
            if (t<=0) {
              refASA  = ARefData->asp*atom[i]->su11;
              DeltaG += refASA;
              if (atomSE)
                atomSE[i] -= refASA;
            } else if (atomSE)
              atomSE[i] -= ASP[t]*refASA;
          }
        }
      }

    for (i=0;i<nASPs;i++)  {
      DeltaG += aspArea[i]*ASP[i];
      aspArea[i] -= sAspArea[i];
    }
    DeltaG -= solvEnergy;

    //  4. Restore the atom structures and vector atomSAS

    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        refASA        = atom[i]->su11 + intSAS[i];
        atom[i]->su11 = atomSAS[i];
        atomSAS[i]    = refASA;
      }


    return DeltaG;

  }


  mmdb::realtype SolvEnergy::calcIntPValue (
                                mmdb::PPAtom    atom,
                                mmdb::rvector   atomSAS,
                                int             nAtoms,
                                int             nIntAtoms,
                                mmdb::realtype  intDeltaG,
                                PMolRefIndex    molRef )  {
  PAtmRefData     ARefData;
  mmdb::rvector   surfEn;
  mmdb::realtype  B;
  int             i,t,nSurfAtoms;

    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        B             = atomSAS[i];
        atomSAS[i]    = atom[i]->su11;
        atom[i]->su11 = B;
      }

    mmdb::GetVectorMemory ( surfEn,nAtoms,0 );

    nSurfAtoms = 0;
    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        ARefData = molRef->getAtmRefData ( atom[i] );
        if (ARefData)  {
          if ((atom[i]->su11>0.0) && (ARefData->refASA>=-0.5))  {
            t = ARefData->aspId;
            if (t==ASP_Special)
              t = specialAspType ( atom[i] );
            if (t>0)
                 surfEn[nSurfAtoms++] = -atom[i]->su11*ASP[t];
            else surfEn[nSurfAtoms++] = -atom[i]->su11*ARefData->asp;
          }
        }
      }

    calcPValue ( surfEn,nSurfAtoms,nIntAtoms,intDeltaG );
    mmdb::FreeVectorMemory ( surfEn,0 );

    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        B             = atom[i]->su11;
        atom[i]->su11 = atomSAS[i];
        atomSAS[i]    = B;
      }

    return PValue;

  }


  void  SolvEnergy::calcPValue ( mmdb::rvector surfEn, int nSurfAtoms,
                                 int nIntAtoms, mmdb::realtype intDeltaG )  {
  //   P-value is estimated in thermodynamics limit here, which
  // corresponds to a negligeably small selection pool (or number
  // of interface atoms) as comparing to the total number of
  // surface atoms.
  mmdb::realtype Emin,Emax,dE,sumP;
  mmdb::rvector  pE;
  int      nBins,nConv,nConv2;
  int      i,k,m;

    if (nSurfAtoms<=0)  {
      PValue    = -1.0;   // indication of an error
      aveDeltaG = 0.0;
      return;
    }
    if (nIntAtoms<=0)  {
      PValue    = -2.0;   // indication of an error
      aveDeltaG = 0.0;
      return;
    }

    Emin =  mmdb::MaxReal;
    Emax = -mmdb::MaxReal;
    aveDeltaG = 0.0;
    for (i=0;i<nSurfAtoms;i++)  {
      if (surfEn[i]<Emin)  Emin = surfEn[i];
      if (surfEn[i]>Emax)  Emax = surfEn[i];
      aveDeltaG += surfEn[i];
    }

    //  Calculate average expected solvation energy gain for interface
    //  of nIntAtoms;
    aveDeltaG /= nSurfAtoms;
    aveDeltaG *= nIntAtoms;

    if (fabs(Emax-Emin)<1.0e-6)  {
      if (intDeltaG>=aveDeltaG)  PValue = 0.0;
                           else  PValue = 1.0;
      return;
    }

    if (intDeltaG>nIntAtoms*Emax)  {
      PValue = 0.0;
      return;
    }
    if (intDeltaG<nIntAtoms*Emin)  {
      PValue = 1.0;
      return;
    }

    //  Calculate optimal number of energy bins
    nConv = 1024*128;  // maximal length of the resulted convolution
    nBins = mmdb::IMin(nSurfAtoms,nConv/(2*nIntAtoms));
    //  Calculate the needed number of convolution bins
    i      = nBins*nIntAtoms;
    nConv2 = 2;
    while (nConv2<i)  nConv2 *= 2;
    nConv  = nConv2*2;  // optimal length of convolution

    //  Calculate the energy probability spectrum in wrap-around order
    //  prepared for FFT
    mmdb::GetVectorMemory ( pE,nConv,0 );
    for (i=0;i<nConv;i++)
      pE[i] = 0.0;

    //  Fill energy bins

    dE = mmdb::RMax(fabs(Emax),fabs(Emin))/nBins;
    for (i=0;i<nSurfAtoms;i++)  {
      k = mmdb::mround(surfEn[i]/dE);
      if (k<0)  {
        k += nConv;
        if (k<nConv2)  k = nConv2;
      } else if (k>=nConv2)
        k = nConv2-1;
      pE[k] += 1.0;
    }

    sumP = 0.0;
    for (i=0;i<nConv;i++)
      sumP += pE[i];
    for (i=0;i<nConv;i++)
      pE[i] /= sumP;

    //  Make m-convolution
    mmdb::math::mConvolve ( pE,nConv,nIntAtoms-1 );

    //  Calculate integral distribution in wrap-around order

    if (pE[nConv2]<0.0)  pE[nConv2] = 0.0;  // just the round-offs
    for (i=nConv2+1;i<nConv;i++)
      if (pE[i]>=0.0)  pE[i] += pE[i-1];
                 else  pE[i]  = pE[i-1];
    if (pE[nConv-1]>=0.0)  pE[0] += pE[nConv-1];
                     else  pE[0]  = pE[nConv-1];
    for (i=1;i<nConv2;i++)
      if (pE[i]>=0.0)  pE[i] += pE[i-1];
                 else  pE[i]  = pE[i-1];

    //  Find P-value
    k = mmdb::mround(intDeltaG/dE);
    if (k<0)  {
      k += nConv;
      if (k<nConv2)  k = nConv2;
    } else if (k>=nConv2)
      k = nConv2-1;

    PValue = pE[k]/pE[nConv2-1]; // probability that solvation energy
                                 // gain for nIntAtoms ineterface atoms
                                 // may be greater than intDeltaG,
                                 // i.e. that a lower DeltaG is achieved

    // Calculate probable interface solvation energy gain. Here it is
    // defined as energy for which pE=0.5 and it should theoretically
    // coincide with average energy.
    if (pE[0]<0.5)  {
      i = 0;
      k = nConv2-1;
    } else  {
      i = nConv2;
      k = nConv-1;
    }
    do {
      m = (i+k)/2;
      if (pE[m]>0.5)      k = m;
      else if (pE[m]<0.5) i = m;
      else  {
        i = m;
        k = m;
      }
    } while (i<k-1);

    if (fabs(pE[i]-0.5)<fabs(pE[k]-0.5)) m = i;
                                    else m = k;
    if (m>=nConv2)  m -= nConv;
    probDeltaG = dE*m;

    mmdb::FreeVectorMemory ( pE,0 );

  }


  //  ==================================================================

  #ifndef  __MMDB_MMCIF__
  #include "mmdb_mmcif.h"
  #endif

  mmdb::cpstr ASPname[nASPs] = {
    "Unknown",
    "AA Carbon",
    "AA Sulphur",
    "AA Neutral NO",
    "AA Charged N",
    "AA Charged O",
    "DNA C_D_S",
    "DNA C_D_B",
    "DNA C_ARO",
    "DNA C_ALI",
    "DNA N_ARO",
    "DNA N_AMI",
    "DNA O_D_P",
    "DNA O_D_B",
    "DNA O_D_S",
    "DNA O_BYL",
    "DNA P_DNA"
  };


  void makeASPCIF ( mmdb::mmcif::PLoop mmCIFLoop )  {
  int  i;

    mmCIFLoop->AddLoopTag ( "id"       );
    mmCIFLoop->AddLoopTag ( "asp_name" );
    mmCIFLoop->AddLoopTag ( "asp"      );
    for (i=0;i<nASPs;i++)  {
      mmCIFLoop->AddInteger ( i          );
      mmCIFLoop->AddString  ( ASPname[i] );
      mmCIFLoop->AddReal    ( ASP[i],10  );
    }

  }

}  // namespace mmdb::Pisa
