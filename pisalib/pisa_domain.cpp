// $Id: pisa_domain.cpp $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_domain <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Domain
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2014
//
// =================================================================
//

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "pisa_domain.h"
#include "pisa_defs.h"

#include "mmdb2/mmdb_seqsuperpose.h"
#include "mmdb2/mmdb_math_linalg.h"
#include "mmdb2/mmdb_tables.h"

#include "ssm/ssm_superpose.h"
#include "ssm/ssm_csia.h"

namespace pisa  {

  //  ========================  Domain  =============================

  Domain::Domain()  {
    InitDomain();
  }

  Domain::Domain ( mmdb::io::RPStream Object )
        : mmdb::io::Stream ( Object )  {
    InitDomain();
  }

  Domain::~Domain()  {
    FreeMemory();
  }

  void Domain::InitDomain()  {
  int i,j;

    range      = NULL;
    typeId     = NULL;

    ncsOpNo    = 0;     // NCS operation 0..nNCSOps, 0<->identity
    ncsParent  = -1;    // parential domain for NCS mates
    selHndSurf = 0;
    nAtoms     = 0;     // total number of atoms in the domain
    nRes       = 0;     // total number of residues in the domain
    nRes0      = 0;     // total number of AA/NT residues in the domain
    nSurfAtoms = 0;     // number of surafce atoms in the domain
    nSurfRes   = 0;     // number of surface residues in the domain
    type       = 0;     // domain type 1,2,...
    dclass     = DCLASS_None;  // domain class (protein/dna/ligand)
    symNumber  = 1;     // symmetry number
    assemble   = true;  // true if to be inculded into assembling
    agent      = false; // true if precipitation agent
    mmdb::Mat4Init ( ncs_m ); // NCS matrix
    for (i=0;i<3;i++)  {
      for (j=0;j<3;j++)  {
        itb[i][j] = 0.0;  // inertia tensor basis
        itm[i][i] = 0.0;  // inertia tensor matrix
      }
      itb[i][i] = 1.0;
      itm[i][i] = 1.0;
      im [i]    = 1.0;    // inertia moments
    }
    surfArea = 0.0;
    DeltaG   = 0.0;   // relative solvation energy
    mx       = 0.0;
    my       = 0.0;
    mz       = 0.0;
    for (i=0;i<frameLen;i++)
      for (j=0;j<3;j++)
        frame[i][j] = -mmdb::MaxReal;  // frame not calculated
    weight   = 0.0;
    R        = -1.0;  // domain dimensions not calculated
    asa      = NULL;  // vector[0..nRes-1] of accesible surface area
                      //   of individual residues
    solvEn   = NULL;  // vector[0..nRes-1] of solvation energy changes
                      //   of individual residues
    ligandRes = NULL; // for protein/rna/dna, list of residues that
                      // were recognized as ligands
  }

  void Domain::FreeMemory()  {
    if (range)  {
      delete[] range;
      range = NULL;
    }
    if (typeId)  {
      delete[] typeId;
      typeId = NULL;
    }
    selHndSurf = 0;
    nAtoms     = 0;     // total number of atoms in the domain
    nRes       = 0;     // total number of residues in the domain
    nRes0      = 0;     // total number of AA/NT residues in the domain
    nSurfAtoms = 0;     // number of surafce atoms in the domain
    nSurfRes   = 0;     // number of surface residues in the domain
    type       = 0;     // domain type 1,2,...
    dclass     = DCLASS_None;  // domain class (protein/dna/ligand)
    symNumber  = 1;     // symmetry number
    assemble   = true;  // true if to be inculded into assembling
    surfArea   = 0.0;
    DeltaG     = 0.0;   // relative solvation energy
    mx         = 0.0;
    my         = 0.0;
    mz         = 0.0;
    weight     = 0.0;
    R          = -1.0;
    if (!ncsOpNo)  {
      //   For NCS-mates, only pointers to asa and solvEn vectors
      // of their parents are retained. These should not be
      // repeatedly allocated/deallocated.
      mmdb::FreeVectorMemory ( asa   ,0 );
      mmdb::FreeVectorMemory ( solvEn,0 );
    }
    if (ligandRes)  {
      delete[] ligandRes;
      ligandRes = NULL;
    }
  }

  void Domain::SetChainRange ( mmdb::cpstr chainID )  {
    if (!chainID)
          mmdb::CreateCopy   ( range,"-"         );
    else if (!chainID[0])
          mmdb::CreateCopy   ( range,"-"         );
    else if (!strcmp(chainID," "))
          mmdb::CreateCopy   ( range,"-"         );
    else  mmdb::CreateCopCat ( range,chainID,":" );
    // terminating colon is important for domain selection functions
  }

  void Domain::SetLigandRange ( mmdb::PResidue res )  {
  char          N[100];
  mmdb::ChainID chainID;
  mmdb::pstr    chID;
    chID = res->GetChainID();
    if (!chID)                   strcpy ( chainID,"-" );
    else if (!chID[0])           strcpy ( chainID,"-" );
    else if (!strcmp(chID," "))  strcpy ( chainID,"-" );
                           else  strcpy ( chainID,chID );
    sprintf ( N,"[%s]%s:%i%s",res->name,chainID,res->seqNum,
                              res->insCode );
    mmdb::CreateCopy ( range,N );
  }

  bool Domain::checkLigandRange ( PDomain D )  {
  int  i;
  bool B;
    if (dclass!=DCLASS_Ligand)     return false;
    if (D->dclass!=DCLASS_Ligand)  return false;
    i = -1;
    do  {
      i++;
      B = (range[i]==D->range[i]);
    } while (B && (range[i]!=']') && range[i] && D->range[i]);
    return B;
  }

  void Domain::SetLigandResidues ( mmdb::cpstr ligRes )  {
    mmdb::CreateCopy ( ligandRes,ligRes );
  }

  void Domain::getCID ( mmdb::pstr chID, mmdb::pstr resName, int & resNo,
                        mmdb::pstr insCode )  {
  mmdb::pstr p,p1;
  char c;

    p = mmdb::FirstOccurence ( range,']' );
    if (!p)  {
      // domain range, just a chain ID
      strcpy ( chID,range  );
      strcpy ( resName,"*" );
      resNo = mmdb::ANY_RES;
      strcpy ( insCode,"*" );
    } else  {
      *p = char(0);
      strcpy ( resName,&(range[1]) );
      *p = ']';
      if (p[1]=='-')  {
        chID[0] = char(0);
        p += 3;
      } else  {
        //  allow for many-character chain IDs here - for future!
        p1 = &(p[1]);
        p  = mmdb::FirstOccurence ( p1,':' );
        if (!p)  {
          // error, should never happen
          chID[0] = *p1;
          chID[1] = char(0);
          p = &(p1[1]);
        } else  {
          *p = char(0);
          strcpy ( chID,p1 );
          *p = ':';
          p++;
        }
      }
      p1 = p;
      if (*p=='-')  p++;
      while ((*p>='0') && (*p<='9'))  p++;
      if (p==p1)  resNo = mmdb::ANY_RES;
      else  {
        c   = *p;
        *p = char(0);
        resNo = atoi(p1);
        *p = c;
      }
      strcpy ( insCode,p );
    }

  }


  int Domain::SelectDomain ( mmdb::PManager       MMDB,
                             mmdb::SELECTION_TYPE selType,
                             int                  modelNo,
                             bool         removeHydrogens )  {
  int selHnd;
    selHnd = MMDB->NewSelection();
    if (SelectDomain(selHnd,MMDB,selType,mmdb::SKEY_NEW,
                     modelNo,removeHydrogens))
         return 0;
    else return selHnd;

  }

  int Domain::SelectDomain ( int                  selHnd,
                             mmdb::PManager       MMDB,
                             mmdb::SELECTION_TYPE selType,
                             mmdb::SELECTION_KEY  selKey,
                             int                  modelNo,
                             bool         removeHydrogens )  {
  mmdb::ChainID chID;
  mmdb::ResName resName;
  mmdb::InsCode insCode;
  int           rc,resNo;
    if (dclass!=DCLASS_Ligand)  {
      rc = MMDB->SelectDomain ( selHnd,range,selType,selKey,modelNo );
      if ((!rc) && ligandRes)
        MMDB->Select ( selHnd,selType,modelNo,
                       "*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                       ligandRes,"*","*","*",mmdb::SKEY_CLR );
    } else  {
      getCID ( chID,resName,resNo,insCode );
      MMDB->Select ( selHnd,selType,modelNo,chID,resNo,insCode,
                     resNo,insCode,resName,"*","*","*",selKey );
      rc = 0;
    }
    if (rc)
      MMDB->DeleteSelection ( selHnd );
    else  {
      if (dclass!=DCLASS_Ligand)
        MMDB->SelectProperty ( selHnd,mmdb::SELPROP_Solvent,selType,
                               mmdb::SKEY_CLR );
      if (removeHydrogens && (selType==mmdb::STYPE_ATOM))
        MMDB->Select ( selHnd,selType,0,
                       "*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                       "*","*","H","*",mmdb::SKEY_CLR );
    }
    return rc;
  }


  const mmdb::realtype frame_shoulder[frameLen] =
                                     { 20.0,10.0,30.0,40.0 };

  void Domain::makeReferenceFrame()  {
  TFrame         p;
  mmdb::realtype r;
  int            i;

    if (frame[0][0]>-mmdb::MaxReal)  return;  // never recalculate

    p[0][0] = -0.5;
    p[1][0] =  0.5;
    p[2][0] =  0.0;
    p[3][0] =  0.0;

    p[0][1] = -sqrt(3.0)/6.0;
    p[1][1] =  p[0][1];
    p[2][1] = -2.0*p[0][1];
    p[3][1] =  0.0;

    p[0][2] = -sqrt(2.0/3.0)/4.0;
    p[1][2] =  p[0][2];
    p[2][2] =  p[0][2];
    p[3][2] = -3.0*p[0][2];

    for (i=0;i<frameLen;i++)  {
      r = itb[0][0]*p[i][0] + itb[0][1]*p[i][1] + itb[0][2]*p[i][2];
      frame[i][0] = frame_shoulder[i]*r + mx;
      r = itb[1][0]*p[i][0] + itb[1][1]*p[i][1] + itb[1][2]*p[i][2];
      frame[i][1] = frame_shoulder[i]*r + my;
      r = itb[2][0]*p[i][0] + itb[2][1]*p[i][1] + itb[2][2]*p[i][2];
      frame[i][2] = frame_shoulder[i]*r + mz;
    }

  }

  void Domain::getReferenceFrame ( RTFrame F, mmdb::mat44 & TM )  {
  int i,j,k;
    for (i=0;i<frameLen;i++)
      for (j=0;j<3;j++)  {
        F[i][j] = TM[j][3];
        for (k=0;k<3;k++)
          F[i][j] += TM[j][k]*frame[i][k];
      }
  }

  void Domain::setReferenceFrame ( RTFrame F )  {
  int i,j;
    for (i=0;i<frameLen;i++)
      for (j=0;j<3;j++)
        frame[i][j] = F[i][j];
  }


  bool CompareReferenceFrames ( RTFrame F1, RTFrame F2,
                                mmdb::realtype eps )  {
  int     i,j;
  bool ok;
    ok = true;
    for (i=0;(i<frameLen) && ok;i++)
      for (j=0;(j<3) && ok;j++)
        ok = fabs(F1[i][j]-F2[i][j])<eps;
    return ok;
  }


  bool ParallelFrames ( RTFrame F1, RTFrame F2, mmdb::realtype eps )  {
  mmdb::vect3    d;
  mmdb::realtype acc;
  int            i,j;
  bool           ok;

    for (j=0;j<3;j++)  {
      d[j] = 0.0;
      for (i=0;i<frameLen;i++)
        d[j] += F1[i][j] - F2[i][j];
      d[j] /= frameLen;
    }

    ok = true;
    for (i=0;(i<frameLen) && ok;i++)  {
      acc = eps*frame_shoulder[i];
      for (j=0;(j<3) && ok;j++)
        ok = (fabs(F1[i][j]-F2[i][j]-d[j])<acc);
    }

    return ok;

  }


  void Domain::getInertiaTensor ( mmdb::mat33 & itn,
                                  mmdb::mat44 & TMatrix,
                                  mmdb::realtype x,
                                  mmdb::realtype y,
                                  mmdb::realtype z )  {
  //  This function calculates the inertia mass tensor of
  // the domain transformed into position TMatrix,
  // relatively to the point of origin (x,y,z). TMatrix
  // operates in absolute coordinates.
  mmdb::mat33    itm1;
  mmdb::realtype dx,dy,dz;
  int            i,j,k;

    //  1. Apply rotation  T*I*T^{-1} to the inertia tensor
    //     in mass center

    for (i=0;i<3;i++)
      for (j=0;j<3;j++)  {
        itm1[i][j] = 0.0;
        for (k=0;k<3;k++)
          itm1[i][j] += itm[i][k]*TMatrix[j][k];
      }

    for (i=0;i<3;i++)
      for (j=i;j<3;j++)  {
        itn[i][j] = 0.0;
        for (k=0;k<3;k++)
          itn[i][j] += TMatrix[i][k]*itm1[k][j];
      }


    //  2. Add translation of mass center T*{mx,my,mz) and bring
    //     the new inertia tensor into new origin (x,y,z)

    dx = mx;
    dy = my;
    dz = mz;
    mmdb::TransformXYZ ( TMatrix,dx,dy,dz );
    dx -= x;
    dy -= y;
    dz -= z;

    itn[0][0] += weight*(dy*dy+dz*dz);
    itn[0][1] -= weight*dx*dy;
    itn[0][2] -= weight*dx*dz;
    itn[1][1] += weight*(dx*dx+dz*dz);
    itn[1][2] -= weight*dy*dz;
    itn[2][2] += weight*(dx*dx+dy*dy);
    itn[1][0]  = itn[0][1];
    itn[2][0]  = itn[0][2];
    itn[2][1]  = itn[1][2];

  }


  void Domain::calcDimensions ( mmdb::PManager MMDB )  {
  mmdb::PPAtom  atom;
  mmdb::rmatrix  A,T;
  mmdb::rvector  Eigen,Aik;
  mmdb::mat33    itm1;
  mmdb::vect3    c;
  mmdb::realtype aw, x,y,z, rx,ry,rz;
  int            selHnd,i,j,k,n;
  mmdb::SELECTION_PROPERTY sp;

    if (R>0.0)  {
      // Never recalculate domain dimensions.
      // Check whether the frame has been calculated.
      makeReferenceFrame();
      return;
    }

    nAtoms = 0;
    nRes   = 0;
    nRes0  = 0;
    mx     = 0.0;
    my     = 0.0;
    mz     = 0.0;
    for (i=0;i<3;i++)
      for (k=0;k<3;k++)
        itm[i][k] = 0.0;
    weight = 0.0;
    R      = -1.0;
    selHnd = SelectDomain ( MMDB,mmdb::STYPE_ATOM,1,false );

    if (selHnd<=0)  return;

    MMDB->GetSelIndex ( selHnd,atom,nAtoms );
    // calculate the number of "relevant" residues - used only in output
    k = MMDB->NewSelection();
    MMDB->Select ( k,mmdb::STYPE_RESIDUE,selHnd,mmdb::SKEY_OR );
    nRes = MMDB->GetSelLength ( k );
    if (dclass!=DCLASS_Ligand)  {
      if ((dclass==DCLASS_DNA) || (dclass==DCLASS_RNA))
            sp = mmdb::SELPROP_Nucleotide;
      else  sp = mmdb::SELPROP_Aminoacid;
      MMDB->SelectProperty ( k,sp,mmdb::STYPE_RESIDUE,mmdb::SKEY_AND );
      nRes0 = MMDB->GetSelLength ( k );
    } else
      nRes0 = nRes;
    MMDB->DeleteSelection ( k );


    //   In the next loop we calculate the inertia tensor
    // relatively to the *absolute* origin of coordinates
    // (0,0,0).

    // kth atom will be used for canonicalising the reference frame
    if (dclass==DCLASS_Ligand)  k = nAtoms/2;
                          else  k = nRes0/3;
    k = -mmdb::IMax(1,k);
    n = 0;
    for (i=0;i<nAtoms;i++)  {
      if (atom[i])  {
        if (!atom[i]->Ter)  {
          aw = mmdb::getMolecWeight ( atom[i]->element );
          x  = atom[i]->x;
          y  = atom[i]->y;
          z  = atom[i]->z;
          if (k<0)  {
            if ((dclass==DCLASS_Ligand) ||
                (!strcmp(atom[i]->name," CA ")))  {
              k++;
              if (k>=0)  {
                rx = x;
                ry = y;
                rz = z;
              }
            }
          }
          mx += aw*x;
          my += aw*y;
          mz += aw*z;
          weight += aw;
          itm[0][0] += aw*(y*y+z*z);
          itm[0][1] -= aw*x*y;
          itm[0][2] -= aw*x*z;
          itm[1][1] += aw*(x*x+z*z);
          itm[1][2] -= aw*y*z;
          itm[2][2] += aw*(x*x+y*y);
          n++;
        }
      }
    }

    if ((weight>0.0) && (n==1))  {

      mx /= weight;
      my /= weight;
      mz /= weight;

      // Domain made of the only atom. Find its radius.
      k = -1;
      for (i=0;(i<nAtoms) && (k<0);i++)
        if (atom[i])  k = i;
      R = mmdb::getVdWaalsRadius ( atom[k]->element );

      // calculate inertia tensor as that of a solid sphere
      // in the center of its mass
      y = 2.0/5.0*weight*R*R;
      for (i=0;i<3;i++)  {
        for (j=0;j<3;j++)  {
          itm[i][j] = 0.0;
          itb[i][j] = 0.0;
        }
        itm[i][i] = y;
        im [i]    = y;
        itb[i][i] = 1.0;
      }

      if (ncsOpNo)  //  This domain is an NCS-mate, so transform
                    // everything into position ncs_m
        mmdb::TransformXYZ ( ncs_m,mx,my,mz );


    } else if (weight>0.0)  {

      mx /= weight;
      my /= weight;
      mz /= weight;

      //   Inertia tensor itm is currently calculated relatively to
      // the absolute coordinate origin (0,0,0). Now bring it into
      // the mass center (mx,my,mz).

      itm[0][0] -= weight*(my*my+mz*mz);
      itm[0][1] += weight*mx*my;
      itm[0][2] += weight*mx*mz;
      itm[1][1] -= weight*(mx*mx+mz*mz);
      itm[1][2] += weight*my*mz;
      itm[2][2] -= weight*(mx*mx+my*my);
      itm[1][0]  = itm[0][1];
      itm[2][0]  = itm[0][2];
      itm[2][1]  = itm[1][2];

      //   Calculate the maximum radius of the domain
      R = 0.0;
      for (i=0;i<nAtoms;i++)
        if (atom[i])  {
          if (!atom[i]->Ter)
            R = mmdb::RMax ( R,atom[i]->GetDist2(mx,my,mz) );
        }
      R = sqrt ( R );

      if (ncsOpNo)  {

        //  This domain is an NCS-mate, so transform everything
        //  into position ncs_m

        x = mx;
        y = my;
        z = mz;
        mmdb::TransformXYZ ( ncs_m,x,y,z );
        mmdb::TransformXYZ ( ncs_m,rx,ry,rz );

        //  Transform the inertia tensor into new position
        //  relatively to the new mass center

        getInertiaTensor ( itm1,ncs_m,x,y,z );
        for (i=0;i<3;i++)
          for (j=0;j<3;j++)
            itm[i][j] = itm1[i][j];

        mx = x;
        my = y;
        mz = z;

      }

      //  Calculate the principal moments of inertia relatively
      // to the mass center.

      mmdb::GetMatrixMemory ( A,3,3,1,1 );
      mmdb::GetMatrixMemory ( T,3,3,1,1 );
      mmdb::GetVectorMemory ( Eigen,3,1 );
      mmdb::GetVectorMemory ( Aik  ,3,1 );

      for (i=1;i<=3;i++)
        for (j=1;j<=3;j++)
          A[i][j] = itm[i-1][j-1];

      mmdb::math::Jacobi ( 3,A,T,Eigen,Aik,i );

      for (i=0;i<3;i++)  {
        for (k=0;k<3;k++)
          itb[i][k] = T[i+1][k+1];
        im[i] = weight*Eigen[i+1];
      }

      // canonicalise the basis
      aw = mmdb::MaxReal;
      k  = -1;
      for (i=0;i<3;i++)  {
        c[i] = rx*itb[0][i] + ry*itb[1][i] + rz*itb[2][i];
        if (fabs(c[i])<aw)  {
          aw = fabs(c[i]);
          k  = i;
        }
      }

      if (aw>0.3)  k = -1;  // all cosines are well defined

      // choose axes along largest positive cosines
      for (i=0;i<3;i++)
        if ((i!=k) && (c[i]<0.0))  {
          itb[0][i] = -itb[0][i];
          itb[1][i] = -itb[1][i];
          itb[2][i] = -itb[2][i];
        }

      if (k>=0)  {
        // One cosine is ill-defined.
        // Choose the remained axis in the right-hand screw direction
        // (cross-product of the previously defined axes)
        if (k==0)  {
          j = 1;  n = 2;
        } else if (k==1)  {
          j = 0;  n = 2;
        } else  {
          j = 0;  n = 1;
        }
        x = itb[1][j]*itb[2][n] - itb[2][j]*itb[1][n]; // yz-zy
        y = itb[2][j]*itb[0][n] - itb[0][j]*itb[2][n]; // zx-xz
        z = itb[0][j]*itb[1][n] - itb[1][j]*itb[0][n]; // xy-yx

        c[k] = x*itb[0][k] + y*itb[1][k] + z*itb[2][k];
        if (c[k]<0.0)  {  // c[k] is either +1 or -1
          itb[0][k] = -itb[0][k];
          itb[1][k] = -itb[1][k];
          itb[2][k] = -itb[2][k];
        }
      }

      mmdb::FreeVectorMemory ( Aik  ,1 );
      mmdb::FreeVectorMemory ( Eigen,1 );
      mmdb::FreeMatrixMemory ( T,3,1,1 );
      mmdb::FreeMatrixMemory ( A,3,1,1 );

    }

    MMDB->DeleteSelection ( selHnd );

    makeReferenceFrame();

  }

  #define  contThresh  5.0

  bool Domain::isContact ( PDomain D )  {
  mmdb::realtype dx,dy,dz,A;
    if ((R>0.0) && (D->R>0.0))  {
      dx = mx - D->mx;
      dy = my - D->my;
      dz = mz - D->mz;
      A  = dx*dx + dy*dy + dz*dz;
      dx = R + D->R + contThresh;  // plus water radii
      return (A<=dx*dx);
    } else
      return false;
  }

  bool Domain::isContact ( PDomain D, mmdb::mat44 & TM )  {

  mmdb::realtype dx,dy,dz,A;
    if ((R>0.0) && (D->R>0.0))  {
      dx = mx - (TM[0][0]*D->mx+TM[0][1]*D->my+TM[0][2]*D->mz+TM[0][3]);
      dy = my - (TM[1][0]*D->mx+TM[1][1]*D->my+TM[1][2]*D->mz+TM[1][3]);
      dz = mz - (TM[2][0]*D->mx+TM[2][1]*D->my+TM[2][2]*D->mz+TM[2][3]);
      A  = dx*dx + dy*dy + dz*dz;
      dx = R + D->R + contThresh;  // plus water radii
      return (A<=dx*dx);
    } else
      return false;
  }

  bool Domain::isContact ( mmdb::mat44 & T, PDomain D, mmdb::mat44 & TM )  {
  mmdb::realtype dx,dy,dz,A;

    if ((R>0.0) && (D->R>0.0))  {

      dx = T[0][0]*mx + T[0][1]*my + T[0][2]*mz + T[0][3];
      dy = T[1][0]*mx + T[1][1]*my + T[1][2]*mz + T[1][3];
      dz = T[2][0]*mx + T[2][1]*my + T[2][2]*mz + T[2][3];

      dx -= TM[0][0]*D->mx + TM[0][1]*D->my + TM[0][2]*D->mz + TM[0][3];
      dy -= TM[1][0]*D->mx + TM[1][1]*D->my + TM[1][2]*D->mz + TM[1][3];
      dz -= TM[2][0]*D->mx + TM[2][1]*D->my + TM[2][2]*D->mz + TM[2][3];

      A  = dx*dx + dy*dy + dz*dz;
      dx = R + D->R + contThresh;  // plus water radii

      return (A<=dx*dx);

    } else
      return false;

  }


  #define  cont_brick_margin 15.0
  #define  cont_brick_size    9.0
  #define  max_contact2      81.0

  void  Domain::MakeContactBricks ( mmdb::PManager MMDB,
                                       mmdb::mat44 & TMatrix,
                                       int & selHnd )  {
  //  makes bricks for 'this' domain in location TMatrix
  mmdb::PPAtom  A;
  mmdb::rvector x,y,z;
  mmdb::mat44   tm;
  int           i,nAtoms;
  bool          unitTM;

    MMDB->RemoveBricks();
    if (!selHnd)
      selHnd = SelectDomain ( MMDB,mmdb::STYPE_ATOM,1,false );

    MMDB->GetSelIndex ( selHnd,A,nAtoms );
    if (A && (nAtoms>0))  {
      mmdb::Mat4Mult ( tm,TMatrix,ncs_m );
      unitTM = mmdb::isMat4Unit ( tm,0.00001,false );
      if (!unitTM)  {
        mmdb::GetVectorMemory ( x,nAtoms,0 );
        mmdb::GetVectorMemory ( y,nAtoms,0 );
        mmdb::GetVectorMemory ( z,nAtoms,0 );
        for (i=0;i<nAtoms;i++)  {
          x[i] = A[i]->x;
          y[i] = A[i]->y;
          z[i] = A[i]->z;
          A[i]->Transform ( tm );
        }
      }
      MMDB->MakeBricks ( A,nAtoms,cont_brick_margin,cont_brick_size );
      if (!unitTM)  {
        for (i=0;i<nAtoms;i++)  {
          A[i]->x = x[i];
          A[i]->y = y[i];
          A[i]->z = z[i];
        }
        mmdb::FreeVectorMemory ( x,0 );
        mmdb::FreeVectorMemory ( y,0 );
        mmdb::FreeVectorMemory ( z,0 );
      }
    }

  }

  #define  specPosThresh2  0.0225

  bool Domain::isContact ( mmdb::mat44 & TM, PDomain D,
                                 mmdb::PManager MMDB,
                                 int & selHnd )  {
  //  This function returns true if domain D is in potential
  // contact with 'this' domain in orientation TM. Unlike other
  // isContact functions, this one checks contacts more thoroughly,
  // on atom level. Atoms of 'this' domain in the original
  // orientation should be bricked in MMDB outside this function.
  // The calling function is also responsible for removing bricks.
  // The function also applies all necessary NCS transformations
  // to both this domain and D.
  //  selHnd is used to minimize the number of selection in
  // external loop. Initially it should be set 0 to initiate
  // selection in this function. The calling function is
  // then responsible for deleting all selections.
  mmdb::PPAtom   A;
  mmdb::PBrick   B;
  mmdb::mat44    tmatrix;
  mmdb::realtype x,y,z,radius, x0,y0,z0,minDist;
  int            nAtoms,i,j,nx,ny,nz, ib,jb,kb;
  bool           contact,unitTM,specPos;

    if ((R<=0.0) || (D->R<=0.0))  return false;  // no atoms at all

    //  first make rough estimate

    x = D->mx - (TM[0][0]*mx+TM[0][1]*my+TM[0][2]*mz+TM[0][3]);
    y = D->my - (TM[1][0]*mx+TM[1][1]*my+TM[1][2]*mz+TM[1][3]);
    z = D->mz - (TM[2][0]*mx+TM[2][1]*my+TM[2][2]*mz+TM[2][3]);
    radius = x*x + y*y + z*z;  // square distance between mass centers

    x = R + D->R + contThresh;  // maximal contact distance between
                                // mass centers

    if (radius>x*x)
      return false;  // no contact for sure

    specPos = (radius<specPosThresh2);  // possible special position

    //  check contacts on atom level

    if (!selHnd)
      selHnd = D->SelectDomain ( MMDB,mmdb::STYPE_ATOM,1,false );

    MMDB->GetSelIndex ( selHnd,A,nAtoms );
    contact = false;

    if (A && (nAtoms>0))  {

      mmdb::Mat4Mult ( tmatrix,TM,ncs_m );
      unitTM  = mmdb::isMat4Unit(tmatrix,0.00001,false);

      for (i=0;(i<nAtoms) && (!contact);i++)  {
        if (D->ncsOpNo>0)
          A[i]->TransformCopy ( D->ncs_m,x,y,z );
        else  {
          x = A[i]->x;
          y = A[i]->y;
          z = A[i]->z;
        }
        // {x,y,z} is absolute location of ith atom of domain D
        MMDB->GetBrickCoor ( x,y,z, nx,ny,nz );
        // {nx,ny,nz} is brick containing ith atom of domain D
        if (nx>=0)  {
          // check +/- 1 brick surroundings
          for (ib=nx-1;(ib<=nx+1) && (!contact);ib++)
            for (jb=ny-1;(jb<=ny+1) && (!contact);jb++)
              for (kb=nz-1;(kb<=nz+1) && (!contact);kb++)  {
                B = MMDB->GetBrick ( ib,jb,kb );
                if (B)
                  for (j=0;(j<B->nAtoms) && (!contact);j++)  {
                    if (unitTM)  {
                      x0 = B->atom[j]->x;
                      y0 = B->atom[j]->y;
                      z0 = B->atom[j]->z;
                    } else
                      B->atom[j]->TransformCopy ( tmatrix,x0,y0,z0 );
                    // {x0,y0,z0} is absolute location of atom
                    // in 'this' domain
                    x0 -= x;
                    x0 *= x0;
                    if (x0<=max_contact2)  {
                      y0 -= y;
                      y0 *= y0;
                      if (y0<=max_contact2)  {
                        z0 -= z;
                        radius  = x0 + y0 + z0*z0;
                        if (radius<=max_contact2)  contact = true;
                      }
                    }
                  }
              }
        }
      }

      if (contact && specPos)  {
        // Check that the structures do not do a perfect overlap
        // in crystallographic special position. If special
        // position is found (minimal distance less than 0.15 A)
        // then no contacts is assumed
        contact = false;
        for (i=0;(i<nAtoms) && (!contact);i++)  {
          if (D->ncsOpNo>0)
            A[i]->TransformCopy ( D->ncs_m,x,y,z );
          else  {
            x = A[i]->x;
            y = A[i]->y;
            z = A[i]->z;
          }
          minDist = mmdb::MaxReal;
          // {x,y,z} is absolute location of ith atom of domain D
          MMDB->GetBrickCoor ( x,y,z, nx,ny,nz );
          // {nx,ny,nz} is brick containing ith atom of domain D
          if (nx>=0)  {
            // check +/- 1 brick surroundings
            for (ib=nx-1;ib<=nx+1;ib++)
              for (jb=ny-1;jb<=ny+1;jb++)
                for (kb=nz-1;kb<=nz+1;kb++)  {
                  B = MMDB->GetBrick ( ib,jb,kb );
                  if (B)
                    for (j=0;j<B->nAtoms;j++)  {
                      if (unitTM)  {
                        x0 = B->atom[j]->x;
                        y0 = B->atom[j]->y;
                        z0 = B->atom[j]->z;
                      } else
                        B->atom[j]->TransformCopy ( tmatrix,x0,y0,z0 );
                      // {x0,y0,z0} is absolute location of atom
                      // in 'this' domain
                      x0 -= x;
                      x0 *= x0;
                      if (x0<=specPosThresh2)  {
                        y0 -= y;
                        y0 *= y0;
                        if (y0<=specPosThresh2)  {
                          z0 -= z;
                          radius  = x0 + y0 + z0*z0;
                          if (radius<minDist)  minDist = radius;
                        }
                      }
                    }
                }
          }
          contact = (minDist>specPosThresh2);
        }
      }

    }

    return contact;

  }


  PROSURF_RC Domain::calcSurface ( mmdb::PManager MMDB,
                                   PProSurf       proSurf,
                                   PMolRefIndex   molRef )  {
  SolvEnergy      SolvEnergy;
  mmdb::PPResidue res;
  mmdb::PPAtom    atom;
  mmdb::rvector   SAS,atomSE;
  mmdb::ivector   index;
  int             selHnd,selHndR;
  int             i,k;
  PROSURF_RC      rc;

    if (ncsOpNo)  return PROSURF_Stop1;

    // select domain without hydrogens
    selHnd     = SelectDomain ( MMDB,mmdb::STYPE_ATOM,1,true );
    nSurfAtoms = 0;
    nSurfRes   = 0;
    surfArea   = 0.0;
    DeltaG     = 0.0;     // relative solvation energy

    if (selHnd<=0)  return PROSURF_Stop2;

    mmdb::FreeVectorMemory  ( asa   ,0 );
    mmdb::FreeVectorMemory  ( solvEn,0 );

    MMDB->GetSelIndex ( selHnd,atom,nAtoms );
    mmdb::GetVectorMemory   ( SAS   ,nAtoms,0    );
    mmdb::GetVectorMemory   ( atomSE,nAtoms,0    );

    rc = proSurf->calcSurface ( atom,nAtoms,molRef,SAS );

    if (rc==PROSURF_Ok)  {

      selHndSurf = proSurf->selHndSurf1;
      nSurfAtoms = MMDB->GetSelLength ( selHndSurf );

      // calculate the number of surface residues
      selHndR    = MMDB->NewSelection();
      MMDB->Select ( selHndR,mmdb::STYPE_RESIDUE,
                     selHndSurf,mmdb::SKEY_NEW );
      nSurfRes   = MMDB->GetSelLength ( selHndR );
      MMDB->DeleteSelection ( selHndR );

      // calculate solvation energy gain
      surfArea   = proSurf->surfArea1;
      DeltaG     = SolvEnergy.calcSurfSolvEnergy ( atom,SAS,nAtoms,
                                                   atomSE,molRef );

      // calculate residue asa and solvation energy gain
      selHndR    = MMDB->NewSelection();
      MMDB->Select ( selHndR,mmdb::STYPE_RESIDUE,selHnd,mmdb::SKEY_NEW );
      MMDB->GetSelIndex ( selHndR,res,nRes );

      mmdb::GetVectorMemory ( asa   ,nRes,0 );
      mmdb::GetVectorMemory ( solvEn,nRes,0 );
      mmdb::GetVectorMemory ( index ,nRes,0 );
      for (i=0;i<nRes;i++)  {
        asa   [i] = 0.0;
        solvEn[i] = 0.0;
        index [i] = res[i]->index;
        res   [i]->index = i;
      }

      for (i=0;i<nAtoms;i++)  {
        k = atom[i]->GetResidue()->index;
        asa   [k] += SAS[i];
        solvEn[k] += atomSE[i];
      }

      for (i=0;i<nRes;i++)
        res[i]->index = index[i];

      mmdb::FreeVectorMemory ( index,0 );
      MMDB->DeleteSelection ( selHndR );

    }

    mmdb::FreeVectorMemory ( SAS   ,0 );
    mmdb::FreeVectorMemory ( atomSE,0 );
    MMDB->DeleteSelection  ( selHnd   );

    return rc;

  }

  mmdb::pstr Domain::getDomainRange ( mmdb::pstr S, int fieldLen )  {
  int k;

    if (range)  {
      if (ncsOpNo>0)  {
        if (dclass==DCLASS_Ligand)
          sprintf ( S,"%i\"%s",ncsOpNo,range );
        else  {
          k = strlen(range)-1;
          if (range[k]!=':')  k++;
          if (fieldLen>0)  sprintf ( S,"%4i\"",ncsOpNo );
                     else  sprintf ( S,"%i\"",ncsOpNo );
          strncat ( S,range,k );
        }
      } else  {
        if (dclass==DCLASS_Ligand)
          strcpy ( S,range );
        else  {
          k = strlen(range)-1;
          if (range[k]!=':')  k++;
          if (fieldLen>0)  strcpy  ( S,"     " );
                     else  S[0] = char(0);
          strncat ( S,range,k );
        }
      }
    } else
      strcpy ( S,"<no range>" );

    k = strlen(S);
    while (k<fieldLen)
      S[k++] = ' ';
    S[k] = char(0);

    return S;

  }

  mmdb::pstr Domain::getFileSafeName ( mmdb::pstr S, int fieldLen )  {
  int i,j;
    if (range)  {
      getDomainRange ( S,fieldLen );
      i = 0;
      j = 0;
      while (S[i])  {
        if ((S[i]!='[') && (S[i]!=' '))  {
          if ((S[i]==']') || (S[i]==':'))
            S[j++] = '_';
          else if (j<i)
            S[j++] = S[i];
        }
        i++;
      }
      S[j] = char(0);
    } else
      strcpy ( S,"unk" );
    return S;
  }


  mmdb::pstr Domain::getDomainClass ( mmdb::pstr S )  {
    switch (dclass)  {
      case DCLASS_Protein : strcpy ( S,"Protein" );  break;
      case DCLASS_DNA     : strcpy ( S,"DNA"     );  break;
      case DCLASS_RNA     : strcpy ( S,"RNA"     );  break;
      case DCLASS_Ligand  : strcpy ( S,"Ligand " );  break;
      default : strcpy ( S,"<error>" );
    }
    return S;
  }


  mmdb::pstr Domain::getLigandRange ( mmdb::pstr S )  {
  int i,k;
    if (range)  {
      if (dclass==DCLASS_Ligand)  {
        i = 0;
        k = 1;
        while (range[k] && (range[k]!=']'))
          S[i++] = range[k++];
        S[i] = char(0);
      } else
        strcpy ( S,"<not lig>" );
    } else
      strcpy ( S,"<no range>" );
    return S;
  }

  mmdb::pstr Domain::getChainID ( mmdb::pstr S )  {
  mmdb::pstr p;
  int  k;
    if (range)  {
      if (dclass==DCLASS_Ligand)  {
        p = mmdb::FirstOccurence ( range,']' );
        if (p)  {
          strcpy ( S,&(p[1]) );
          p = mmdb::FirstOccurence ( S,':' );
          if (p)  *p = char(0);
        } else
          strcpy ( S,"<wrong lig id>" );
      } else if (range[0]!='-')  {
        strcpy ( S,range );
        k = strlen(S)-1;
        if (S[k]==':')  S[k] = char(0);
      } else
        strcpy ( S,"-" );
    } else
      strcpy ( S,"<no chain id>" );
    return S;
  }


  mmdb::pstr Domain::getDomainID ( mmdb::pstr S )  {
  int  k;
    if (range)  {
      if (range[0]!='-')  {
        strcpy ( S,range );
        k = strlen(S)-1;
        if (S[k]==':')  S[k] = char(0);
      } else
        strcpy ( S,"-" );
    } else
      strcpy ( S,"<no range>" );
    return S;
  }


  mmdb::pstr Domain::getDomainRange_html ( mmdb::pstr S,
                                           bool       separate )  {
  int  k;
    if (range)  {
      if (ncsOpNo>0)  {
        if (separate) sprintf ( S,"<sup>&nbsp;%i</sup>",ncsOpNo );
                 else sprintf ( S,"<sup>%i</sup>",ncsOpNo );
      } else
        S[0] = char(0);
      if (dclass==DCLASS_Ligand)
        strcat ( S,range );
      else if (range[0]!='-')  {
        strcat ( S,range );
        k = strlen(S)-1;
        if (S[k]==':')  S[k] = char(0);
      } else
        strcat ( S,"{-}" );
    } else
      strcpy ( S,"<no range>" );
    return S;
  }

  void Domain::CopyParentData ( PPDomain domain )  {
  PDomain     D;
  mmdb::mat33 tb;
  int         i,j,k;

    if (ncsOpNo>0)  {

      D = domain[ncsParent];
      if (D)  {

        mmdb::CreateCopy ( range ,D->range  );
        mmdb::CreateCopy ( typeId,D->typeId );
        selHndSurf = D->selHndSurf;
        nAtoms     = D->nAtoms;
        nRes       = D->nRes;
        nRes0      = D->nRes0;
        nSurfAtoms = D->nSurfAtoms;
        nSurfRes   = D->nSurfRes;
        type       = D->type;
        dclass     = D->dclass;
        symNumber  = D->symNumber;
        assemble   = D->assemble;
        surfArea   = D->surfArea;
        DeltaG     = D->DeltaG;
        mx         = D->mx;
        my         = D->my;
        mz         = D->mz;

        mmdb::TransformXYZ ( ncs_m,mx,my,mz );
        for (i=0;i<frameLen;i++)  {
          for (j=0;j<3;j++)
            frame[i][j] = D->frame[i][j];
          mmdb::TransformXYZ ( ncs_m,frame[i][0],frame[i][1],frame[i][2] );
        }
        weight     = D->weight;
        R          = D->R;
        asa        = D->asa;     // copy pointer only for NCS-mate
        solvEn     = D->solvEn;  // copy pointer only for NCS-mate

        D->getInertiaTensor ( itm,ncs_m,mx,my,mz );

        for (i=0;i<3;i++)  {
          im[i] = D->im[i];
          for (j=0;j<3;j++)  {
            tb[i][j] = 0.0;
            for (k=0;k<3;k++)
              tb[i][j] += D->itb[i][k]*ncs_m[j][k];
          }
        }
        for (i=0;i<3;i++)
          for (j=0;j<3;j++)  {
            itb[i][j] = 0.0;
            for (k=0;k<3;k++)
              itb[i][j] += ncs_m[i][k]*tb[k][j];
          }

        mmdb::CreateCopy ( ligandRes,D->ligandRes );

      }

    }

  }


  void Domain::write ( mmdb::io::RFile f )  {
  int        i,j,dc;
  mmdb::byte Version=4;

    f.WriteByte   ( &Version    );

    f.WriteInt    ( &ncsOpNo    );
    f.WriteInt    ( &ncsParent  );
    f.CreateWrite ( typeId      );

    if (!ncsOpNo)  {

      f.CreateWrite ( range     );
      f.WriteInt  ( &selHndSurf );
      f.WriteInt  ( &nAtoms     );
      f.WriteInt  ( &nRes       );
      f.WriteInt  ( &nRes0      );
      f.WriteInt  ( &nSurfAtoms );
      f.WriteInt  ( &nSurfRes   );
      f.WriteInt  ( &type       );
      dc = dclass;
      f.WriteInt  ( &dc         );
      f.WriteInt  ( &symNumber  );
      f.WriteBool ( &assemble   );
      if (Version>2)
        f.WriteBool ( &agent );
      f.WriteReal ( &surfArea   );
      f.WriteReal ( &DeltaG     );
      f.WriteReal ( &mx         );
      f.WriteReal ( &my         );
      f.WriteReal ( &mz         );
      f.WriteReal ( &weight     );
      f.WriteReal ( &R          );

      if (asa) {
        j = 1;
        f.WriteInt ( &j );
        for (i=0;i<nRes;i++)  {
          f.WriteReal ( &(asa   [i]) );
          f.WriteReal ( &(solvEn[i]) );
        }
      } else  {
        j = 0;
        f.WriteInt ( &j );
      }

      for (i=0;i<frameLen;i++)
        for (j=0;j<3;j++)
          f.WriteReal ( &(frame[i][j]) );

      for (i=0;i<3;i++)  {
        for (j=0;j<3;j++)  {
          f.WriteReal ( &(itb[i][j]) );
          f.WriteReal ( &(itm[i][j]) );
        }
        f.WriteReal ( &(im[i]) );
      }

      f.CreateWrite ( ligandRes );

    } else  {

      for (i=0;i<3;i++)
        for (j=0;j<4;j++)
          f.WriteReal ( &(ncs_m[i][j]) );

    }

  }

  void Domain::read ( mmdb::io::RFile f )  {
  int        i,j,dc;
  mmdb::byte Version;

    f.ReadByte ( &Version   );

    f.ReadInt  ( &ncsOpNo   );
    f.ReadInt  ( &ncsParent );
    if (Version>3)
      f.CreateRead ( typeId  );
    else mmdb::CreateCopy ( typeId,range );

    if (!ncsOpNo)  {

      f.CreateRead ( range     );
      f.ReadInt  ( &selHndSurf );
      f.ReadInt  ( &nAtoms     );
      f.ReadInt  ( &nRes       );
      f.ReadInt  ( &nRes0      );
      f.ReadInt  ( &nSurfAtoms );
      f.ReadInt  ( &nSurfRes   );
      f.ReadInt  ( &type       );
      f.ReadInt  ( &dc         );
      dclass = (DOMAIN_CLASS)dc;
      f.ReadInt  ( &symNumber  );
      if (Version>1)  f.ReadBool ( &assemble );
                else  assemble = true;
      if (Version>2)  f.ReadBool ( &agent );
                else  agent    = false;
      f.ReadReal ( &surfArea   );
      f.ReadReal ( &DeltaG     );
      f.ReadReal ( &mx         );
      f.ReadReal ( &my         );
      f.ReadReal ( &mz         );
      f.ReadReal ( &weight     );
      f.ReadReal ( &R          );

      mmdb::Mat4Init ( ncs_m );

      mmdb::FreeVectorMemory ( asa   ,0 );
      mmdb::FreeVectorMemory ( solvEn,0 );

      f.ReadInt ( &j );
      if (j>0)  {
        mmdb::GetVectorMemory ( asa   ,nRes,0 );
        mmdb::GetVectorMemory ( solvEn,nRes,0 );
        for (i=0;i<nRes;i++)  {
          f.ReadReal ( &(asa   [i]) );
          f.ReadReal ( &(solvEn[i]) );
        }
      }

      for (i=0;i<frameLen;i++)
        for (j=0;j<3;j++)
          f.ReadReal ( &(frame[i][j]) );

      for (i=0;i<3;i++)  {
        for (j=0;j<3;j++)  {
          f.ReadReal ( &(itb[i][j]) );
          f.ReadReal ( &(itm[i][j]) );
        }
        f.ReadReal ( &(im[i]) );
      }

      f.CreateRead ( ligandRes );

    } else  {

      for (i=0;i<3;i++)
        for (j=0;j<4;j++)
          f.ReadReal ( &(ncs_m[i][j]) );
      ncs_m[3][0] = 0.0;
      ncs_m[3][1] = 0.0;
      ncs_m[3][2] = 0.0;
      ncs_m[3][3] = 1.0;

    }

  }

  MakeStreamFunctions(Domain)


}  // namespace pisa
