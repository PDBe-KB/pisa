// $Id: pisa_engine.cpp $
// =================================================================
//
//    29.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  asm_engine <implementation>
//       ~~~~~~~~~
//  **** Project :  Protein interfaces
//       ~~~~~~~~~
//  **** Classes :  pisa::Assembler
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2004-2013
//
// =================================================================
//

#include <math.h>
#include <string.h>

#include "pisa_engine.h"
#include "pisa_types.h"

namespace pisa  {

  // =========================  Assembler  ==========================

  #define RSTAT_Ok        0
  #define RSTAT_Finish    1
  #define RSTAT_QuitLoop  2
  #define RSTAT_TimeLimit 3

  static int maxNofETypes = 260;

  void SetMaxNofETypes ( int maxNETypes )  {
    maxNofETypes = maxNETypes;
  }

  Assembler::Assembler()  {
    InitAssembler();
    startCPUclock = 0;
  }

  Assembler::~Assembler()  {
    DeleteMonomers     ( M,nMonomers,nMonAlloc );
    DeletePAGraph      ();
    DeleteMultimers    ( U,nMultimers,nMultAlloc );
    DeleteComplex      ();
    DeleteRecStack     ();
    DeleteMultimerSets ( multSet,nMultSets,nMultSetAlloc );
    RemoveBricks       ();
  }

  void Assembler::InitAssembler()  {

    interface     = NULL;  // allocated and resorted
    M             = NULL;  // monomers
    ucv           = NULL;  // unit cell translations
    eType         = NULL;  // interface enumeration types
    parMon        = NULL;  // matrix of parallel monomers
    U             = NULL;  // vector of multimers
    recStack      = NULL;  // recursion stack
    multSet       = NULL;  // sets of resulting multimers
    Complex       = NULL;  // multimer for single-complex analysis
    ligKey        = LIGANDS_Auto; // ligand processing key
    retcode       = ASSMB_Void;   // return code
    nSymOps       = 0;     // number of symmetry operations
    nInterfaces   = 0;     // calculated number of interfaces
    nInterfaces0  = 0;     // original number of interfaces
    nMonomers     = 0;     // number of monomers in unit cell
    nCellLayers   = 3;     // number of unit cell layers
    nETypes       = 0;     // number of interface enumeration types
    nMonAlloc     = 0;     // allocation monomer number
    nMultimers    = 0;     // current number of multimers
    nMultSets     = 0;     // number of resulting sets
    nMultAlloc    = 0;     // allocation multimer number
    nRSAlloc      = 0;     // alloc. num. for disposal of recursion stack
    nMultSetAlloc = 0;     // number of allocated multimer sets
    Crystal       = true;  // true if crystal is given
    fullSearch    = false; // true if full search is requested

    //   enumeration parameters to control the enumeration length
    //   for large structures
    dissThreshold  = 0.0;  // dissociation threshold
    maxStableSize  = 0;
    noImproveCnt   = 0;
    noImproveLimit = 10;   // enumeration limit for large structures
    runStatus      = RSTAT_Ok;
    timeLimit      = 3*3600;  // // CPU time limit in seconds

    //   the ligands are fixed (permanently attached to macromolecules)
    //   if their size is less than fixLigandSize or if there are
    //   more than fixLigandNum such ligands in ASU
    fixLigandSize  = 80;  // maximal size of ligands to fix compulsory
    fixLigandNum   = 2;   // maximal number of non-fixed ligands

    MB         = NULL;

    nx_bricks  = 0;
    ny_bricks  = 0;
    nz_bricks  = 0;

    brick_size = 100.0;
    x1_brick   = 0.0;
    x2_brick   = 0.0;
    y1_brick   = 0.0;
    y2_brick   = 0.0;
    z1_brick   = 0.0;
    z2_brick   = 0.0;

  }


  void Assembler::DeletePAGraph()  {
  int i,j;

    if (interface)  {
      delete[] interface;
      interface = NULL;
    }
    nInterfaces  = 0;
    nInterfaces0 = 0;

    if (ucv)  {
      for (i=-nCellLayers;i<=nCellLayers;i++)
        if (ucv[i])  {
          for (j=-nCellLayers;j<=nCellLayers;j++)
            if (ucv[i][j])  {
              ucv[i][j] -= nCellLayers;
              delete[] ucv[i][j];
            }
          ucv[i] -= nCellLayers;
          delete[] ucv[i];
        }
      ucv -= nCellLayers;
      delete[] ucv;
      ucv = NULL;
    }

    mmdb::FreeVectorMemory ( eType,0 );
    mmdb::FreeMatrixMemory ( parMon,nMonAlloc,0,0 );
    nMonAlloc = 0;

  }

  void Assembler::DeleteMonomers ( PPMonomer & Mon,
                                   int & nM, int & nMA )  {
  int i;

    if (Mon)  {
      for (i=0;i<nMA;i++)
        if (Mon[i])  delete Mon[i];
      delete[] Mon;
      Mon = NULL;
    }
    nM  = 0;
    nMA = 0;

  }

  void Assembler::DeleteMultimers ( PPMultimer & Mult,
                                    int & nM, int & nMA )  {
  int i;
    if (Mult)  {
      for (i=0;i<nMA;i++)
        if (Mult[i])  delete Mult[i];
      delete[] Mult;
      Mult = NULL;
    }
    nM  = 0;
    nMA = 0;
  }

  void Assembler::DeleteComplex ()  {
    if (Complex)  {
      delete Complex;
      Complex = NULL;
    }
  }

  void Assembler::DeleteRecStack()  {
  int i;
    if (recStack)  {
      for (i=0;i<nRSAlloc;i++)
        if (recStack[i])  delete recStack[i];
      delete[] recStack;
      recStack = NULL;
    }
    nRSAlloc = 0;
  }

  void Assembler::DeleteMultimerSets ( PPMultimerSet & mSet,
                                        int & nM, int & nMA )  {
  int i;
    if (mSet)  {
      for (i=0;i<nMA;i++)
        if (mSet[i])  delete mSet[i];
      delete[] mSet;
      mSet = NULL;
    }
    nM  = 0;
    nMA = 0;
  }


  void Assembler::Reset()  {
    retcode = ASSMB_Void;
    DeleteMonomers     ( M,nMonomers,nMonAlloc );
    DeletePAGraph      ();
    DeleteMultimers    ( U,nMultimers,nMultAlloc );
    DeleteComplex      ();
    DeleteRecStack     ();
    DeleteMultimerSets ( multSet,nMultSets,nMultSetAlloc );
  }

  int Assembler::GetRC()  {
    return retcode;
  }

  void Assembler::SetDissThreshold ( mmdb::realtype threshold )  {
    dissThreshold = threshold;
  }

  void Assembler::SetLigandKey ( int ligandKey )  {
    ligKey = ligandKey;
  }

  void Assembler::SetFullSearch ( bool fSearch )  {
    fullSearch = fSearch;
  }

  ASSMB_RC Assembler::MakeInterfaces ( PInterfaces PI,
                                       PPDomain  Domain )  {
  PInterface Intf;
  int         i,j;
  bool    getIntf;

    nInterfaces0 = PI->getNofInterfaces();
    if (nInterfaces0<=0)  return  ASSMB_noInterfaces;

    //   Make a local copy of interface properties vector, ordered by
    // decreasing the stabilization energy. All assembly calculations
    // are done with thus resorted interface vector, but the original
    // interface numbering is restored when results are being copied
    // into resulting structures (see CMolecule::_copy(..)).

    interface   = new PInterface[nInterfaces0];
    nInterfaces = 0;
    j           = nInterfaces0;
    for (i=0;i<nInterfaces0;i++)  {
      Intf = PI->getInterface(i);
      Intf->stabEn = Intf->getStabEn();
                    // function call is needed for fitting procedure
      getIntf = (Domain[Intf->domain1]->assemble &&
                 Domain[Intf->domain2]->assemble);
      if ((!Crystal) && getIntf)
        getIntf = (Intf->symOpNo<=0);
      if (getIntf)  interface[nInterfaces++] = Intf;
              else  interface[--j]           = Intf;
    }

    if (nInterfaces<=0)  return  ASSMB_noInterfaces;

    for (i=0;i<nInterfaces;i++)
      for (j=i+1;j<nInterfaces;j++)
        if (interface[j]->stabEn<interface[i]->stabEn)  {
          Intf         = interface[i];
          interface[i] = interface[j];
          interface[j] = Intf;
        }

    return ASSMB_Ok;

  }


  bool Assembler::traceMonomerChain ( PMonomer M1, int engType,
                                      int i1, int j1, int k1,
                                      PMonRef L, int cnt )  {
  PMonomer M2;
  int      i,j, i2,j2,k2;
  bool     B;

    if (!M1->L[0])  return false;  // no interfaces, nothing to trace

    M1->id2 = 2;  // signal "traced"

    B = false;    // will be true if engaged and fixed interfaces
                  // lead to monomer marked "1"

    // loop over all engaged and fixed interfaces made by monomer M1
    for (i=0;(i<nInterfaces) && (!B);i++)
      if (M1->L[0][i].M &&
          ((eType[i]==engType) || interface[i]->fixedLigand))  {
        // check whether any of ith interfaces leads to monomer
        // marked "1"
        j = 0;
        while ((j<_max_n_int) && (!B))
          if (M1->L[j])  {
            M2 = M1->L[j][i].M;
            if (M2)  {
              i2 = i1 + M1->L[j][i].i;
              j2 = j1 + M1->L[j][i].j;
              k2 = k1 + M1->L[j][i].k;
              if (M2->id2==1)  {
                // the mark is there, return "true" only if unit cell
                // is right
                B = (i2==L->i) && (j2==L->j) && (k2==L->k);
              } else if ((!M2->id2) && (cnt<3))  {
                B = traceMonomerChain ( M2,engType,i2,j2,k2,L,cnt+1 );
              }
            }
            j++;
          } else
            j = _max_n_int;
      }

    M1->id2 = 0;  // restore

    return B;

  }


  void Assembler::CalcEnumerationTypes()  {
  //
  //   This function calculates the vector of interface enumeration
  // types: eType[j] gives the type ID of jth interface. interface
  // between C and D is assigned the same enumeration type as that
  // between A and B if:
  //
  //   1) The interfaces are identified as chemically equivalent,
  //      which means that:
  //     a) A is highly similar or identical to C
  //     b) B is highly similar or identical to D
  //     c) A and B are found in the same relative position as C and D.
  //
  //   2) One of the interfaces should engage automatically upon
  //      engagement of the other, e.g. due to a ligand fix:
  //
  //        ////// A //////    - ligand B is fixed to C, then
  //        ---------------      interface B-A should be automatically
  //        ---------- B --      engaged if interface C-A is engaged.
  //        //// C ////////
  //
  //      Such interface pairs represent a first-order interface
  //      induction.
  //
  //
  //   The assembly-building procedure makes identical treatment of
  // interfaces belonging to the same enumeration type, even though
  // they may not be crystallographically (and in case of induced
  // interfaces, even chemically) identical. This, e.g., applies to
  // the case when the structure has been solved in a lower-symmetry
  // space group resulting in ASU having several identical copies of
  // the molecules.
  //
  int     i,j,k,n,m,et;
  bool done;

    mmdb::GetVectorMemory ( eType,nInterfaces,0 );

    //  1. Find all chemically-equivalent interfaces first

    nETypes = 0;

    for (i=0;i<nInterfaces;i++)
      eType[i] = -1;

    for (i=0;i<nInterfaces;i++)
      if (eType[i]<0)  {
        eType[i] = nETypes;
        for (j=i+1;j<nInterfaces;j++)
          if (eType[j]<0)  {
            if (interface[i]->type==interface[j]->type)
              eType[j] = nETypes;
          }
        nETypes++;
      }

    //  2. Find the first-order induction interfaces:
    //     simply check that domains in jth interface get
    //     connected if ith interface and all fixed-ligand
    //     interfaces are engaged
    for (i=0;i<nInterfaces;i++)
      if (!interface[i]->fixedLigand)
        do {
          done = true;
          for (j=i+1;j<nInterfaces;j++)  {
            et = 1;
            for (m=0;(m<=i) && et;m++)
              if (eType[j]==eType[m])  et = 0;
            if ((!interface[j]->fixedLigand) && et)
              for (k=0;k<nMonomers;k++)  {
                n = 0;
                while (n<_max_n_int)
                  if (M[k]->L[n])  {
                    if (M[k]->L[n][j].M)  {
                      M[k]->L[n][j].M->id2 = 1; // termination id
                      if (traceMonomerChain(M[k],eType[i],0,0,0,
                                            &(M[k]->L[n][j]),0))  {
                        // jth interface engages automatically if
                        // ith interface is engaged. Assign the same
                        // engagement type (which is not, in general,
                        // equivalent to chemical type!) to all
                        // interfaces currently having jth engagement
                        // type
                        et = eType[j];
                        for (m=0;m<nInterfaces;m++)
                          if (eType[m]==et)  eType[m] = eType[i];
                        done = false;
                      }
                      M[k]->L[n][j].M->id2 = 0;   // restore id
                    }
                    n++;
                  } else
                    n = _max_n_int;
              }
          }
        } while (!done);


    //  3. Renumber engagement types so that they go sequentially

    for (i=0;i<nInterfaces;i++)
      eType[i] = -(eType[i]+1);

    nETypes = 0;
    for (i=0;i<nInterfaces;i++)
      if (eType[i]<0)  {
        for (j=i+1;j<nInterfaces;j++)
          if (eType[j]==eType[i])  eType[j] = nETypes;
        eType[i] = nETypes;
        nETypes++;
      }

  }


  #define _mon_eps  0.1

  ASSMB_RC Assembler::MakeCrystalMonomers()  {
  // Generates monomers in case of in-crystal calculations
  PDomain        domain,domain2;
  TFrame         frame1,frame2;
  mmdb::mat44    uct,tm;
  mmdb::realtype tx,ty,tz, x1,y1,z1, x2,y2,z2;
  int            i,j,k,m,n;
  int            ic,jc,kc;
  bool           ok;
  ASSMB_RC       rc;

    //  1. Check that crystallographic information is present

    if (MMDB->CrystReady()>=0)  nSymOps = MMDB->GetNumberOfSymOps();
                          else  nSymOps = 0;

    if (nSymOps<=0)      return  ASSMB_noSymOps;
    if (D->nDomains<=0)  return  ASSMB_noDomains;
    rc = ASSMB_Ok;


    //  2. Calculate translation shifts to all unit cells in
    //     consideration

    n = 2*nCellLayers + 1;
    ucv  = new ppvect3[n];
    ucv += nCellLayers;

    MMDB->GetROMatrix ( uct );

    for (i=-nCellLayers;i<=nCellLayers;i++)  {
      ucv[i]  = new pvect3[n];
      ucv[i] += nCellLayers;
      for (j=-nCellLayers;j<=nCellLayers;j++)  {
        ucv[i][j]  = new mmdb::vect3[n];
        ucv[i][j] += nCellLayers;
        for (k=-nCellLayers;k<=nCellLayers;k++)  {
          ucv[i][j][k][0] = i*uct[0][0] + j*uct[0][1] + k*uct[0][2];
          ucv[i][j][k][1] = i*uct[1][0] + j*uct[1][1] + k*uct[1][2];
          ucv[i][j][k][2] = i*uct[2][0] + j*uct[2][1] + k*uct[2][2];
        }
      }
    }


    //  3. Allocate all monomers in unit cell 333

    nMonAlloc = D->nDomains*nSymOps;

    M = new PMonomer[nMonAlloc];
    for (i=0;i<nMonAlloc;i++)
      M[i] = NULL;

    //  4. Calculate transformation matrices for all monomers

    nMonomers = 0;  //  monomer count

    for (i=0;(i<D->nDomains) && (rc==ASSMB_Ok);i++)  {
      domain = D->domain[i];
      if (domain->assemble)  {
        for (m=0;(m<nSymOps) && (rc==ASSMB_Ok);m++)
          if (MMDB->GetUCTMatrix(uct,m,
              D->domain[domain->ncsParent]->mx,
              D->domain[domain->ncsParent]->my,
              D->domain[domain->ncsParent]->mz,0,0,0)==mmdb::SYMOP_Ok) {

            if (isMat4Rot(uct,1.0e-2))  {

              // check that there is no overlapping with same-type
              // symmetry mates (PDB is full of crap)
              ok = true;

              if (domain->nAtoms>3)  {
                // large molecules, will compare frames
                domain->getReferenceFrame ( frame1,uct );
              } else  {
                // small molecules without symmetry, will compare
                // only mass centers
                x1 = domain->mx;
                y1 = domain->my;
                z1 = domain->mz;
                mmdb::TransformXYZ ( uct,x1,y1,z1 );
              }

              for (j=0;(j<nMonomers) && ok;j++)  {
                domain2 = D->domain[M[j]->ncsParent];
                if (domain2->type==domain->type)  {
                  // compare frames of two monomers of the same type,
                  // checking all symmetry mates in the neighbouring
                  // unit cells
                  tx = M[j]->uct[0][3];
                  ty = M[j]->uct[1][3];
                  tz = M[j]->uct[2][3];
                  for (ic=-nCellLayers;(ic<=nCellLayers) && ok;ic++)
                    for (jc=-nCellLayers;(jc<=nCellLayers) && ok;jc++)
                      for (kc=-nCellLayers;(kc<=nCellLayers) &&
                                                              ok;kc++) {
                        M[j]->uct[0][3] = tx + ucv[ic][jc][kc][0];
                        M[j]->uct[1][3] = ty + ucv[ic][jc][kc][1];
                        M[j]->uct[2][3] = tz + ucv[ic][jc][kc][2];
                        if (domain->nAtoms>3)  {
                          // large molecules, compare frames
                          domain2->getReferenceFrame (frame2,M[j]->uct);
                          ok = !CompareReferenceFrames (
                                                 frame1,frame2,1.0e-7 );
                        } else  {
                          // small molecules without symmetry, compare
                          // only mass centers
                          x2 = domain2->mx;
                          y2 = domain2->my;
                          z2 = domain2->mz;
                          mmdb::TransformXYZ ( M[j]->uct,x2,y2,z2 );
                          ok = ((fabs(x1-x2)>_mon_eps) ||
                                (fabs(y1-y2)>_mon_eps) ||
                                (fabs(z1-z2)>_mon_eps));
                        }
                      }
                  M[j]->uct[0][3] = tx;
                  M[j]->uct[1][3] = ty;
                  M[j]->uct[2][3] = tz;
                }
              }
              if (ok)  {
                // transformation matrix refers to the original unit
                // rather than to an NCS-mate
                mmdb::Mat4Mult ( tm,uct,domain->ncs_m );
                M[nMonomers] = new Monomer ( i,domain->ncsParent,
                                             domain->dclass,
                                             tm,nInterfaces );
                M[nMonomers]->symOpNo = m;
                M[nMonomers]->ix      = nMonomers;
                nMonomers++;
              }
            } else  {
              rc = ASSMB_improperSymOp;
            }
          } else  {
            rc = ASSMB_noSymOp; // no symmetry operation although it
                                // should be there
          }
      }
    }

    return rc;

  }


  ASSMB_RC Assembler::MakeComplexMonomers()  {
  // Generates monomers for given complex analysis
  PDomain        domain,domain2;
  TFrame         frame2;
  mmdb::realtype x2,y2,z2;
  int            i,j;
  bool           ok;

    if (D->nDomains<=0)  return  ASSMB_noDomains;

    nSymOps   = 0;
    nMonAlloc = D->nDomains;

    M = new PMonomer[nMonAlloc];
    for (i=0;i<nMonAlloc;i++)
      M[i] = NULL;

    nMonomers = 0;  //  monomer count

    for (i=0;i<D->nDomains;i++)  {
      domain = D->domain[i];
      if (domain->assemble)  {
        // check that there is no overlapping with same-type
        // symmetry mates (PDB is full of crap)
        ok = true;
        for (j=0;(j<nMonomers) && ok;j++)  {
          domain2 = D->domain[M[j]->ncsParent];
          if (domain2->type==domain->type)  {
            // compare frames of two monomers of the same type,
            // checking all symmetry mates in the neighbouring
            // unit cells
            if (domain->nAtoms>3)  {
              // large molecules, compare frames
              domain2->getReferenceFrame   ( frame2,M[j]->uct );
              ok = !CompareReferenceFrames ( domain->frame,
                                             frame2,1.0e-7 );
            } else  {
              // small molecules without symmetry, compare
              // only mass centers
              x2 = domain2->mx;
              y2 = domain2->my;
              z2 = domain2->mz;
              mmdb::TransformXYZ ( M[j]->uct,x2,y2,z2 );
              ok = ((fabs(domain->mx-x2)>_mon_eps) ||
                    (fabs(domain->my-y2)>_mon_eps) ||
                    (fabs(domain->mz-z2)>_mon_eps));
            }
          }
        }
        if (ok)  {
          M[nMonomers] = new Monomer ( i,domain->ncsParent,
                                       domain->dclass,
                                       domain->ncs_m,
                                       nInterfaces );
          M[nMonomers]->symOpNo = 0;
          M[nMonomers]->ix      = nMonomers;
          nMonomers++;
        }
      }
    }

    return ASSMB_Ok;

  }


  void Assembler::RemoveBricks()  {
  int ix,iy,iz;

    if (MB)  {
      for (ix=0;ix<nx_bricks;ix++)
        if (MB[ix])  {
          for (iy=0;iy<ny_bricks;iy++)
            if (MB[ix][iy])  {
              for (iz=0;iz<nz_bricks;iz++)
                if (MB[ix][iy][iz])
                  delete MB[ix][iy][iz];
              delete[] MB[ix][iy];
            }
          delete[] MB[ix];
        }
      delete[] MB;
      MB = NULL;
    }

    nx_bricks = 0;
    ny_bricks = 0;
    nz_bricks = 0;

    x1_brick  = 0.0;
    x2_brick  = 0.0;
    y1_brick  = 0.0;
    y2_brick  = 0.0;
    z1_brick  = 0.0;
    z2_brick  = 0.0;

  }

  void Assembler::MakeBricks()  {
  mmdb::rvector  x,y,z;
  mmdb::realtype margin;
  int            i,j, ix,iy,iz;

    RemoveBricks();

    x1_brick =  mmdb::MaxReal;
    x2_brick = -mmdb::MaxReal;
    y1_brick =  mmdb::MaxReal;
    y2_brick = -mmdb::MaxReal;
    z1_brick =  mmdb::MaxReal;
    z2_brick = -mmdb::MaxReal;

    mmdb::GetVectorMemory ( x,nMonomers,0 );
    mmdb::GetVectorMemory ( y,nMonomers,0 );
    mmdb::GetVectorMemory ( z,nMonomers,0 );

    for (i=0;i<nMonomers;i++)  {
      M[i]->getMassCenter ( D,x[i],y[i],z[i] );
      if (x[i]<x1_brick)  x1_brick = x[i];
      if (x[i]>x2_brick)  x2_brick = x[i];
      if (y[i]<y1_brick)  y1_brick = y[i];
      if (y[i]>y2_brick)  y2_brick = y[i];
      if (z[i]<z1_brick)  z1_brick = z[i];
      if (z[i]>z2_brick)  z2_brick = z[i];
    }

    brick_size = 2.0*D->Dmax;
    margin = 2.1*brick_size;
    x1_brick -= margin;
    x2_brick += margin;
    y1_brick -= margin;
    y2_brick += margin;
    z1_brick -= margin;
    z2_brick += margin;

    nx_bricks = mmdb::ifloor ( (x2_brick-x1_brick)/brick_size ) + 1;
    ny_bricks = mmdb::ifloor ( (y2_brick-y1_brick)/brick_size ) + 1;
    nz_bricks = mmdb::ifloor ( (z2_brick-z1_brick)/brick_size ) + 1;

    for (i=0;i<nMonomers;i++)  {
      ix = mmdb::ifloor ( (x[i]-x1_brick)/brick_size );
      iy = mmdb::ifloor ( (y[i]-y1_brick)/brick_size );
      iz = mmdb::ifloor ( (z[i]-z1_brick)/brick_size );
      if (!MB)  {
        MB = new PPPBrick[nx_bricks];
        for (j=0;j<nx_bricks;j++)
          MB[j] = NULL;
      }
      if (!MB[ix])  {
        MB[ix] = new PPBrick[ny_bricks];
        for (j=0;j<ny_bricks;j++)
          MB[ix][j] = NULL;
      }
      if (!MB[ix][iy])  {
        MB[ix][iy] = new PBrick[nz_bricks];
        for (j=0;j<nz_bricks;j++)
          MB[ix][iy][j] = NULL;
      }
      if (!MB[ix][iy][iz])
        MB[ix][iy][iz] = new Brick();
      MB[ix][iy][iz]->AddObject ( i );
    }

    mmdb::FreeVectorMemory ( z,0 );
    mmdb::FreeVectorMemory ( y,0 );
    mmdb::FreeVectorMemory ( x,0 );

  }

  void Assembler::getBricks ( PPBrick brick, mmdb::mat44 & uct,
                               PDomain dom )  {
  mmdb::realtype x,y,z;
  int      i,i1,i2, j,j1,j2, k,k1,k2, n;

    for (i=0;i<27;i++)
      brick[i] = NULL;

    if (!MB)  return;

    x = uct[0][0]*dom->mx + uct[0][1]*dom->my +
        uct[0][2]*dom->mz + uct[0][3];
    i = mmdb::ifloor ( (x-x1_brick)/brick_size );
    if ((i<0) || (i>=nx_bricks))  return;

    y = uct[1][0]*dom->mx + uct[1][1]*dom->my +
        uct[1][2]*dom->mz + uct[1][3];
    j = mmdb::ifloor ( (y-y1_brick)/brick_size );
    if ((j<0) || (j>=ny_bricks))  return;

    z = uct[2][0]*dom->mx + uct[2][1]*dom->my +
        uct[2][2]*dom->mz + uct[2][3];
    k = mmdb::ifloor ( (z-z1_brick)/brick_size );
    if ((k<0) || (k>=nz_bricks))  return;

    n  = 0;
    i1 = mmdb::IMax ( i-1,0 );
    i2 = mmdb::IMin ( i+2,nx_bricks );
    j1 = mmdb::IMax ( j-1,0 );
    j2 = mmdb::IMin ( j+2,ny_bricks );
    k1 = mmdb::IMax ( k-1,0 );
    k2 = mmdb::IMin ( k+2,nz_bricks );
    for (i=i1;i<i2;i++)
      if (MB[i])
        for (j=j1;j<j2;j++)
          if (MB[i][j])
            for (k=k1;k<k2;k++)
              if (MB[i][j][k])  {
                if (MB[i][j][k]->nObjects>0)
                  brick[n++] = MB[i][j][k];
              }

  }


  void Assembler::CalcCrystalInterfaces()  {
  //   This function calculates interfaces between all monomers,
  // treating the outer unit cells in a wrap-around mode. This
  // function should be called only after MakeMonomers() completes.
  PBrick brick[27];
  PDomain    D1,D2;
  mmdb::mat44       ucti;
  DFrame      dframe;
  int         i,j,k, n,n1,n2, ic,jc,kc, ji,iNo;

    MakeBricks();

    for (i=0;i<nMonomers;i++)  {
      n1 = M[i]->ncsParent;
      for (j=0;j<4;j++)
        for (k=0;k<4;k++)
          ucti[j][k] = M[i]->uct[j][k];
      for (ic=-nCellLayers;ic<=nCellLayers;ic++)
        for (jc=-nCellLayers;jc<=nCellLayers;jc++)
          for (kc=-nCellLayers;kc<=nCellLayers;kc++)  {
            ucti[0][3] = M[i]->uct[0][3] + ucv[ic][jc][kc][0];
            ucti[1][3] = M[i]->uct[1][3] + ucv[ic][jc][kc][1];
            ucti[2][3] = M[i]->uct[2][3] + ucv[ic][jc][kc][2];
            getBricks ( brick,ucti,D->domain[n1] );
            for (k=0;(k<27) && brick[k];k++)
              for (j=0;j<brick[k]->nObjects;j++)  {
                n = brick[k]->n[j];
                if ((ic || jc || kc || (n!=i)) && (n>=i))  {
                  n2 = M[n]->ncsParent;
                  if (D->domain[n1]->isContact(ucti,
                                 D->domain[n2],M[n]->uct))  {
                    MakeDFrame ( dframe,D->domain[n1],ucti,
                                        D->domain[n2],M[n]->uct );
                    iNo = -1;
                    for (ji=0;(ji<nInterfaces) && (iNo<0);ji++)  {
                      D1 = D->domain[interface[ji]->domain1];
                      D2 = D->domain[interface[ji]->domain2];
                      if ((D1->ncsParent==n1) && (D2->ncsParent==n2)) {
                        if (interface[ji]->checkDFrame(dframe,true))
                          iNo = ji;
                      }
                      if ((iNo<0) && (D1->ncsParent==n2) &&
                                     (D2->ncsParent==n1))  {
                        if (interface[ji]->checkDFrame(dframe,false))
                          iNo = ji;
                      }
                    }
                    if (iNo>=0)  {
                      M[i]->AddInterface ( M[n],iNo,-ic,-jc,-kc,
                                                nInterfaces );
                      if (i!=n)
                        M[n]->AddInterface ( M[i],iNo,ic,jc,kc,
                                                  nInterfaces );
                    }
                  }
                }
              }
          }
    }

    RemoveBricks();

  }


  void Assembler::CalcComplexInterfaces()  {
  //   This function calculates interfaces between all monomers
  // of given complex, neglecting a possible crystal environment.
  // This function should be called only after MakeComplexMonomers()
  // completes.
  PBrick  brick[27];
  PDomain     D1,D2;
  DFrame       dframe;
  int          i,j,k, n,n1,n2, ji,iNo;

    MakeBricks();

    for (i=0;i<nMonomers;i++)  {
      n1 = M[i]->ncsParent;
      getBricks ( brick,M[i]->uct,D->domain[n1] ); // use parents only
      for (k=0;(k<27) && brick[k];k++)
        for (j=0;j<brick[k]->nObjects;j++)  {
          n = brick[k]->n[j];
          if (n>i)  {
            n2 = M[n]->ncsParent;  // again, this is parent (no NCS)
            if (D->domain[n1]->isContact(M[i]->uct,
                                 D->domain[n2],M[n]->uct))  {
              MakeDFrame ( dframe,D->domain[n1],M[i]->uct,
                                  D->domain[n2],M[n]->uct );
              iNo = -1;
              for (ji=0;(ji<nInterfaces) && (iNo<0);ji++)  {
                D1 = D->domain[interface[ji]->domain1];
                D2 = D->domain[interface[ji]->domain2];
                if ((D1->ncsParent==n1) && (D2->ncsParent==n2)) {
                  if (interface[ji]->checkDFrame(dframe,true))
                    iNo = ji;
                }
                if ((iNo<0) && (D1->ncsParent==n2) &&
                               (D2->ncsParent==n1))  {
                  if (interface[ji]->checkDFrame(dframe,false))
                    iNo = ji;
                }
              }
              if (iNo>=0)  {
                M[i]->AddInterface ( M[n],iNo,0,0,0,nInterfaces );
                if (i!=n)
                  M[n]->AddInterface ( M[i],iNo,0,0,0,nInterfaces );
              }
            }
          }
        }
    }

    RemoveBricks();

  }

  int Assembler::FixInterface ( int   iNo, mmdb::ovector fixed,
                                PMonRef G, mmdb::ovector mmol ) {
  //   This function attempts to fix interface iNo and returns true
  // in case of success. Returns:
  //    0 - interface iNo, same-type interfaces and all induced
  //        interfaces are fixed
  //    1 - interface iNo cannot be fixed due to crystallographic
  //        considerations (periodic loop over unit cell)
  //    2 - interface iNo cannot be fixed because macromolecules
  //        get linked through a ligand
  //    3 - interface iNo cannot be fixed because macromolecules
  //        get linked immediately
  int  i,j,k,m,n,jj,kk, nMM, i1,j1,k1, t, lgroup, rc;
  bool done;

    if (nMonomers<=1)  return 1;

    // 1. Remember current fixing for possible rollback and
    //    fix all interfaces of given type

    for (i=0;i<nInterfaces;i++)  {
      fixed[i] = interface[i]->fixedLigand;
      if (interface[i]->type==interface[iNo]->type)
        interface[i]->fixedLigand = true;
    }

    // 2. Identify and make induced fixings

    rc = 0;

    do { // iterate until no more new fixings is done
         // or fixing goes wrong

      done = true;

      for (i=0;i<nMonomers;i++)
        M[i]->id2 = 0;

      lgroup = 0;
      for (i=0;i<nMonomers;i++)
        if ((!M[i]->id2) && (M[i]->symOpNo==0) && M[i]->L[0] && (!rc)) {

          //  2.1 Identify a group of monomers connected by
          //      fixed interfaces

          lgroup++;
          m      = 0;  // group size
          nMM    = 0;  // number of macromolecules in the group

          M[i]->id2 = lgroup;  // group id
          if (mmol[M[i]->ix])  nMM++;
          G[m].M = M[i];
          G[m].i = 0;  // so that we start with a monomer in principal
          G[m].j = 0;  // asu (symOpNo==0) of the principal unit cell
          G[m].k = 0;  // (all cell shifts are zero)
          m++;

          // iteratively include all monomers that are connected
          // by fixed interfaces
          j = 0;
          while ((j<m) && (!rc))  {
            for (k=0;(k<nInterfaces) && (!rc);k++)
              if (interface[k]->fixedLigand)  {
                n  = 0;
                while ((n<_max_n_int) && (!rc))
                  if (G[j].M->L[n])  {
                    if (G[j].M->L[n][k].M)  {
                      // make a candidate group member
                      G[m].M = G[j].M->L[n][k].M;
                      G[m].i = G[j].i + G[j].M->L[n][k].i;
                      G[m].j = G[j].j + G[j].M->L[n][k].j;
                      G[m].k = G[j].k + G[j].M->L[n][k].k;
                      // check whether this monomer is already
                      // in the group
                      kk = 0;
                      for (jj=0;(jj<m) && (!kk);jj++)
                        if ((G[jj].M->ix==G[m].M->ix) &&
                            (G[jj].i==G[m].i) &&
                            (G[jj].j==G[m].j) &&
                            (G[jj].k==G[m].k))  kk = 1;
                      if (!kk)  {
                        // new monomer to the group, check whether
                        // this is a macromolecule
                        if (mmol[G[m].M->ix])  {
                          nMM++;
                          if (nMM>1)  rc = 2; // not allowed!
                        } else  {
                          // commit
                          if (G[m].M->symOpNo==0) // this is to skip it
                            G[m].M->id2 = lgroup; // on next i
                          m++;
                          if (m>=nMonomers)  rc = 1; // over UC size
                        }
                      }
                    }
                    n++;
                  } else
                    n = _max_n_int;
              }
            j++;
          }

          //  6.2 Fix all interfaces within identified group

          for (j=0;(j<m) && (!rc);j++)
            for (k=0;k<nInterfaces;k++)
              if (!interface[k]->fixedLigand)  {
                // find out whether this interface links two monomers
                // from the group
                n  = 0;
                while ((n<_max_n_int) && (!rc))
                  if (G[j].M->L[n])  {
                    if (G[j].M->L[n][k].M)  {
                      // the interface is linking if G[j].M->L[n][k].M
                      // is in the group
                      i1 = G[j].i + G[j].M->L[n][k].i;
                      j1 = G[j].j + G[j].M->L[n][k].j;
                      k1 = G[j].k + G[j].M->L[n][k].k;
                      kk = 0;
                      for (jj=0;(jj<m) && (!kk);jj++)
                        if ((G[jj].M->ix==G[j].M->L[n][k].M->ix) &&
                            (G[jj].i==i1) &&
                            (G[jj].j==j1) &&
                            (G[jj].k==k1))  kk = 1;
                      if (kk)  {
                        if (mmol[G[j].M->L[n][k].M->ix] &&
                            mmol[G[j].M->ix])  {
                          // induced fixing binds macromolecules,
                          // for which we do not allow
                          rc = 3;
                        } else  {
                          // fix all interfaces of k's type
                          for (t=0;t<nInterfaces;t++)
                            if (interface[t]->type==interface[k]->type)
                              interface[t]->fixedLigand = true;
                          done = false;  // keep iterating
                          n = _max_n_int; // quit for this interface
                        }
                      }
                    }
                    n++;
                  } else
                    n = _max_n_int;
              }

        }

    } while ((!rc) && (!done));

    if (rc)  // rollback
      for (i=0;i<nInterfaces;i++)
        interface[i]->fixedLigand = fixed[i];

    return rc;

  }

  void Assembler::FixLigands()  {
  //   This function identifies ligands that may be tied up to
  // macromolecules and other ligands. These include unconditionally
  // the molecules making covalent linking, and others excluding
  // those making crystallographically-related interfaces. Eligible
  // ligands are then fixed to macromolecules and to other ligand by
  // raising the "fixed" flag in the corresponding interfaces.
  PMonRef       G;
  mmdb::ivector nd;
  mmdb::ovector fixed,mmol;
  int           i,j,k;

    // Identify crystallographically-related interfaces (just in case)

    for (i=0;i<nInterfaces;i++)
      interface[i]->Xrel = false;

    for (i=0;i<nMonomers;i++)
      if (M[i]->L[0])  {
        // the above condition means that ith monomer makes
        // at least one interface. Loop over all interfaces
        // made by ith monomer
        for (j=0;j<nInterfaces;j++)
          if ((!interface[j]->Xrel) && M[i]->L[0][j].M)  {
            if (M[i]->L[1])  {
              if (M[i]->L[1][j].M)  {
                // crystallographically-related interfaces
                for (k=0;k<nInterfaces;k++)
                  if (interface[k]->type==interface[j]->type)
                    interface[k]->Xrel = true;
              }
            }
          }
      }

    if (nMonomers<=1)  return;

    //  1. Allocated working memory

    mmdb::GetVectorMemory ( fixed,nInterfaces,0 );
    mmdb::GetVectorMemory ( mmol ,nMonomers  ,0 );
    G = new MonRef[nMonomers];

    //  2. Initially, set all interfaces unfixed

    for (i=0;i<nInterfaces;i++)
      interface[i]->fixedLigand = false;

    //  3. Firstly, fix covalent linkage between ligands and
    //     macromolecules

    //  Set macromolecuar flags. They could be set to none on
    //  this stage if covalent linking between macromolecules
    //  may be trusted
    for (i=0;i<nMonomers;i++)
      mmol[i] = (M[i]->dclass!=DCLASS_Ligand);

    for (i=0;i<nInterfaces;i++)
      if ((!interface[i]->fixedLigand) && (interface[i]->nCovBonds>0) &&
          (((interface[i]->dclass1==DCLASS_Ligand) &&
            (interface[i]->dclass2!=DCLASS_Ligand)) ||
           ((interface[i]->dclass1!=DCLASS_Ligand) &&
            (interface[i]->dclass2==DCLASS_Ligand))))
        FixInterface ( i,fixed,G,mmol );

    //  4. Secondly, fix covalent linkage between ligands

    for (i=0;i<nInterfaces;i++)
      if ((!interface[i]->fixedLigand) && (interface[i]->nCovBonds>0) &&
          (interface[i]->dclass1==DCLASS_Ligand) &&
          (interface[i]->dclass2==DCLASS_Ligand))
        FixInterface ( i,fixed,G,mmol );

    //  5. Finally, fix all other fixable interfaces

    if (ligKey!=LIGANDS_FreeAll)  {

      if (ligKey==LIGANDS_Auto)  {
        //   Mark monomers that should not be fixed: simply declare
        // them macromolecules

        mmdb::GetVectorMemory ( nd,D->nDomains,0 );

        for (i=0;i<D->nDomains;i++)  {
          nd[i] = 0;
          for (j=0;j<D->nDomains;j++)
            if (D->isEquivalent(i,j))  nd[i]++;
        }

        for (i=0;i<nMonomers;i++)
          if (M[i]->dclass==DCLASS_Ligand)  {
            j = M[i]->ncsParent;
            mmol[i] = ((D->domain[j]->nAtoms>=fixLigandSize) &&
                       (nd[j]<=fixLigandNum));
          }

        mmdb::FreeVectorMemory ( nd,0 );

      }

      for (i=0;i<nInterfaces;i++)
        if ((!interface[i]->fixedLigand) &&
            ((interface[i]->dclass1==DCLASS_Ligand) ||
             (interface[i]->dclass2==DCLASS_Ligand)))
          FixInterface ( i,fixed,G,mmol );

    }

    //  6. Release working memory

    delete[] G;
    mmdb::FreeVectorMemory ( mmol ,0 );
    mmdb::FreeVectorMemory ( fixed,0 );

    //  7. Reset monomer indicators

    for (i=0;i<nMonomers;i++)
      M[i]->id2 = 0;

  }

  void Assembler::CalcParallelMonomers()  {
  //   This function builds matrix parMon such that if
  // parMon[i][j]==true then molecules i and j are equivalent
  // and parallel.
  TFrame   f1,f2;
  mmdb::realtype eps;
  int      i,j;

    mmdb::GetMatrixMemory ( parMon,nMonAlloc,nMonAlloc,0,0 );

    //   Set low tolerance if NCS mates are present. This is
    // because: a) NCS-mates are supposed to be not
    // crystallographically equivalent, therefore even close
    // orientations should be taken as different b) NCS mates
    // are normally used in viral structures and often they
    // are many as in, e.g. 1w8x (900 NCS mates). Simply
    // because of large numbers, some orientations may be
    // accidentally close.
    //   Without NCS mates, we want to allow for greater
    // deviations, which are required due to imperfectness
    // of both data and crystal. eps=0.075 has been chosen in
    // extreme case of 1kqc, where a tetramer may be wrongly
    // formed by two almost parallel dimers.
    if (D->nDomains>D->nNCSParents)  eps = 0.001;
                               else  eps = 0.075;

    for (i=0;i<nMonomers;i++)  {
      parMon[i][i] = true;
      for (j=i+1;j<nMonomers;j++)  {
        if (M[i]->id==M[j]->id)  {
          // identical monomers are to be checked simply through
          // comparison of their orientations
          parMon[i][j] = M[i]->isParallel ( M[j] );
        } else if (D->isEquivalent(M[i]->id,M[j]->id))  {
          // molecules are not crystallographically identical,
          // but they are physically equivalent
          D->domain[M[i]->ncsParent]->getReferenceFrame ( f1,M[i]->uct );
          D->domain[M[j]->ncsParent]->getReferenceFrame ( f2,M[j]->uct );
          parMon[i][j] = ParallelFrames ( f1,f2,eps );
        } else  // molecules are not parallel or identical/equivalent
          parMon[i][j] = false;
        parMon[j][i] = parMon[i][j];
      }
    }

  }


  ASSMB_RC Assembler::MakePAGraph ( mmdb::PManager MMDBManager,
                                    PDomains       Domains,
                                    PInterfaces    PI )  {
  //  This function builds a Protein Assembly Graph, that is
  // a linked list of molecules and interfaces. Both hierarchical
  // (unit cell -> asymmetric unit -> molecule) and serial
  // indexes are created.
  int  i,j,k;
  bool ok;

    Reset();

    //  1. Retain pointers to domains and structure

    D    = Domains;
    MMDB = MMDBManager;

    //   Get orthogonalisation matrix. It is needed for calculating
    // the transformation matrices, which put monomers onto their
    // 3D positions in multimers.
    if (Crystal)  MMDB->GetROMatrix ( rom );
            else  mmdb::Mat4Init ( rom );

    //  2. Generate monomers in the unit cell

    if (Crystal)  retcode = MakeCrystalMonomers();
            else  retcode = MakeComplexMonomers();
    if (retcode)  return retcode;


    //  3. Check interfaces

    if (!PI) return ASSMB_noInterfaces;

    if ((!MMDBManager) || (!Domains))
      return ASSMB_incompleteData;

    k = PI->AnalyseInterfaces ( Domains );
    if ((k==INTS_notCalculated) || (k==INTS_noInterfaces))
                             return ASSMB_noInterfaces;
    if (k & INTS_mmOverlap)  return ASSMB_Overlap;

    //  4. Make a local copy of interface properties vector
    //
    //   The vector is ordered by decreasing the stabilization energy.
    // All assembly calculations are done with thus resorted interface
    // vector, but the original interface numbering is restored when
    // results are being copied into resulting structures
    // (see CMolecule::_copy(..)).

    retcode = MakeInterfaces ( PI,D->domain );
    if (retcode)  return retcode;


    //  5. Calculate matrix of interfaces between all monomers

    repeated_assignment = false;
    if (Crystal)  CalcCrystalInterfaces();
            else  CalcComplexInterfaces();
    if (repeated_assignment)  return ASSMB_repeatedAssignment;
  //  if (Crystal)  CalInterfaceStats ( PI );


    //  6. Fix ligands if they are too many

    FixLigands();

    //  7. Calculate interface enumeration types

    CalcEnumerationTypes();

    //  8. Identify parallel monomers.

    CalcParallelMonomers();

    //  9. Non-compulsory but cheap checks (may be removed)

    ok = true;
    for (i=0;i<nMonomers;i++)
      if (M[i]->L[0])
        for (j=i+1;j<nMonomers;j++)
          if (M[j]->L[0] && (M[i]->id==M[j]->id))  {
            for (k=0;k<nInterfaces;k++)  {
              if ((( M[i]->L[0][k].M) && (!M[j]->L[0][k].M)) ||
                  ((!M[i]->L[0][k].M) && ( M[j]->L[0][k].M)))
                ok = false;
              if (M[i]->L[1] && M[j]->L[1])  {
                if ((( M[i]->L[1][k].M) && (!M[j]->L[1][k].M)) ||
                    ((!M[i]->L[1][k].M) && ( M[j]->L[1][k].M)))
                  ok = false;
              }
            }
        }
    if (!ok)  return  ASSMB_brokenComposition;

    for (i=0;i<nMonomers;i++)
      if (M[i]->L[0])
        for (j=0;j<nInterfaces;j++)  {
          if (M[i]->L[0][j].c && M[i]->L[0][j].M) {
            if (M[i]->L[0][j].M->L[0] && M[i]->L[0][j].M->L[1])  {
              if ((M[i]->L[0][j].M->L[0][j].M!=M[i]) &&
                  (M[i]->L[0][j].M->L[1][j].M!=M[i]))
                ok = false;
            }
            if (M[i]->L[1])  {
              if ((M[i]->L[0][j].M==M[i]->L[1][j].M) &&
                  (M[i]->L[0][j].i==M[i]->L[1][j].i) &&
                  (M[i]->L[0][j].j==M[i]->L[1][j].j) &&
                  (M[i]->L[0][j].k==M[i]->L[1][j].k))
                ok = false;
            }
          }
          if (M[i]->L[1])  {
            if (M[i]->L[1][j].c && M[i]->L[1][j].M)  {
              if (M[i]->L[1][j].M->L[0] && M[i]->L[1][j].M->L[1]) {
                if ((M[i]->L[1][j].M->L[0][j].M!=M[i]) &&
                    (M[i]->L[1][j].M->L[1][j].M!=M[i]))
                  ok = false;
              }
            }
          }
        }

    if (!ok)  retcode = ASSMB_brokenComplementarity;

    return retcode;

  }


  void Assembler::SeedMultimers()  {
  //   The function generates an initial set of multimeric clusters for
  //  further combinatorial search.
  int  i;

    DeleteMultimers ( U,nMultimers,nMultAlloc );
    DeleteComplex   ();

    //  1. Allocate the theoretically maximal number of all
    //     multimers in the simulation, that is the total
    //     number of monomers.

    nMultAlloc = nMonomers;
    nMultimers = nMonomers;
    U = new PMultimer[nMultAlloc];

    //  2. Generate initial multimers

    for (i=0;i<nMonomers;i++)  {
      U[i] = new Multimer ( nMonomers );
      U[i]->MakeMultimer ( M[i] );
    }

  }


  int  Assembler::makeInterfaceFlags ( mmdb::ivector iflag )  {
  //  Set up the initial interface engagement flags
  mmdb::ovector s;
  int     i,nopen;

    //  1. Mark all interface types as available for engaging

    mmdb::GetVectorMemory ( s,nETypes,0 );

    for (i=0;i<nETypes;i++)  {
      iflag[i] = INTF_Undefined;
      s[i]     = true;
    }

    //  2. Mark fixed ligand and non-affine/casual interface types.
    //     Fixed interfaces will be engaged beforeahead of assembly
    //     search and will never be attempted to dissociate. Casual
    //     interfaces can not start a branch of the search tree.
    //     The algorithm will, however, engage casual interfaces
    //     as usual in all other branches.

    for (i=0;i<nInterfaces;i++)  {
      if (interface[i]->fixedLigand)
        iflag[eType[i]] = INTF_Fixed;
      else if ((iflag[eType[i]]==INTF_Undefined) &&
               ((interface[i]->stabEn>=0.0) || interface[i]->casual))
        iflag[eType[i]] = INTF_Off;
      else
        iflag[eType[i]] = INTF_Open;
      if ((D->domain[interface[i]->domain1]->dclass!=DCLASS_Ligand) ||
          (D->domain[interface[i]->domain2]->dclass!=DCLASS_Ligand))
        s[eType[i]] = false;
    }

    nopen = 0;
    for (i=0;i<nETypes;i++)  {
      if ((iflag[i]==INTF_Open) && s[i])
        iflag[i] = INTF_Off;
      if (iflag[i]==INTF_Open)  nopen++;
    }

    mmdb::FreeVectorMemory ( s,0 );

    return nopen;

  }

  int  Assembler::MakeRecStack()  {
  int i,nopen;

    DeleteRecStack();

    nRSAlloc = nInterfaces + 1;
    recStack = new PIntfRecStack[nRSAlloc];
    for (i=0;i<nRSAlloc;i++)
      recStack[i] = new IntfRecStack ( nETypes,nMultimers );

    nopen = makeInterfaceFlags ( recStack[0]->iflag );

    for (i=0;i<nMultimers;i++)  {
      recStack[0]->mSize[i] = U[i]->mSize;
      recStack[0]->U[i]     = U[i];
    }
    recStack[0]->nMultimers = nMultimers;

    return nopen;

  }


  ASSMB_RC Assembler::CalcAssemblies()  {
  int i,r;

    //  1. Initialize assembly seeds - assemblies made of just
    //     a single molecule

    SeedMultimers();

    //  2. Initialize recursion stack

    r = MakeRecStack();
    if (r>maxNofETypes)  return ASSMB_tooBigSystem;

    //  3. Engage fixed interfaces

    r = 0;  // initial recursion level
    for (i=0;i<nETypes;i++)
      if (recStack[r]->iflag[i]==INTF_Fixed)  {
        if (EngageInterface(i,r))
              r++;
        else  recStack[r]->iflag[i] = INTF_Open;
      }


    //  4. Run the recursive interface enumeration

    nMultSets         = 0;     // number of resulting multimer sets
    multimer_overflow = false;
    maxStableSize     = 0;

    runStatus         = RSTAT_Ok;
    NoteMultimer ( r );

    startCPUclock = clock();

    for (i=0;(i<nETypes) && (runStatus==RSTAT_Ok);i++)
      if (recStack[r]->iflag[i]==INTF_Open)  {
        noImproveCnt  = 0;
        EnumerateInterfaces ( r,i );
        if (runStatus==RSTAT_Ok)        runStatus = RSTAT_Finish;
        if (runStatus==RSTAT_QuitLoop)  runStatus = RSTAT_Ok;
      }
    CalcEqIds();


    //  4. Roll back fixed interfaces

    while (r>0)  {
      RollBack ( r );
      r--;
    }

    //  5. Delete recursion stack

    DeleteRecStack();

    if (runStatus==RSTAT_TimeLimit)
                            return ASSMB_timeLimit;
    if (multimer_overflow)  return ASSMB_assemblyOverflow;

    return ASSMB_Ok;

  }

  ASSMB_RC Assembler::CalcAssemblies ( mmdb::PManager MMDBManager,
                                       PDomains       Domains,
                                       PInterfaces    PI )  {
    Crystal = true;  // analyse crystal
    retcode = MakePAGraph ( MMDBManager,Domains,PI );
    if (retcode==ASSMB_Ok)
      retcode = CalcAssemblies();
    return retcode;
  }


  ASSMB_RC Assembler::AnalyseComplex ( mmdb::PManager MMDBManager,
                                       PDomains       Domains,
                                       PInterfaces    PI )  {
  //    Analyse complex given in coordinate part of input file.
  //  The coordinate part includes also NCS-mates.
  //    This function may be called ONLY after CalcAssemblies().
  PPMonomer     Mon0;
  PPMultimer    U0;
  PPMultimerSet multSet0; // temporary buffer
  int           nMons,nMonA,nMults,nMultA,nMSets;//,nMSetsA;
  int           i,j,id;

    DeleteComplex();

    // retain current structures for assigning complex id
    // after it was computed

    nMons     = nMonomers;
    nMonA     = nMonAlloc;
    Mon0      = M;
    M         = NULL;
    nMonomers = 0;
    nMonAlloc = 0;

    nMults    = nMultimers;
    nMultA    = nMultAlloc;
    U0        = U;
    U         = NULL;
    nMultimers = 0;
    nMultAlloc = 0;

    nMSets    = nMultSets;
  //  nMSetsA   = nMultSetAlloc;
    multSet0  = multSet;
    multSet   = NULL;
    nMultSets = 0;
    nMultSetAlloc = 0;

    Crystal = false;  // analyse given assembly
    
    retcode = MakePAGraph ( MMDBManager,Domains,PI );

    if (retcode<ASSMB_Overlap)  {

      Complex = new Multimer ( 0 );
      Complex->MakeBoundMultimer ( M,nMonomers );
      if (Complex->mSize<=0)  {
        //   The complex is either a single monomer or all
        // monomers are separated. Make the following call only
        // for safety reasons (otherwise pisa will crash), however
        // all the results will be meaningless.
        Complex->MakeMultimer ( M,nMonomers );
      }

      Complex->CalcInternalInterfaces ( nInterfaces );
      Complex->CalcProperties ( D,interface,eType,nInterfaces,
                                nETypes,&SNum,rom );

      
      if (multimer_overflow)  retcode = ASSMB_assemblyOverflow;

      Complex->type = -1;

      id = -1;
      for (i=0;(i<nMSets) && (Complex->type<0);i++)
        for (j=0;(j<multSet0[i]->nMultimers) &&
                 (Complex->type<0);j++)  {
          id = mmdb::IMax ( id,multSet0[i]->U[j]->type );
          if (multSet0[i]->U[j]->isIdentical(Complex,D,rom,false))
            Complex->type = multSet0[i]->U[j]->type;
        }

      if (Complex->type<0)
        Complex->type = id;

    }

    DeleteMonomers     ( Mon0   ,nMons    ,nMonA         );
    DeleteMultimers    ( U0     ,nMults   ,nMultA        );
    DeleteMultimerSets ( multSet,nMultSets,nMultSetAlloc );

    return retcode;

  }


  bool Assembler::EngageInterface ( int p, int r )  {
  //   This function attemts to engage interfaces of type p in the
  // crystal. An interface can not be engaged if doing so results
  // in an assembly containing two alike molecules in parallel
  // orientations. Engaging an interface may make other interfaces
  // internal to some assemblies, therefore the function seeks such
  // induced interfaces and tries to engage them all until no new
  // inductions can be made. The function sets iflag=p for all
  // interface types that have been attempted to engage or identified
  // as induced.
  //   The function returns true if all interfaces of type p,
  // as well as all induced interfaces, have been engaged (and
  // assemblies merged). If interface p or any induced
  // interface cannot be engaged, no assemblies get merged
  // (function RollBack(..) used).
  int      i,j,k,p1,r1;
  bool  success;


    //  1. Initiate interface engagement flags on the next recursion
    //     level and stabilization energy collectors

    r1 = r + 1;
    for (i=0;i<nETypes;i++)
      recStack[r1]->iflag[i] = recStack[r]->iflag[i];

    //  2. Engage interface number p and any interfaces equivalent to p

    success = true;
    for (i=0;i<nInterfaces;i++)
      if (eType[i]==p)  {
        for (j=0;(j<recStack[r]->nMultimers) && success;j++)
          if (recStack[r]->U[j]->mSize>0)  {
            k = recStack[r]->U[j]->EngageInterface (
                            i,parMon,interface[i]->intArea );
            if (k<0) success = false;
          }
      }

    recStack[r1]->iflag[p] = p;

    //  2. Regardless of the success, check for all possible
    //     induced interfaces, including multiple induction

    p1 = 0;  // to be the induced interface number
    while (p1>=0)  {

      //  2.1 Look for an induced internal interface

      p1 = -1;
      for (i=0;(i<nETypes) && (p1<0);i++)
        if (recStack[r1]->iflag[i]<=INTF_Open)  { // INTF_Off also
                                                  //     checked here!
          for (j=0;(j<recStack[r]->nMultimers) && (p1<0);j++)
            if (recStack[r]->U[j]->mSize>0)  {
              if (recStack[r]->U[j]->isInternal(i,eType,nInterfaces)) {
                // induced interface type 'i' has been found,
                // so quit the loop
                if (recStack[r1]->iflag[i]==INTF_Void)  {
                  // but this interface cannot be engaged, therefore
                  // set success off and keep looking for all other
                  // induced interfaces
                  success = false;
                } else  {
                  // mark this interface as induced by interface 'p'
                  p1 = i;
                }
              }
            }
        }

      if (p1>=0)  {  // induced interface type found

        //  2.2  Engage the induced interface and all interfaces
        //       equivalent to it

        for (i=0;i<nInterfaces;i++)
          if (eType[i]==p1)  {
            for (j=0;(j<recStack[r]->nMultimers) && success;j++)
              if (recStack[r]->U[j]->mSize>0)  {
                k = recStack[r]->U[j]->EngageInterface (
                           i,parMon,interface[i]->intArea );
                if (k<0)  success = false;
              }
          }
        recStack[r1]->iflag[p1] = p; // p1!=p <=> induced interface

      }

    }

    if (success)  {

      //  3. interface 'p' and all induced interfaces have been
      //     successfully engaged, prepare stack for the next
      //     recursion level

      r1 = r+1;
      j  = 0; // j counts assemblies in stack position r1
      for (i=0;i<recStack[r]->nMultimers;i++)
        if (recStack[r]->U[i]->mSize>0)  {
          recStack[r1]->U    [j] = recStack[r]->U[i];
          recStack[r1]->mSize[j] = recStack[r1]->U[j]->mSize;
          j++;
        }
      recStack[r1]->nMultimers = j;

    } else  {

      //  4. Engaging the interface 'p' or one of the induced
      //     interfaces was not sucessful. Roll back
      //     (unmerge assemblies).
      RollBack ( r );

      //  Mark all interfaces of type p as non-engageable in this
      // and all further configurations
      recStack[r]->iflag[p] = INTF_Void;

    }

    return success;

  }


  void Assembler::RollBack ( int r )  {
  int i;

    for (i=0;i<recStack[r]->nMultimers;i++)  {
      if (recStack[r]->mSize[i]<=recStack[r]->U[i]->mSize)
        recStack[r]->U[i]->mSize = recStack[r]->mSize[i];
      else  {
        recStack[r]->U[i]->mSize = recStack[r]->mSize[i];
        recStack[r]->U[i]->SetMultimerReferences();
      }
      recStack[r]->U[i]->CalcInternalInterfaces ( nInterfaces );
    }

  }


  void Assembler::EnumerateInterfaces ( int r, int i0 )  {
  //   Recursive procedure making all possible multimers in the crystal,
  // r is the recursion level. The assemblies are build by sequential
  // engaging (glueing) interfaces between them in all possible
  // combinations. Due to the crystal symmetry, engaging an interface
  // between two multimers means that all other interfaces of the same
  // type (i.e. those between particular molecules in particular relative
  // position) in crystal must be engaged. An interface cannot be
  // engaged if doing so yields a multimer containing alike monomers
  // in parallel orientations. In a correctly described crystal, this
  // also means that multimer size cannot exceed the size of a unit cell.
  int i,j,k,r1,setNo;

    setNo = nMultSets-1;

    if (setNo>=0)  {

      r1 = r+1;
      for (i=i0;(i<nETypes) && (runStatus==RSTAT_Ok);i++)  {
        if (multSet[setNo]->getMaxEType(i,dissThreshold,eType,
                                        interface,nInterfaces)<i)
          break;
        if (recStack[r]->iflag[i]==INTF_Open)  {
          if (EngageInterface(i,r))  {

            // The interface has been engaged together with all
            // interfaces induced by this engagement. iflags[] of all
            // engaged interfaces have been set equal to 'i'. Now notice
            // the result as a candidate assembly and terminate the
            // branch if the result is a repetitive one

            if (NoteMultimer(r1)>=0)  {

              // Proceed to the next recursion level.
              //  a) count engageable interfaces
              k = 0;
              for (j=0;j<nETypes;j++)
                if (recStack[r1]->iflag[j]==INTF_Open)  k++;

              //  b) move up the recursion if any engageable interfaces
              //     left
              if (k>0)  EnumerateInterfaces ( r1,i+1 );
            }

            //  Roll back the whole branch (unmerge assemblies)

            RollBack ( r );

          }

        }

        if (startCPUclock>0)  {
          if (int(clock()-startCPUclock)/CLOCKS_PER_SEC>timeLimit)
            runStatus = RSTAT_TimeLimit;
        }

      }

    }

  }


  void Assembler::CalcEqIds()  {
  int i,j,k,m,id;

    for (i=0;i<nMultSets;i++)
      for (j=0;j<multSet[i]->nMultimers;j++)
        multSet[i]->U[j]->type = -1;

    id = 0;
    for (i=0;i<nMultSets;i++)
      for (j=0;j<multSet[i]->nMultimers;j++)
        if (multSet[i]->U[j]->type<0)  {
          multSet[i]->U[j]->type = id;
          for (k=0;k<nMultSets;k++)
            for (m=0;m<multSet[k]->nMultimers;m++)
              if (multSet[k]->U[m]->type<0)  {
                if (multSet[i]->U[j]->isIdentical(
                                        multSet[k]->U[m],D,rom,false))
                  multSet[k]->U[m]->type = id;
              }
          id++;
        }

  }


  int Assembler::NoteMultimer ( int r )  {
  PPMultimerSet US;
  int            i,k;

    //  1. Check that obtained set of multimers is a new one

    k = -1;
    for (i=0;(i<nMultSets) && (k<0);i++)
      if (multSet[i]->isEqual(recStack[r],nETypes))  k = i;

    //  2. Add the set to results if it is a new one

    if (k<0)  {

      if (nMultSets>=nMultSetAlloc)  {
        nMultSetAlloc += 100;
        US = new PMultimerSet[nMultSetAlloc];
        for (i=0;i<nMultSets;i++)
          US[i] = multSet[i];
        for (i=nMultSets;i<nMultSetAlloc;i++)
          US[i] = NULL;
        delete[] multSet;
        multSet = US;
      }

      if (!multSet[nMultSets])  multSet[nMultSets] = new MultimerSet();
      multSet[nMultSets]->makeSet ( recStack[r],nETypes,nInterfaces,
                                    nMonomers,interface,D,rom,eType );
      multSet[nMultSets]->CalcProperties ( D,interface,nInterfaces,
                                           eType,nETypes,&SNum,rom );
      if (multSet[nMultSets]->maxStableSize>maxStableSize)  {
        maxStableSize = multSet[nMultSets]->maxStableSize;
        noImproveCnt  = 0;
      } else  {
        noImproveCnt++;
        if ((noImproveCnt>noImproveLimit) && (!fullSearch))  {
          // this will terminate only current branch from the very top
          // of it; other branches will be attempted as usual
          runStatus = RSTAT_QuitLoop;
        }
      }

      k = nMultSets;
      nMultSets++;

    } else
      k = -(k+1);

    return k;

  }

  void Assembler::GetAssemblies ( RPAssemblies assemblies )  {
  int  i,j,serNo;

    if (assemblies)  {
      delete assemblies;
      assemblies = NULL;
    }

    if ((retcode==ASSMB_Ok) && (nMultSets>0))  {

      //  1. Make Assemblies from stable multimers

      assemblies = new Assemblies();
      assemblies->nCrystSplits = nMultSets;
      assemblies->nInterfaces  = nInterfaces0;
      assemblies->crystSplit       = new PCrystSplit[nMultSets];

      for (i=0;i<nMultSets;i++)  {
        assemblies->crystSplit[i] = new CrystSplit();
        assemblies->crystSplit[i]->Copy ( multSet[i],interface,
                                          nInterfaces,nInterfaces0,
                                          D->domain,rom );
      }

  //    assemblies->Sort             (); // now intentionally external
      assemblies->makeChainMapping ();
      assemblies->makeOrientations ( D,MMDB );
      assemblies->Orth2Frac        ( MMDB   );
      assemblies->checkOriginalOrientations ( D );

      //  2. Calculate assembly symmetry numbers

      serNo = 1;
      for (i=0;i<assemblies->nCrystSplits;i++)
        for (j=0;j<assemblies->crystSplit[i]->nAssemblies;j++)  {
          assemblies->crystSplit[i]->A[j]->symNumber =
            SNum.getSymNumber ( assemblies->crystSplit[i]->A[j],D );
          assemblies->crystSplit[i]->A[j]->serNo = serNo++;
        }

    }

  }

  void Assembler::GetComplex ( RPAssembly assembly )  {

    if (Complex)  {
      if (!assembly)  assembly = new Assembly();
      assembly->Copy ( Complex,interface,nInterfaces,nInterfaces0,
                       D->domain,rom );
      assembly->makeOrientation ( D,MMDB );
      assembly->symNumber = SNum.getSymNumber ( assembly,D );
      assembly->serNo     = 0;
    } else if (assembly)  {
      delete assembly;
      assembly = NULL;
    }
  }

}  // namespace pisa
