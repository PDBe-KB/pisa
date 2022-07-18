// $Id: pisa_domains.cpp $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_domains <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Domains
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2014
//
// =================================================================
//

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "pisa_domains.h"
#include "pisa_defs.h"

#include "mmdb2/mmdb_seqsuperpose.h"
#include "mmdb2/mmdb_math_linalg.h"
#include "mmdb2/mmdb_tables.h"

#include "ssm/ssm_superpose.h"
#include "ssm/ssm_csia.h"

namespace pisa  {

  //  =======================  Domains  ============================

  void DAlign::Init()  {
    Qscore = -1.0;
    rmsd   = -1.0;
    seqId  = -1.0;
    mmdb::Mat4Init ( TMatrix );
    Nalign = -1;
    equal  = false;
  }

  #define Qthreshold   0.65
  #define SIthreshold  0.90

  void DAlign::calcEqual()  {
    equal = (Qscore>Qthreshold) && (seqId>SIthreshold);
  }

  void DAlign::Copy ( RDAlign DA, bool direct )  {
    Qscore = DA.Qscore;
    rmsd   = DA.rmsd;
    seqId  = DA.seqId;
    Nalign = DA.Nalign;
    equal  = DA.equal;
    if (direct)  mmdb::Mat4Copy    ( DA.TMatrix,TMatrix );
           else  mmdb::Mat4Inverse ( DA.TMatrix,TMatrix );
  }

  void DAlign::write ( mmdb::io::RFile f )  {
  int i,j;
    f.WriteReal ( &Qscore );
    f.WriteReal ( &rmsd   );
    f.WriteReal ( &seqId  );
    for (i=0;i<3;i++)
      for (j=0;j<4;j++)
        f.WriteReal ( &(TMatrix[i][j]) );
    f.WriteInt  ( &Nalign );
    f.WriteBool ( &equal  );
  }

  void DAlign::read ( mmdb::io::RFile f )  {
  int i,j;
    f.ReadReal ( &Qscore );
    f.ReadReal ( &rmsd   );
    f.ReadReal ( &seqId  );
    for (i=0;i<3;i++)
      for (j=0;j<4;j++)
        f.ReadReal ( &(TMatrix[i][j]) );
    TMatrix[3][0] = 0.0;
    TMatrix[3][1] = 0.0;
    TMatrix[3][2] = 0.0;
    TMatrix[3][3] = 1.0;
    f.ReadInt  ( &Nalign );
    f.ReadBool ( &equal  );
  }


  //  -----------------------------------------------------------------

  Domains::Domains() : mmdb::io::Stream()  {
    InitDomains();
  }

  Domains::Domains ( mmdb::io::RPStream Object )
         : mmdb::io::Stream ( Object )  {
    InitDomains();
  }

  Domains::~Domains()  {
    FreeMemory  ();
    RemoveBricks();
  }

  void Domains::InitDomains()  {

    domain      = NULL;
    DA          = NULL;
    G           = NULL;
    Dmax        = 200.0;    // maximal domain radius for bricking
    nDomains    = 0;
    nDTypes     = 0;
    nNCSOps     = 0;
    nNCSParents = 0;

    DB         = NULL;

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

  void Domains::FreeMemory()  {
  int i;
    if (domain)  {
      for (i=0;i<nDomains;i++)
        if (domain[i])  delete domain[i];
      delete[] domain;
      domain = NULL;
    }
    FreeAligns();
    FreeGraphs();
    nDomains    = 0;
    nDTypes     = 0;
    nNCSOps     = 0;
    nNCSParents = 0;
  }

  void Domains::FreeGraphs()  {
  int i;
    if (G)  {
      for (i=0;i<nNCSParents;i++)
        if (G[i])  delete G[i];
      delete[] G;
      G = NULL;
    }
  }

  void Domains::FreeAligns()  {
  int i;
    if (DA)  {
      for (i=0;i<nNCSParents-1;i++)
        if (DA[i])  {
          DA[i] += (i+1);
          delete[] DA[i];
        }
      delete[] DA;
      DA = NULL;
    }
  }


  bool Mat4Compare ( mmdb::mat44 & t1, mmdb::mat44 & t2 )  {
  mmdb::realtype d,dr,dt;
  int      i;
    dr = 0.0;
    dt = 0.0;
    for (i=0;i<3;i++)  {
      d   = t1[i][0] - t2[i][0];
      dr += d*d;
      d   = t1[i][1] - t2[i][1];
      dr += d*d;
      d   = t1[i][2] - t2[i][2];
      dr += d*d;
      d   = t1[i][3] - t2[i][3];
      dt += d*d;
    }
    return (dr<0.0003) && (dt<1.0);
  }


  bool isMat4Rot ( mmdb::mat44 & t, mmdb::realtype eps )  {
  mmdb::realtype a;
  int      i,j,k;
  bool  r;
    r = true;
    for (i=0;(i<3) && r;i++)
      for (j=0;(j<3) && r;j++)  {
        a = 0.0;
        for (k=0;k<3;k++)
          a += t[i][k]*t[j][k];
        if (i!=j)  r = (fabs(a)<=eps);
             else  r = (fabs(1.0-a)<=eps);
      }
    if (r)  {
      a = t[0][0]*t[1][1]*t[2][2] +
          t[0][2]*t[1][0]*t[2][1] +
          t[2][0]*t[0][1]*t[1][2] -
          t[0][2]*t[1][1]*t[2][0] -
          t[0][0]*t[1][2]*t[2][1] -
          t[2][2]*t[1][0]*t[0][1];
      r = (fabs(1.0-a)<=eps);
    }
    return r;
  }


  void makeResWarning ( mmdb::pstr & warnings, mmdb::pstr resName )  {
  mmdb::pstr p,q;
    if (!warnings)
      mmdb::CreateCopy ( warnings,resName );
    else {
      p = NULL;
      mmdb::CreateCopCat ( p,"@",warnings,"@" );
      q = p;
      while (*q)  {
        if (*q==' ')  *q = '@';
        q++;
      }
      q = NULL;
      mmdb::CreateCopCat ( q,"@",resName,"@" );
      if (!strstr(p,q))  mmdb::CreateConcat ( warnings," ",resName );
      if (p)  delete[] p;
      if (q)  delete[] q;
    }
  }

  void addLigandRes ( mmdb::pstr & ligandRes, mmdb::pstr resName )  {
  mmdb::pstr p,q;
    if (!ligandRes)
      mmdb::CreateCopy ( ligandRes,resName );
    else {
      p = NULL;
      mmdb::CreateCopCat ( p,",",ligandRes,"," );
      q = NULL;
      mmdb::CreateCopCat ( q,",",resName,"," );
      if (!strstr(p,q))  mmdb::CreateConcat ( ligandRes,",",resName );
      if (p)  delete[] p;
      if (q)  delete[] q;
    }
  }

  void Domains::allocateDomains ( int & nAlloc )  {
  PPDomain D1;
  int         i;
    if (nDomains>=nAlloc)  {
      nAlloc = nDomains + 100;
      D1 = new PDomain[nAlloc];
      for (i=0;i<nDomains;i++)
        D1[i] = domain[i];
      for (i=nDomains;i<nAlloc;i++)
        D1[i] = NULL;
      if (domain)  delete[] domain;
      domain = D1;
    }
  }


  //  Maximum tracktable number of ligands in unit cell.
  static int maxNofLigandsUC  = 750;
  static int maxNofLigandsASU = 500;

  void SetMaxNofLigands ( int maxNLigandsUC, int maxNLigandsASU )  {
    maxNofLigandsUC  = maxNLigandsUC;
    maxNofLigandsASU = maxNLigandsASU;
  }


  void Domains::MakeDomains ( mmdb::PManager MMDB,
                              PMolRefIndex   molRef,
                              mmdb::cpstr    agentsRef,
                              mmdb::pstr   & warnUnk,
                              mmdb::pstr   & warnExcl )  {
  //   Create list of domains, that is, each AA or NA chain is a domain.
  // If a ligand comes with AA or NA chain ID, consider that it makes
  // a covalent bonding with the chain (perhaps wrong), otherwise
  // each ligand residue makes a domain.
  mmdb::PPChain   chain;
  mmdb::PPResidue res;
  PMolRefEntry    MRE;
  mmdb::mat44     ncs_m,tm,TMatrix;;
  mmdb::pstr      ligRes;
  char            S[100];
  mmdb::realtype  eps, mjx,mjy,mjz, dx,dy,dz;
  int             i,j,nOfChains,nRes,iGiven;
  int             nAAAtoms,nAAResidues,nNAAtoms,nNAResidues,nRNARes;
  int             ic,jc,kc,m,k,nAlloc, nSymOps,nNCSMates;
  int             nLigands,symFactor;

    FreeMemory();

    MMDB->GetChainTable ( 1,chain,nOfChains );

    ligRes   = NULL;
    nLigands = 0;

    nAlloc   = 0;
    for (i=0;i<nOfChains;i++)
      if (chain[i])  {
        nAAAtoms    = 0;
        nAAResidues = 0;
        nNAAtoms    = 0;
        nNAResidues = 0;
        nRNARes     = 0;
        chain[i]->GetResidueTable ( res,nRes );
        for (j=0;j<nRes;j++)
          if (res[j])  {
            MRE = molRef->getMolRefEntry ( res[j]->name );
            if (!MRE)
              makeResWarning ( warnUnk,res[j]->name );
            else  {
              switch (MRE->classId)  {

                case MRCLASS_Aminoacid     :
                case MRCLASS_ModAminoacid  :
                           nAAAtoms += res[j]->nAtoms;
                           nAAResidues++;
                         break;

                case MRCLASS_Solvent       : break;

                case MRCLASS_Nucleotide    :
                case MRCLASS_ModNucleotide :
                           nNAAtoms += res[j]->nAtoms;
                           if (res[j]->isDNARNA()==2)  nRNARes++;
                           nNAResidues++;
                         break;

                case MRCLASS_Ligand        :
                           allocateDomains ( nAlloc );
                           addLigandRes    ( ligRes  ,res[j]->name );
                           domain[nDomains] = new Domain();
                           domain[nDomains]->SetLigandRange ( res[j] );
                           domain[nDomains]->ncsParent = nDomains;
                           domain[nDomains]->dclass    = DCLASS_Ligand;
                           domain[nDomains]->symNumber = MRE->symNumber;
                           if (agentsRef)  {
                             sprintf ( S,",%s,",res[j]->name );
                             domain[nDomains]->agent =
                                          (strstr(agentsRef,S)!=NULL);
//                             domain[nDomains]->assemble
                           }
                           nDomains++;
                           nLigands++;
                         break;

                default                    :
                case MRCLASS_Unknown       :
                           makeResWarning ( warnUnk,res[j]->name );

              }
            }
          }

        if (nAAResidues>nNAResidues)  {
          //  make a protein domain
          if (nAAAtoms>=5*nAAResidues)  {
            // complete aminoacid chain, take it as a domain
            allocateDomains ( nAlloc );
            domain[nDomains] = new Domain();
            domain[nDomains]->SetChainRange ( chain[i]->GetChainID() );
            domain[nDomains]->ncsParent = nDomains;
            domain[nDomains]->dclass    = DCLASS_Protein;
            domain[nDomains]->symNumber = 1;
            domain[nDomains]->SetLigandResidues ( ligRes );
            nDomains++;
          }
        } else if (nNAResidues>nAAResidues)  {
          // make a RNA/DNA domain
          if (nNAAtoms>=17*nNAResidues)  {
            // complete nucleotide chain, take it as a domain
            allocateDomains ( nAlloc );
            domain[nDomains] = new Domain();
            domain[nDomains]->SetChainRange ( chain[i]->GetChainID() );
            domain[nDomains]->ncsParent = nDomains;
            if (nRNARes>=nNAResidues-nRNARes)
                  domain[nDomains]->dclass = DCLASS_RNA;
            else  domain[nDomains]->dclass = DCLASS_DNA;
            domain[nDomains]->symNumber = 1;
            domain[nDomains]->SetLigandResidues ( ligRes );
            nDomains++;
          }
        }

        if (ligRes)  {
          delete[] ligRes;
          ligRes = NULL;
        }

      }

    //  Generate NCS-mates

    nNCSOps   = MMDB->GetNumberOfNCSMatrices();
    nNCSMates = MMDB->GetNumberOfNCSMates   ();
    nSymOps   = MMDB->GetNumberOfSymOps     ();

    if (nSymOps>0)  {
      symFactor = (1+nNCSMates)*nSymOps;
      if (nLigands*symFactor>maxNofLigandsUC)
        excludeLigands ( symFactor,maxNofLigandsUC,warnExcl );
    } else  {
      symFactor = 1+nNCSMates;
      if (nLigands*symFactor>maxNofLigandsASU)
        excludeLigands ( symFactor,maxNofLigandsASU,warnExcl );
    }

    nNCSParents = nDomains;

    for (m=0;(m<nSymOps) && (nNCSMates>0);m++)
      if (MMDB->GetTMatrix(tm,m,0,0,0)==mmdb::SYMOP_Ok)  {
        if (!isMat4Rot(tm,1.0e-2))  nNCSMates = -1;
      } else
        nNCSMates = -1;

    if (nNCSMates>0)  {

      //   Checking on the number of NCS-mates (NCS operations for
      // which there is no coordinates in the file, iGiven==0), allows
      // to avoid this branch and save some CPU on the web-server

      eps     = 100.0*mmdb::floatMachEps;
      calcDimensions ( MMDB );

      for (i=0;i<nNCSOps;i++)
        if (MMDB->GetNCSMatrix(i,ncs_m,iGiven))  {
          if ((!mmdb::isMat4Unit(ncs_m,eps,false)) &&
                isMat4Rot(ncs_m,1.0e-2))  {
            // we don't believe in iGiven as it comes from the PDB
            for (j=0;j<nNCSParents;j++)  {
              // check that ith NCS mate of jth NCS parent
              // does not clash with already existing
              // NCS mates and parents
              iGiven = 0;
              for (m=0;(m<nSymOps) && (!iGiven);m++)
                for (ic=-2;(ic<=2) && (!iGiven);ic++)
                  for (jc=-2;(jc<=2) && (!iGiven);jc++)
                    for (kc=-2;(kc<=2) && (!iGiven);kc++)
                      if (MMDB->GetTMatrix(TMatrix,m,ic,jc,kc)
                                                   ==mmdb::SYMOP_Ok)  {
                        mmdb::Mat4Mult ( tm,TMatrix,ncs_m );
                        mjx = domain[j]->mx;
                        mjy = domain[j]->my;
                        mjz = domain[j]->mz;
                        mmdb::TransformXYZ ( tm,mjx,mjy,mjz );
                        for (k=0;(k<nDomains) && (!iGiven);k++)  {
                          dx = domain[k]->mx - mjx;
                          if (fabs(dx)<1.0)  {
                            dy = domain[k]->my - mjy;
                            if (fabs(dy)<1.0)  {
                              dz = domain[k]->mz - mjz;
                              if (dx*dx+dy*dy+dz*dz<1.0)  iGiven = 1;
                            }
                          }
                       }
                      }
              if (!iGiven)  {
                // there's no clash, allocate domain
                allocateDomains ( nAlloc );
                domain[nDomains] = new Domain();
                domain[nDomains]->ncsOpNo   = i+1;
                domain[nDomains]->ncsParent = j;
                mmdb::Mat4Copy ( ncs_m,domain[nDomains]->ncs_m );
                domain[nDomains]->CopyParentData(domain);
                nDomains++;
              }
            }
          }
        }

    }


  /*
    if (nNCSMates>0)  {
      //   Checking on the number of NCS-mates (NCS operations for
      // which there is no coordinates in the file, iGiven==0), allows
      // to avoid this branch and save some CPU on the web-server
      eps     = 100.0*floatMachEps;
      calcDimensions ( MMDB );
      for (i=0;i<nNCSOps;i++)
        if (MMDB->GetNCSMatrix(i,ncs_m,iGiven))  {
          if ((!iGiven) && (!mmdb::isMat4Unit(ncs_m,eps,false)) &&
              isMat4Rot(ncs_m,1.0e-2))  {
            for (m=0;(m<nSymOps) && (!iGiven);m++)
              if (MMDB->GetTMatrix(TMatrix,m,0,0,0)==mmdb::SYMOP_Ok)  {
                if (isMat4Rot(TMatrix,1.0e-2))  {
                  for (ic=-2;(ic<=2) && (!iGiven);ic++)
                    for (jc=-2;(jc<=2) && (!iGiven);jc++)
                      for (kc=-2;(kc<=2) && (!iGiven);kc++)
                        if (MMDB->GetTMatrix(TMatrix,m,ic,jc,kc)
                                                     ==mmdb::SYMOP_Ok)  {
                          mmdb::Mat4Mult ( tm,TMatrix,ncs_m );
                          for (j=1;(j<nNCSParents) && (!iGiven);j++)
                            for (k=0;(k<nDomains) && (!iGiven);k++)  {
                              dx = tm[0][0]*domain[j]->mx +
                                   tm[0][1]*domain[j]->my +
                                   tm[0][2]*domain[j]->mz +
                                   tm[0][3] - domain[k]->mx;
                              if (fabs(dx)<4.0)  {
                                dy = tm[1][0]*domain[j]->mx +
                                     tm[1][1]*domain[j]->my +
                                     tm[1][2]*domain[j]->mz +
                                     tm[1][3] - domain[k]->my;
                                if (fabs(dy)<4.0)  {
                                  dz = tm[2][0]*domain[j]->mx +
                                       tm[2][1]*domain[j]->my +
                                       tm[2][2]*domain[j]->mz +
                                       tm[2][3] - domain[k]->mz;
                                  if (dx*dx+dy*dy+dz*dz<16.0) iGiven = 1;
                                }
                              }
                            }
                        }
                }
              }
            if (!iGiven)
              for (j=0;j<nNCSParents;j++)  {
                allocateDomains ( nAlloc );
                domain[nDomains] = new Domain();
                domain[nDomains]->ncsOpNo   = i+1;
                domain[nDomains]->ncsParent = j;
                mmdb::Mat4Copy ( ncs_m,domain[nDomains]->ncs_m );
                domain[nDomains]->CopyParentData ( Domain );
                nDomains++;
              }
          }
        }

    }
  */

  }

  void Domains::excludeLigands ( int symFactor, int maxNofLigands,
                                  mmdb::pstr & warnExcl ) {
  // This function excludes some ligands if their total number is
  // too high. The function should be called only from MakeDomains(..).
  mmdb::ivector nLig;
  char    S[100];
  int     i,j,k,n;

    mmdb::GetVectorMemory ( nLig,nDomains,0 );

    for (i=0;i<nDomains;i++)
      nLig[i] = 0;

    n = 0;
    for (i=0;i<nDomains;i++)
      if ((domain[i]->dclass==DCLASS_Ligand) && (!nLig[i]))  {
        nLig[i] = 1;
        for (j=i+1;j<nDomains;j++)
          if (!nLig[j])  {
            if (domain[i]->checkLigandRange(domain[j]))  {
              nLig[i]++;
              nLig[j] = -(i+1);
            }
          }
        nLig[i] *= symFactor;
        n += nLig[i];
      }

    while (n>maxNofLigands)  {
      j = -1;
      k = 0;
      for (i=0;i<nDomains;i++)
        if (nLig[i]>k)  {
          k = nLig[i];
          j = i;
        }
      if (j>=0)  {
        makeResWarning ( warnExcl,domain[j]->getLigandRange(S) );
        n -= nLig[j];
        delete domain[j];
        domain[j] = NULL;
        nLig[j]   = 0;
        for (i=j+1;i<nDomains;i++)
          if (nLig[i]==-(j+1))  {
            delete domain[i];
            domain[i] = NULL;
            nLig[i]   = 0;
          }
      } else
        n = -1;  // abnormal quit
    }

    k = 0;
    for (i=0;i<nDomains;i++)
      if (domain[i])  {
        if (i>k)  {
          domain[k] = domain[i];
          domain[i] = NULL;
          domain[k]->ncsParent = k;
        }
        k++;
      }

    nDomains = k;

    mmdb::FreeVectorMemory ( nLig,0 );

  }


  void Domains::getNofDomains ( int & nAADomains, int & nNADomains,
                                int & nLigands ) {
  int i;
    nAADomains = 0;
    nNADomains = 0;
    nLigands   = 0;
    for (i=0;i<nDomains;i++)
      if (domain[i])  {
        switch (domain[i]->dclass)  {
          case DCLASS_Protein :  nAADomains++; break;
          case DCLASS_DNA     :  nNADomains++; break;
          case DCLASS_RNA     :  nNADomains++; break;
          case DCLASS_Ligand  :  nLigands++;   break;
          default : ;
        }
      }
  }


  void Domains::calcDimensions ( mmdb::PManager MMDB )  {
  int i;

    Dmax = 0.0;
    for (i=0;i<nNCSParents;i++)
      if (domain[i])  {
        domain[i]->calcDimensions ( MMDB );
        if (domain[i]->R>Dmax)  Dmax = domain[i]->R;
      }

    randomizeReferenceFrames();

    for (i=nNCSParents;i<nDomains;i++)
      if (domain[i])
        domain[i]->CopyParentData ( domain );


  }

  PROSURF_RC Domains::CalcSurfaces ( mmdb::PManager MMDB,
                                     PProSurf       proSurf,
                                     PMolRefIndex   molRef )  {
  // Calculate domain surfaces for the identification of
  // potential interfaces. Domain dimension calculations
  // must be done before calling this function
  int        i;
  PROSURF_RC rc,rc1;

    Dmax = 0.0;
    rc   = PROSURF_Ok;
    for (i=0;i<nNCSParents;i++)
      if (domain[i])  {
        rc1 = domain[i]->calcSurface ( MMDB,proSurf,molRef );
        if (rc1!=PROSURF_Ok)  rc = rc1;
        if (domain[i]->R>Dmax)  Dmax = domain[i]->R;
      }

    if (rc==PROSURF_Ok)
      for (i=nNCSParents;i<nDomains;i++)
        if (domain[i])
          domain[i]->CopyParentData ( domain );

    return rc;

  }


  void Domains::makeSSGraphs ( mmdb::PManager MMDB,
                                mmdb::ivector & selHnd )  {
  //  Makes SS graphs for each domain. The function does not
  //  deallocate vector GG.
  ssm::PGraph GG;
  int         i,rc;

    FreeGraphs();
    mmdb::FreeVectorMemory ( selHnd,0 );

    GG = new ssm::Graph();
    rc = GG->MakeGraph ( MMDB );
    if (!rc)  {
      G = new ssm::PGraph[nNCSParents];
      mmdb::GetVectorMemory ( selHnd,nNCSParents,0 );
      for (i=0;i<nNCSParents;i++)  {
        G[i] = NULL;
        if (domain[i])  {
          if (domain[i]->dclass==DCLASS_Protein)  {
            selHnd[i] = domain[i]->SelectDomain(MMDB,mmdb::STYPE_ATOM,1,false);
            if (selHnd[i]>0)  {
              G[i] = new ssm::Graph();
              G[i]->Copy ( GG );
              G[i]->LeaveVertices ( selHnd[i],MMDB );
              G[i]->calcVTypes();
            }
          }
        }
      }
    }

    delete GG;

  }

  void Domains::completeSSGraphs ( mmdb::PManager MMDB,
                                      mmdb::ivector & selHnd )  {

  ssm::PGraph GG;
  int         i;

    if (!G)  {
      makeSSGraphs ( MMDB,selHnd );
      return;
    }

    GG = NULL;
    if (!selHnd)  {
      mmdb::GetVectorMemory ( selHnd,nNCSParents,0 );
      for (i=0;i<nNCSParents;i++)
        selHnd[i] = 0;
    }

    for (i=0;i<nNCSParents;i++)
      if (!selHnd[i])
        selHnd[i] = domain[i]->SelectDomain ( MMDB,mmdb::STYPE_ATOM,1,false );

    for (i=0;i<nNCSParents;i++)
      if (domain[i])  {
        if (domain[i]->dclass==DCLASS_Protein)  {
          if ((selHnd[i]>0) && (!G[i]))  {
            if (!GG)  {
              GG = new ssm::Graph();
              if (GG->MakeGraph(MMDB)) break;
            }
            G[i] = new ssm::Graph();
            G[i]->Copy ( GG );
            G[i]->LeaveVertices ( selHnd[i],MMDB );
            G[i]->calcVTypes();
          }
        }
      }

    if (GG)  delete GG;

  }


  void Domains::AlignSimilar ( mmdb::PManager MMDB )  {
  //  Alignes similar domains to each other using SeqSuperpose.
  // This function works faster than SSM for similar chains,
  // at the same quality of alignment. Performance is crucial
  // as similar chains need to be identified on the front-end
  // machine. AlignSimilar(..) destroys any previously existing
  // alignments.
  mmdb::SeqSuperpose SeqSup;
  mmdb::PPAtom       a1,a2;
  mmdb::ivector      selHnd,C;
  mmdb::mat44        TM;
  mmdb::vect3        v,vi;
  mmdb::realtype     rmsd,x,y,z;
  int                i,j,k,m,n,na1,na2,rc;

    FreeAligns();

    if (nNCSParents<=1)  {
      if (nNCSParents==1)  {
        domain[0]->type = 1;
        nDTypes = 1;
      } else
        nDTypes = 0;
      return;
    }

    //  1. Select all non-ligand domains

    mmdb::GetVectorMemory ( selHnd,nNCSParents,0 );
    for (i=0;i<nNCSParents;i++)  {
      selHnd[i] = 0;
      if (domain[i])  {
        if (domain[i]->dclass!=DCLASS_Ligand)  {
          selHnd[i] = domain[i]->SelectDomain(MMDB,mmdb::STYPE_ATOM,1,false);
          if ((domain[i]->dclass==DCLASS_DNA) ||
              (domain[i]->dclass==DCLASS_RNA))
            // this would work also with phosphorus
            MMDB->Select ( selHnd[i],mmdb::STYPE_ATOM,0,"*",
                           mmdb::ANY_RES,"*",mmdb::ANY_RES,"*","*","C5*,C5'","C","*",
                           mmdb::SKEY_AND );
          else if (domain[i]->dclass==DCLASS_Protein)
            MMDB->Select ( selHnd[i],mmdb::STYPE_ATOM,0,"*",
                           mmdb::ANY_RES,"*",mmdb::ANY_RES,"*","*","CA","C","*",
                           mmdb::SKEY_AND );
        }
      }
    }

    //  2. Make all-against-all 3D alignments

    AllocateAlignments();

    //  2.1 Align similar aminoacid and nucleotide chains

    for (i=0;i<nNCSParents-1;i++)
      if (selHnd[i]>0)  {
        MMDB->GetSelIndex ( selHnd[i],a1,na1 );
        if (na1>0)  {
          MMDB->RemoveBricks();  // will be rebuilt by SeqSuperpose once
                                 // for all comparisons in the next loop
          for (j=i+1;j<nNCSParents;j++)  {
            if (selHnd[j]>0)  {
              if (domain[i]->dclass==domain[j]->dclass)  {
                MMDB->GetSelIndex ( selHnd[j],a2,na2 );
                if (na2>0)  {
                  rc = SeqSup.Superpose ( MMDB,a2,na2,a1,na1,0.8,true );
                  if (rc==mmdb::SEQSP_Ok)  {
                    mmdb::Mat4Inverse ( SeqSup.TMatrix,DA[i][j].TMatrix );
                    DA[i][j].rmsd   = SeqSup.rmsd;
                    DA[i][j].Nalign = SeqSup.Nalign;
                    DA[i][j].seqId  = SeqSup.seqId;
                    DA[i][j].Qscore = SeqSup.Q;
                  }
                }
              }
              DA[i][j].calcEqual();
            }
          }
        }
      }

    MMDB->RemoveBricks();

    //  2.2 Align ligands

    C = NULL;
    n = 0;

    for (i=0;i<nNCSParents;i++)  {
      if (selHnd[i]>0)  MMDB->DeleteSelection ( selHnd[i] );
      selHnd[i] = 0;
      if (domain[i])  {
        if (domain[i]->dclass==DCLASS_Ligand)  {
          selHnd[i] = domain[i]->SelectDomain(MMDB,mmdb::STYPE_ATOM,1,false);
          n = mmdb::IMax ( n,MMDB->GetSelLength(selHnd[i]) );
        }
      }
    }

    if (n>0)  {

      mmdb::GetVectorMemory ( C,n,0 );

      for (i=0;i<nNCSParents-1;i++)
        if (selHnd[i])  {
          MMDB->GetSelIndex ( selHnd[i],a1,na1 );
          if (na1>=3)  {
            // superpose all atoms for larger ligands
            for (j=i+1;j<nNCSParents;j++)
              if (selHnd[j])  {
                if (domain[i]->checkLigandRange(domain[j]))  {
                  MMDB->GetSelIndex ( selHnd[j],a2,na2 );
                  if (na2>0)  {
                    // make the corresponding vector C in anticipation
                    // that ligand descriptions may miss atoms
                    for (k=0;k<na1;k++)  {
                      C[k] = -1;
                      for (m=0;(m<na2) && (C[k]<0);m++)
                        if (!strcmp(a1[k]->name,a2[m]->name))
                          C[k] = m;
                    }
                    rc = SuperposeAtoms ( TM,a1,na1,a2,C );
                    if (rc==mmdb::SPOSEAT_Ok)  {
                      rmsd = 0.0;
                      n    = 0;
                      for (k=0;k<na1;k++)
                        if (C[k]>=0)  {
                          n++;  // alignment length
                          a1[k]->TransformCopy ( TM,x,y,z );
                          rmsd += a2[C[k]]->GetDist2 ( x,y,z );
                        }
                      mmdb::Mat4Copy ( TM,DA[i][j].TMatrix );
                      if (n>0) DA[i][j].rmsd = sqrt(rmsd/n);
                          else DA[i][j].rmsd = 0.0;
                      DA[i][j].Nalign = n;
                      DA[i][j].seqId  = 1.0;
                      DA[i][j].Qscore = n*n /
                                        ((1.0+rmsd/(n*9.0))*na1*na2);
                    }
                  }
                  DA[i][j].calcEqual();
                }
              }
          } else  {
            // small ligands are aligned by their principal axes
            vi[0] = domain[i]->mx;
            vi[1] = domain[i]->my;
            vi[2] = domain[i]->mz;
            for (j=i+1;j<nNCSParents;j++)
              if (selHnd[j])  {
                if (domain[i]->checkLigandRange(domain[j]))  {
                  for (k=0;k<3;k++) {
                    v[k] = 0.0;
                    for (m=0;m<3;m++)  {
                      TM[k][m] = 0.0;
                      for (n=0;n<3;n++)
                        TM[k][m] +=
                              domain[j]->itb[k][n]*domain[i]->itb[m][n];
                        v[k] += TM[k][m]*vi[m];
                    }
                    TM[3][k] = 0.0;
                  }
                  TM[3][3] = 1.0;
                  TM[0][3] = domain[j]->mx - v[0];
                  TM[1][3] = domain[j]->my - v[1];
                  TM[2][3] = domain[j]->mz - v[2];
                  mmdb::Mat4Copy ( TM,DA[i][j].TMatrix );
                  DA[i][j].rmsd   = 0.0;
                  DA[i][j].Nalign = domain[i]->nAtoms;
                  DA[i][j].seqId  = 1.0;
                  DA[i][j].Qscore = 1.0;
                  DA[i][j].calcEqual();
                }
              }
          }
        }

      mmdb::FreeVectorMemory ( C,0 );

    }


    //  3. Calculate domain types

    for (i=0;i<nNCSParents;i++)  {
      if (selHnd[i]>0)  MMDB->DeleteSelection ( selHnd[i] );
      domain[i]->type = 0;
    }

    mmdb::FreeVectorMemory ( selHnd,0 );

    CalcTypes();

  }


  int Domains::getNnotAligned()  {
  //  This function returns the total number of yet not aligned
  // pairs of domains. The function accounts for all previous
  // alignments but neglects non-aligneable pairs of domains
  // without secondary structure.
  int i,j,n;

    if (!DA)  return nNCSParents*(nNCSParents-1)/2;

    n = 0;
    for (i=0;i<nNCSParents;i++)
      if (!DA[i])
        n += nNCSParents - i - 1;
      else
        for (j=i+1;j<nNCSParents;j++)
          if (DA[i][j].Qscore<0.0)
            n++;

    return n;

  }

  void Domains::AllocateAlignments()  {
  //   If alignment matrix DA does not exist, then it is
  // allocated and initialized in full. If DA exists then
  // the function checks for unallocated rows, and if such
  // are found then it allocates and initialize them. No
  // pre-existing rows and alignments are destroyed.
  int i,j;

    if (!DA)  {
      DA = new PDAlign[nNCSParents-1];
      for (i=0;i<nNCSParents-1;i++)
        DA[i] = NULL;
    }

    for (i=0;i<nNCSParents-1;i++)
      if (!DA[i]) {
        DA[i]  = new DAlign[nNCSParents-i-1];
        DA[i] -= (i+1);
        for (j=i+1;j<nNCSParents;j++)
          DA[i][j].Init();
      }

  }

  void Domains::writeAlignments ( mmdb::io::RFile f )  {
  //  This function writes alignments DA[i][j] for which
  // TMatrix[3][0]=mmdb::MaxReal, and clears this flag.
  int i,j;
    if (DA)
      for (i=0;i<nNCSParents;i++)
        if (DA[i])
          for (j=i+1;j<nNCSParents;j++)
            if (DA[i][j].TMatrix[3][0]==mmdb::MaxReal)  {
              DA[i][j].TMatrix[3][0] = 0.0;
              f.WriteInt ( &i );
              f.WriteInt ( &j );
              DA[i][j].write ( f );
            }
  }

  void Domains::AlignAll ( mmdb::PManager MMDB,
                              int startNo, int endNo )  {
  //   Alignes all domains to each other using SSM. This function
  // makes only alignments that were not done before (e.g. by
  // calling AlignSimilar(..)). If all alignments should be done,
  // call FreeAligns() right before this function and prior
  // getNnotAligned() if alignments are to be spread over a number
  // of CPUs (cf. below).
  //   startNo and endNo (inclusive) may be used in order to spread
  // alignment calculations over a number of CPUs. Use function
  // getNnotAligned() to get the total number of yet not aligned
  // pairs of domains, then define startNo and endNo for each CPU.
  // After execution of AlignAll(MMDB,startNo,endNo) in each CPU,
  // write alignments into result files using function
  // writeAlignments(f). See writeAlignments(..) for how to retrieve
  // alignments in the result collector. Note that, unless
  // endNo<startNo (which is an input for making *all* alignments),
  // the function sets DA[i][j].TMatrix[3][0]=mmdb::MaxReal to indicate
  // that the alignment was done or attempted. writeAlignment(f)
  // clears this flag, but otherwise it should be removed by
  // application.
  //   This function is not accomplished by calculating unique
  // domain types. Call CalcTypes() explicitely for doing that.
  ssm::GraphMatch U;
  ssm::Superpose  Superpose;
  ssm::PPMatch    Match;
  mmdb::ivector         F1,F2;
  mmdb::ivector         selHnd;
  mmdb::rvector         dist1;
  mmdb::ivector         Ca1,Ca2;
  mmdb::realtype        Q,Q1,ncombs;
  int             nres1,nres2,ngaps,nmd;
  int             i,j,k,nMatches,nm, n1,n2,n;

    //  1. Make all structural graphs

    selHnd = NULL;
    completeSSGraphs ( MMDB,selHnd );
    if (!G)  {
      if (selHnd)  {
        for (i=0;i<nNCSParents;i++)
          if (selHnd[i]>0)  MMDB->DeleteSelection ( selHnd[i] );
        mmdb::FreeVectorMemory ( selHnd,0 );
      }
      return;
    }

    //  check whether the alignments should be done at all

    if (nNCSParents<=1)  {
      if (selHnd)  {
        for (i=0;i<nNCSParents;i++)
          if (selHnd[i]>0)  MMDB->DeleteSelection ( selHnd[i] );
        mmdb::FreeVectorMemory ( selHnd,0 );
      }
      if (nNCSParents==1)  {
        domain[0]->type = 1;
        nDTypes = 1;
      } else
        nDTypes = 0;
      return;
    }

    //  build all graphs

    for (i=0;i<nNCSParents;i++)
      if (G[i])  {
        if (!G[i]->isBuild())
          G[i]->BuildGraph();
      }

    //  2. Make all-against-all 3D alignments

    ssm::SetMatchPrecision    ( ssm::PREC_Normal      );
    ssm::SetConnectivityCheck ( ssm::CONNECT_Flexible );

    U.SetUniqueMatch ( true );
    U.SetBestMatch   ( true );

    Ca1   = NULL;
    Ca2   = NULL;
    dist1 = NULL;

    AllocateAlignments();

    if (endNo<startNo)  {
      n1 = 0;
      n2 = mmdb::MaxInt;
    } else  {
      n1 = startNo;
      n2 = endNo;
    }
    n = 0;  // pair counter

    for (i=0;(i<nNCSParents-1) && (n<=n2);i++)  {
      if (G[i])  {
        for (j=i+1;(j<nNCSParents) && (n<=n2);j++)
          if (DA[i][j].Qscore<-0.6)  {
            if (n>=n1)  {
              // alignment was not done before
              Q = -0.5;  // Q-score of alignment
              if (G[j])  {
                // low similarity matching takes CPU but gives
                // no clue at all; limit it to 33% here
                nm = mmdb::IMax(1,mmdb::IMin(G[i]->GetNofVertices(),
                                             G[j]->GetNofVertices())/3 );
                U.MatchGraphs ( G[i],G[j],nm );
                U.GetMatches  ( Match,nMatches );
                for (k=0;k<nMatches;k++)
                  if (Match[k])  {
                    Match[k]->GetMatch ( F1,F2,nm );
                    Superpose.SuperposeCalphas ( G[i],G[j],F1,F2,nm,
                                                 MMDB,MMDB,selHnd[i],
                                                 selHnd[j] );
                    Q1 = Superpose.GetCalphaQ();
                    if ((Q1>0.0) && (Q1>Q))  {
                      Q = Q1;
                      Superpose.GetSuperposition ( Ca1,dist1,nres1,
                                                   Ca2,nres2,
                                                   DA[i][j].TMatrix,
                                                   DA[i][j].rmsd,
                                                   DA[i][j].Nalign,ngaps,
                                                   DA[i][j].seqId,nmd,
                                                   ncombs );
                    }
                  }
              }
              DA[i][j].Qscore = Q;
              DA[i][j].calcEqual();
              if (n2<mmdb::MaxInt)  DA[i][j].TMatrix[3][0] = mmdb::MaxReal;
            }
            n++;
          }
      } else
        n += nNCSParents - i - 1;
    }

    mmdb::FreeVectorMemory ( Ca1  ,0 );
    mmdb::FreeVectorMemory ( Ca2  ,0 );
    mmdb::FreeVectorMemory ( dist1,0 );

    if (selHnd)  {
      for (i=0;i<nNCSParents;i++)
        if (selHnd[i]>0)  MMDB->DeleteSelection ( selHnd[i] );
      mmdb::FreeVectorMemory ( selHnd,0 );
    }

  }


  void Domains::CalcTypes()  {
  mmdb::ResName  mmTID,ligTID;
  int            i,j,k, i1,j1,ki,kj;
  bool           changed;

    //  1. Make closed groups of equivalent structures, so that
    //     if A~B and B~C then *make* A~C if it is not already so

    do {
      changed = false;
      for (i=0;i<nNCSParents;i++)
        for (j=i+1;j<nNCSParents;j++)
          if (DA[i][j].equal)
            for (k=0;k<nNCSParents;k++)
              if ((k!=i) && (k!=j))  {
                if (i<k)  {
                  i1 = i;
                  ki = k;
                } else  {
                  i1 = k;
                  ki = i;
                }
                if (j<k)  {
                  j1 = j;
                  kj = k;
                } else  {
                  j1 = k;
                  kj = j;
                }
                if (DA[i1][ki].equal && (!DA[j1][kj].equal))  {
                  changed          = true;
                  DA[j1][kj].equal = true;
                } else if ((!DA[i1][ki].equal) && DA[j1][kj].equal)  {
                  changed          = true;
                  DA[i1][ki].equal = true;
                }
              }
    } while (changed);


    //  2. Renumber all domain types as 1,2..nDTypes

    strcpy ( mmTID ,"A" );
    strcpy ( ligTID,"a" );

    nDTypes = 0;
    for (i=0;i<nNCSParents;i++)
      if (domain[i]->type<=0)  {
        nDTypes++;
        domain[i]->type = nDTypes;
        if (domain[i]->dclass==DCLASS_Ligand)  {
          mmdb::CreateCopy ( domain[i]->typeId,ligTID );
          ligTID[0]++;
        } else  {
          mmdb::CreateCopy ( domain[i]->typeId,mmTID );
          mmTID[0]++;
        }
        for (j=0;j<nNCSParents;j++)
          if (domain[j]->type<=0)  {
            if (isEquivalent(i,j))  {
              domain[j]->type = nDTypes;
              mmdb::CreateCopy ( domain[j]->typeId,domain[i]->typeId );
            }
          }
      }

    for (i=nNCSParents;i<nDomains;i++)  {
      domain[i]->type = domain[domain[i]->ncsParent]->type;
      mmdb::CreateCopy ( domain[i]->typeId,
                         domain[domain[i]->ncsParent]->typeId );
    }

  }

  void Domains::SortDomains()  {
  //   This function sorts domains by their types, so that
  // AlignSimilarDomains should be called before.
  PDomain       D;
  ssm::PGraph   sg;
  PPDAlign      DA1,DA2;
  mmdb::ivector ix;
  int           i,j,k, i1,j1;
  bool          lt,direct;

    if (!domain)         return;
    if (nNCSParents<=1)  return;

    mmdb::GetVectorMemory ( ix,nNCSParents,0 );
    for (i=0;i<nNCSParents;i++)
      ix[i] = i;

    // Sort domains by increasing their type Id

    for (i=0;i<nNCSParents;i++)
      for (j=i+1;j<nNCSParents;j++)  {
        if (domain[j]->dclass==domain[i]->dclass)  {
          if (domain[j]->type==domain[i]->type)
               lt = (strcasecmp(domain[j]->range,domain[i]->range)<0);
          else lt = (domain[j]->type<domain[i]->type);
        } else
          lt = (domain[j]->dclass<domain[i]->dclass);
        if (lt)  {
          D = domain[i];
          domain[i] = domain[j];
          domain[j] = D;
          k     = ix[i];
          ix[i] = ix[j];
          ix[j] = k;
          if (G)  {
            sg   = G[i];
            G[i] = G[j];
            G[j] = sg;
          }
        }
      }

    // Now make new alignment matrix

    if (DA)  {
      DA1 = DA;
      DA  = NULL;
      AllocateAlignments();
      for (i=0;i<nNCSParents-1;i++)
        for (j=i+1;j<nNCSParents;j++)  {
          if (ix[i]<ix[j])  {
            i1 = ix[i];
            j1 = ix[j];
            direct = true;
          } else  {
            i1 = ix[j];
            j1 = ix[i];
            direct = false;
          }
          DA[i][j].Copy ( DA1[i1][j1],direct );
        }
      DA2 = DA;
      DA  = DA1;
      FreeAligns();  // erase old alignment matrix
      DA  = DA2;     // put new alignment matrix instead
    }

    // Set correct parential references in NCS-mates

    for (i=0;i<nNCSParents;i++)
      domain[i]->ncsParent = i;

    for (i=nNCSParents;i<nDomains;i++)  {
      k = -1;
      for (j=0;(j<nNCSParents) && (k<0);j++)
        if (ix[j]==domain[i]->ncsParent)  k = j;
      if (k>=0)  domain[i]->ncsParent = k;
    }

    mmdb::FreeVectorMemory ( ix,0 );

    // Reassign domain types so that they are numbered properly

    j = 0;
    k = 0;
    for (i=0;i<nNCSParents;i++)  {
      if (j!=domain[i]->type)  {
        j = domain[i]->type;
        k++;
      }
      domain[i]->type = k;
    }

    for (i=nNCSParents;i<nDomains;i++)
      domain[i]->type = domain[domain[i]->ncsParent]->type;

  }


  void Domains::randomizeReferenceFrames()  {
  //  This function makes pseudo-randomize orientations of symmetric
  // 1-atom domains.
  TFrame   F;
  mmdb::mat44    trm;
  mmdb::realtype dAngle,angle,dx,dy,dz;
  int      i,j,k,i1;

    k  = 0;
    i1 = -1;
    for (i=0;i<nNCSParents;i++)
      if (domain[i]->nAtoms==1)  {
        if (i1<0)  i1 = i;
        k++;
      }
    if (k>1)  {
      dAngle = 2.0*mmdb::Pi/(k+0.14159265);
      angle  = 0.0;
      for (i=0;i<nNCSParents;i++)
        if (domain[i]->nAtoms==1)  {
          if (i!=i1)  {
            mmdb::GetVecTMatrix ( trm,angle,1.0,1.0,1.0,
                                  domain[i1]->mx,domain[i1]->my,
                                  domain[i1]->mz );
            domain[i1]->getReferenceFrame ( F,trm );
            dx = domain[i]->mx - domain[i1]->mx;
            dy = domain[i]->my - domain[i1]->my;
            dz = domain[i]->mz - domain[i1]->mz;
            for (k=0;k<frameLen;k++)  {
              F[k][0] += dx;
              F[k][1] += dy;
              F[k][2] += dz;
            }
            domain[i]->setReferenceFrame ( F );
            for (j=0;j<3;j++)
              for (k=0;k<3;k++)
                domain[i]->itb[j][k] = trm[j][k];
          }
          angle += dAngle;
        }
    }

  }

  void Domains::correctReferenceFrames()  {
  //   This function corrects frame orientations using the results of
  //  structure alignment.
  mmdb::ovector b;
  mmdb::mat44   tm;
  int           i,j;

    mmdb::GetVectorMemory ( b,nDomains,0 );
    for (i=0;i<nDomains;i++)
      b[i] = false;

    for (i=0;i<nNCSParents;i++)
      if (!b[i])  {
        for (j=i+1;j<nDomains;j++)
          if (!b[j])  {
            if (isEquivalent(i,j))  {
              getTMatrix ( tm,i,j );
              domain[i]->getReferenceFrame ( domain[j]->frame,tm );
              b[j]= true;
            }
          }
        b[i] = true;
      }

    mmdb::FreeVectorMemory ( b,0 );

  }

  bool Domains::isAligned ( int domain1, int domain2 )  {
  int d1,d2;
    if (!DA)    return (nDomains==1) && (domain1==0) && (domain2==0);
    d1 = domain1;
    d2 = domain2;
    if (d1>=nNCSParents)  d1 = domain[d1]->ncsParent;
    if (d2>=nNCSParents)  d2 = domain[d2]->ncsParent;
    if (d1==d2) return true;
    if (d1<d2)  return (DA[d1][d2].Qscore>=0.0);
    return (DA[d2][d1].Qscore>=0.0);
  }

  bool Domains::isEquivalent ( int domain1, int domain2 )  {
  int d1,d2;
    if (!DA)    return (nDomains==1) && (domain1==0) && (domain2==0);
    d1 = domain1;
    d2 = domain2;
    if (d1>=nNCSParents)  d1 = domain[d1]->ncsParent;
    if (d2>=nNCSParents)  d2 = domain[d2]->ncsParent;
    if (d1==d2) return true;
    if (d1<d2)  return DA[d1][d2].equal;
    return DA[d2][d1].equal;
  }


  mmdb::realtype Domains::getQscore ( int domain1, int domain2 )  {
  int d1,d2;
    if (!DA)  {
      if ((nDomains==1) && (domain1==0) && (domain2==0))
           return  1.0;
      else return -1.0;
    }
    d1 = domain1;
    d2 = domain2;
    if (d1>=nNCSParents)  d1 = domain[d1]->ncsParent;
    if (d2>=nNCSParents)  d2 = domain[d2]->ncsParent;
    if (d1==d2) return 1.0;
    if (d1<d2)  return DA[d1][d2].Qscore;
    return DA[d2][d1].Qscore;
  }

  mmdb::realtype Domains::getRMSD ( int domain1, int domain2 )  {
  int d1,d2;
    if (!DA)  {
      if ((nDomains==1) && (domain1==0) && (domain2==0))
           return  0.0;
      else return -1.0;
    }
    d1 = domain1;
    d2 = domain2;
    if (d1>=nNCSParents)  d1 = domain[d1]->ncsParent;
    if (d2>=nNCSParents)  d2 = domain[d2]->ncsParent;
    if (d1==d2) return 0.0;
    if (d1<d2)  return DA[d1][d2].rmsd;
    return DA[d2][d1].rmsd;
  }

  mmdb::realtype Domains::getSeqId ( int domain1, int domain2 )  {
  int d1,d2;
    if (!DA)  {
      if ((nDomains==1) && (domain1==0) && (domain2==0))
           return  1.0;
      else return -1.0;
    }
    d1 = domain1;
    d2 = domain2;
    if (d1>=nNCSParents)  d1 = domain[d1]->ncsParent;
    if (d2>=nNCSParents)  d2 = domain[d2]->ncsParent;
    if (d1==d2) return 1.0;
    if (d1<d2)  return DA[d1][d2].seqId;
    return DA[d2][d1].seqId;
  }

  int Domains::getNalign ( int domain1, int domain2 )  {
  int d1,d2;
    if (!DA)  {
      if ((nDomains==1) && (domain1==0) && (domain2==0))
          return  domain[domain1]->nRes0;
     else return -1;
    }
    d1 = domain1;
    d2 = domain2;
    if (d1>=nNCSParents)  d1 = domain[d1]->ncsParent;
    if (d2>=nNCSParents)  d2 = domain[d2]->ncsParent;
    if (d1==d2) return domain[d1]->nRes0;
    if (d1<d2)  return DA[d1][d2].Nalign;
    return DA[d2][d1].Nalign;
  }

  void Domains::getTMatrix ( mmdb::mat44 & TMatrix, int domain1,
                                int domain2 )  {
  //  Returns matrix TMatrix such that TMatrix*domain1 superposes
  // on domain2
  mmdb::mat44 t1,t2;
  int   d1,d2;
    if ((!DA) || (domain1==domain2))
      mmdb::Mat4Init ( TMatrix );
    else  {
      d1 = domain1;
      d2 = domain2;
      if (d1>=nNCSParents)  d1 = domain[d1]->ncsParent;
      if (d2>=nNCSParents)  d2 = domain[d2]->ncsParent;
      if (d1<d2)  {
        mmdb::Mat4Inverse ( domain[domain1]->ncs_m,t1 );
        mmdb::Mat4Mult    ( t2,DA[d1][d2].TMatrix,t1 );
        mmdb::Mat4Mult    ( TMatrix,domain[domain2]->ncs_m,t2 );
      } else if (d1>d2)  {
        mmdb::Mat4Mult    ( t1,domain[domain1]->ncs_m,DA[d2][d1].TMatrix );
        mmdb::Mat4Inverse ( t1,t2 );
        mmdb::Mat4Mult    ( TMatrix,domain[domain2]->ncs_m,t2 );
      } else  {
        mmdb::Mat4Inverse ( domain[domain1]->ncs_m,t1 );
        mmdb::Mat4Mult    ( TMatrix,domain[domain2]->ncs_m,t1 );
      }
    }
  }

  mmdb::pstr Domains::getDomainRange ( mmdb::pstr S, int domainNo, int fieldLen ) {
    if (domain)
      return domain[domainNo]->getDomainRange ( S,fieldLen );
    strcpy ( S,"<no domains>" );
    return S;
  }

  int  Domains::getDomainClass ( int domainNo )  {
    if (domain)
      return domain[domainNo]->dclass;
    return DCLASS_None;
  }

  mmdb::pstr Domains::getChainID ( mmdb::pstr S, int domainNo )  {
    if (domain)
      return domain[domainNo]->getChainID ( S );
    strcpy ( S,"<no domains>" );
    return S;
  }

  int  Domains::getNAtoms ( int domainNo )  {
    if (domain)  return domain[domainNo]->nAtoms;
           else  return 0;
  }


  void Domains::getLigandList ( mmdb::PResName & ligName,
                                int            & nLigNames )  {
  char S[300];
  int  i,j;

    if (ligName)  {
      delete[] ligName;
      ligName = NULL;
    }
    nLigNames = 0;

    if (domain)
      for (i=0;i<nDomains;i++)
        if (domain[i]->dclass==DCLASS_Ligand)  {
          if (!ligName)
            ligName = new mmdb::ResName[nDomains];
          domain[i]->getLigandRange ( S );
          j = 0;
          while (j<nLigNames)
            if (!strcmp(S,ligName[j]))  break;
                                  else  j++;
          if (j>=nLigNames)
            strcpy ( ligName[nLigNames++],
                     domain[i]->getLigandRange(S) );

        }

    if (ligName)
      for (i=nLigNames;i<nDomains;i++)
        ligName[i][0] = char(0);

  }

  int Domains::getMonTypeOccurence ( int monType )  {
  // Returns the occurence number of specified monomer type in the
  // crystal.
  int i,n;

    n = 0;
    for (i=0;i<nDomains;i++)
      if (domain[i]->type==monType)
        n++;

    return n;

  }

  bool Domains::getAssembleStatus ( mmdb::pstr ligandName )  {
  char  S[300];
  int   i;
  bool  assemble;
    assemble = false;
    if (domain)
      for (i=0;i<nDomains;i++)  {
        domain[i]->getLigandRange ( S );
        if (!strcmp(S,ligandName))  {
          assemble = domain[i]->assemble;
          break;
        }
      }
    return assemble;
  }

  bool Domains::setAssembleStatus ( mmdb::pstr ligandName,
                                    bool       assemble )  {
  //  Returns true if status has changed
  char    S[300];
  int     i;
  bool changed;
    changed = false;
    if (domain)
      for (i=0;i<nDomains;i++)  {
        domain[i]->getLigandRange ( S );
        if ((!strcmp(S,ligandName)) &&
            ((assemble && (!domain[i]->assemble)) ||
             ((!assemble) && domain[i]->assemble)))  {
            changed = true;
            domain[i]->assemble = assemble;
        }
      }
    return changed;
  }

  bool Domains::assembleAll()  {
  int     i;
  bool assemble;
    assemble = true;
    if (domain)
      for (i=0;(i<nDomains) && assemble;i++)
        assemble = domain[i]->assemble;
    return assemble;
  }

  void Domains::excludeLigands ( mmdb::pstr ligList )  {
  // ligList must contain comma-separated list of ligand names
  // to be excluded from oligomeric state analysis. The list
  // must start and finish with commas, and contain no spaces.
  //   Example:  ",SO4,Zn,"
  //   Note:     ",!SO4,!Zn," will force inclusion
  char  S[300];
  int   i,k;
  bool  exclAgents;

    if (domain)  {

      if (!ligList)  {

        for (i=0;i<nDomains;i++)
          if (domain[i]->dclass==DCLASS_Ligand)
            domain[i]->assemble = true;

      } else  {

        if (strstr(ligList,lig_excl_all_key))  {

          for (i=0;i<nDomains;i++)
            if (domain[i]->dclass==DCLASS_Ligand)
              domain[i]->assemble = false;

          exclAgents = false;

        } else
          exclAgents = (strstr(ligList,lig_excl_agents_key)!=NULL);

        for (i=0;i<nDomains;i++)
          if (domain[i]->dclass==DCLASS_Ligand)  {
            if (domain[i]->assemble)  {
              if (exclAgents)
                domain[i]->assemble = !domain[i]->agent;
              if (domain[i]->assemble)  {
                domain[i]->getLigandRange ( S );
                k = strlen(S) + 2;
                S[k--] = char(0);
                S[k--] = ',';
                while (k>0)  {
                  S[k] = S[k-1];
                  k--;
                }
                S[0] = ',';
                domain[i]->assemble = (strstr(ligList,S)==NULL);
              }
            }
            if (!domain[i]->assemble)  {
              domain[i]->getLigandRange ( S );
              k = strlen(S) + 3;
              S[k--] = char(0);
              S[k--] = ',';
              while (k>1)  {
                S[k] = S[k-2];
                k--;
              }
              S[1] = '!';
              S[0] = ',';
              domain[i]->assemble = (strstr(ligList,S)!=NULL);
            }

          }

      }

    }

  }


  void Domains::setLigandsAssemble ( bool On )  {
  int i;
    if (domain)
      for (i=0;i<nDomains;i++)
        if (domain[i]->dclass==DCLASS_Ligand)
          domain[i]->assemble = On;
  }

  void Domains::getExclLigandList ( mmdb::PResName & ligName,
                                    int            & nLigNames )  {
  char S[300];
  int  i,j;

    if (ligName)  {
      delete[] ligName;
      ligName = NULL;
    }
    nLigNames = 0;

    if (domain)
      for (i=0;i<nDomains;i++)
        if ((domain[i]->dclass==DCLASS_Ligand) &&
            (!domain[i]->assemble))  {
          if (!ligName)
            ligName = new mmdb::ResName[nDomains];
          domain[i]->getLigandRange ( S );
          j = 0;
          while (j<nLigNames)
            if (!strcmp(S,ligName[j]))  break;
                                  else  j++;
          if (j>=nLigNames)
            strcpy ( ligName[nLigNames++],
                     domain[i]->getLigandRange(S) );

        }

    if (ligName)
      for (i=nLigNames;i<nDomains;i++)
        ligName[i][0] = char(0);

  }


  mmdb::pstr Domains::getDomainRange_html ( mmdb::pstr S, int domainNo )  {
    if (domain)
      return domain[domainNo]->getDomainRange_html ( S );
    strcpy ( S,"<no domains>" );
    return S;
  }

  PDomain Domains::getDomain ( int domainNo )  {
    if (domain)  {
      if ((domainNo>=0) && (domainNo<nDomains))
        return domain[domainNo];
    }
    return NULL;
  }

  int Domains::getNCSParent ( int domainNo )  {
    if (domain)  {
      if ((domainNo>=0) && (domainNo<nDomains))
        return domain[domainNo]->ncsParent;
    }
    return 0;
  }


  Brick::Brick()  {
    n        = NULL;
    nObjects = 0;
    nAlloc   = 0;
  }

  Brick::~Brick()  {
    mmdb::FreeVectorMemory ( n,0 );
  }

  void Brick::AddObject ( int objectNo )  {
  mmdb::ivector n1;
  int           i;
    if (nObjects>=nAlloc)  {
      nAlloc = nObjects + 20;
      mmdb::GetVectorMemory ( n1,nAlloc,0 );
      for (i=0;i<nObjects;i++)
        n1[i] = n[i];
      mmdb::FreeVectorMemory ( n,0 );
      n = n1;
    }
    n[nObjects++] = objectNo;
  }


  void Domains::RemoveBricks()  {
  int ix,iy,iz;

    if (DB)  {
      for (ix=0;ix<nx_bricks;ix++)
        if (DB[ix])  {
          for (iy=0;iy<ny_bricks;iy++)
            if (DB[ix][iy])  {
              for (iz=0;iz<nz_bricks;iz++)
                if (DB[ix][iy][iz])
                  delete DB[ix][iy][iz];
              delete[] DB[ix][iy];
            }
          delete[] DB[ix];
        }
      delete[] DB;
      DB = NULL;
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

  void Domains::MakeBricks()  {
  mmdb::realtype margin;
  int      i,j, ix,iy,iz;

    RemoveBricks();

    Dmax     =  0.0;
    x1_brick =  mmdb::MaxReal;
    x2_brick = -mmdb::MaxReal;
    y1_brick =  mmdb::MaxReal;
    y2_brick = -mmdb::MaxReal;
    z1_brick =  mmdb::MaxReal;
    z2_brick = -mmdb::MaxReal;

    for (i=0;i<nDomains;i++)  {
      if (domain[i]->R>Dmax)       Dmax     = domain[i]->R;
      if (domain[i]->mx<x1_brick)  x1_brick = domain[i]->mx;
      if (domain[i]->mx>x2_brick)  x2_brick = domain[i]->mx;
      if (domain[i]->my<y1_brick)  y1_brick = domain[i]->my;
      if (domain[i]->my>y2_brick)  y2_brick = domain[i]->my;
      if (domain[i]->mz<z1_brick)  z1_brick = domain[i]->mz;
      if (domain[i]->mz>z2_brick)  z2_brick = domain[i]->mz;
    }

    Dmax      += 5.0;  // plus somewhat bigger than water radius
    brick_size = 2.0*Dmax;
    margin     = 2.1*brick_size;
    x1_brick  -= margin;
    x2_brick  += margin;
    y1_brick  -= margin;
    y2_brick  += margin;
    z1_brick  -= margin;
    z2_brick  += margin;

    nx_bricks = mmdb::ifloor ( (x2_brick-x1_brick)/brick_size ) + 1;
    ny_bricks = mmdb::ifloor ( (y2_brick-y1_brick)/brick_size ) + 1;
    nz_bricks = mmdb::ifloor ( (z2_brick-z1_brick)/brick_size ) + 1;

    for (i=0;i<nDomains;i++)  {
      ix = mmdb::ifloor ( (domain[i]->mx-x1_brick)/brick_size );
      iy = mmdb::ifloor ( (domain[i]->my-y1_brick)/brick_size );
      iz = mmdb::ifloor ( (domain[i]->mz-z1_brick)/brick_size );
      if (!DB)  {
        DB = new PPPBrick[nx_bricks];
        for (j=0;j<nx_bricks;j++)
          DB[j] = NULL;
      }
      if (!DB[ix])  {
        DB[ix] = new PPBrick[ny_bricks];
        for (j=0;j<ny_bricks;j++)
          DB[ix][j] = NULL;
      }
      if (!DB[ix][iy])  {
        DB[ix][iy] = new PBrick[nz_bricks];
        for (j=0;j<nz_bricks;j++)
          DB[ix][iy][j] = NULL;
      }
      if (!DB[ix][iy][iz])
        DB[ix][iy][iz] = new Brick();
      DB[ix][iy][iz]->AddObject ( i );
    }

  }

  void Domains::getBricks ( PPBrick brick,
                             PDomain D, mmdb::mat44 & TMatrix )  {
  mmdb::realtype x,y,z;
  int      ix,iy,iz, i,i1,i2, j,j1,j2, k,k1,k2, n;

    for (i=0;i<27;i++)
      brick[i] = NULL;

    if (!DB)  return;

    x = TMatrix[0][0]*D->mx + TMatrix[0][1]*D->my +
        TMatrix[0][2]*D->mz + TMatrix[0][3];
    ix = mmdb::ifloor ( (x-x1_brick)/brick_size );
    if ((ix<0) || (ix>=nx_bricks))  return;

    y = TMatrix[1][0]*D->mx + TMatrix[1][1]*D->my +
        TMatrix[1][2]*D->mz + TMatrix[1][3];
    iy = mmdb::ifloor ( (y-y1_brick)/brick_size );
    if ((iy<0) || (iy>=ny_bricks))  return;

    z = TMatrix[2][0]*D->mx + TMatrix[2][1]*D->my +
        TMatrix[2][2]*D->mz + TMatrix[2][3];
    iz = mmdb::ifloor ( (z-z1_brick)/brick_size );
    if ((iz<0) || (iz>=nz_bricks))  return;

    n  = 0;
    i1 = mmdb::IMax ( ix-1,0 );
    i2 = mmdb::IMin ( ix+2,nx_bricks );
    j1 = mmdb::IMax ( iy-1,0 );
    j2 = mmdb::IMin ( iy+2,ny_bricks );
    k1 = mmdb::IMax ( iz-1,0 );
    k2 = mmdb::IMin ( iz+2,nz_bricks );
    for (i=i1;i<i2;i++)
      if (DB[i])
        for (j=j1;j<j2;j++)
          if (DB[i][j])
            for (k=k1;k<k2;k++)
              if (DB[i][j][k])  {
                if (DB[i][j][k]->nObjects>0)
                  brick[n++] = DB[i][j][k];
              }

  }


  void Domains::write ( mmdb::io::RFile f )  {
  int        i,j;
  mmdb::byte k;
  mmdb::byte Version=1;

    f.WriteByte ( &Version     );

    f.WriteReal ( &Dmax        );
    f.WriteInt  ( &nDomains    );
    f.WriteInt  ( &nDTypes     );
    f.WriteInt  ( &nNCSOps     );
    f.WriteInt  ( &nNCSParents );
    if (nDomains>0)  {
      k = 0x00;
      if (domain)  k |= 0x01;
      if (DA)      k |= 0x02;
      if (G)       k |= 0x04;
      f.WriteByte ( &k );
      if (domain)
        for (i=0;i<nDomains;i++)
          StreamWrite ( f,domain[i] );
      if (DA)
        for (i=0;i<nNCSParents-1;i++)
          for (j=i+1;j<nNCSParents;j++)
            DA[i][j].write ( f );
      if (G)
        for (i=0;i<nNCSParents;i++)
          StreamWrite ( f,G[i] );
    }
  }

  void Domains::read ( mmdb::io::RFile f )  {
  int        i,j;
  mmdb::byte Version,k;

    FreeMemory();

    f.ReadByte ( &Version     );

    f.ReadReal ( &Dmax        );
    f.ReadInt  ( &nDomains    );
    f.ReadInt  ( &nDTypes     );
    f.ReadInt  ( &nNCSOps     );
    f.ReadInt  ( &nNCSParents );

    if (nDomains>0)  {
      f.ReadByte ( &k );
      if (k & 0x01)  {
        domain = new PDomain[nDomains];
        for (i=0;i<nDomains;i++)  {
          domain[i] = NULL;
          StreamRead ( f,domain[i] );
          if ((i>=nNCSParents) && domain[i])
            domain[i]->CopyParentData ( domain );
        }
      }
      if (k & 0x02)  {
        DA = new PDAlign[nNCSParents-1];
        for (i=0;i<nNCSParents-1;i++)  {
          DA[i]  = new DAlign[nNCSParents-i-1];
          DA[i] -= (i+1);
          for (j=i+1;j<nNCSParents;j++)
            DA[i][j].read ( f );
        }
      }
      if (k & 0x04)  {
        G = new ssm::PGraph[nNCSParents];
        for (i=0;i<nNCSParents;i++)  {
          G[i] = NULL;
          StreamRead ( f,G[i] );
        }
      }
    }

  }


  MakeStreamFunctions(Domains)

}  // namespace pisa
