// $Id: pisa_assembly.cpp $
// =================================================================
//
//    24.03.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_assembly <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Assembly
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2004-2014
//
// =================================================================
//

#include <string.h>

#include "pisa_assembly.h"
#include "pisa_defs.h"

namespace pisa  {


  // =========================  Assembly  ===========================


  Assembly::Assembly()  {
    M           = NULL;  // list of AsmUnits (their ids) in assembly
    serNo       = 0;     // assembly serial number
    asmSize     = 0;     // assembly size (number of AsmUnits)
    mmSize      = 0;     // macromolecular size (number of macromols)
    freeSize    = 0;     // number of free (not fixed) monomers
    dissEn      = 0.0;   // maximal enthalpy of dissociation
    entropy     = 0.0;   // entropy change at dissociation
    freeEn      = 0.0;   // maximal free energy of dissociation
    dissIntArea = 0.0;   // dissociation interface area
    dissEn0     = 0.0;   // ground-level enthalpy of dissociation
    entropy0    = 0.0;   // ground-level entropy change at dissociation
    freeEn0     = 0.0;   // ground-level free energy of dissociation
    asa         = 0.0;   // accessible surface area
    bsa         = 0.0;   // buried surface area
    seGain      = 0.0;   // solvation energy gain
    nUC         = 0;     // number of this assemblies in unit cells
    type        = 0;     // assembly type id in assembly set
    nDiss       = 0;     // number of dissociating parts
    symNumber   = 0;     // assembly symmetry number
    stable      = false; // true if assembly is stable
  }

  Assembly::~Assembly()  {
    FreeMemory();
  }

  void Assembly::FreeMemory()  {
  int i;
    if (M)  {
      for (i=0;i<asmSize;i++)
        if (M[i])  delete M[i];
      delete[] M;
      M = NULL;
    }
    asmSize  = 0;
    mmSize   = 0;
    freeSize = 0;
  }

  mmdb::pstr Assembly::getFormula ( mmdb::pstr & S, PDomains D,
                                    int maxLen, mmdb::cpstr leadStr,
                                    bool html )  {
  mmdb::ivector n;
  char          N[30];
  int           i,j,k,m,m0,m1,l0;

    mmdb::GetVectorMemory ( n,D->nDomains,0 );
    for (i=0;i<D->nDomains;i++)
      n[i] = 0;

    for (i=0;i<asmSize;i++)
      n[M[i]->id]++;

    mmdb::CreateCopy ( S,"" );

    if (leadStr)  l0 = strlen ( leadStr );
            else  l0 = 0;
    m  = l0;
    m0 = 0;
    for (i=0;i<D->nDomains;i++)
      if (n[i]>0)  {
        k = n[i];
        for (j=i+1;j<D->nDomains;j++)
          if (n[j]>0)  {
            if (D->isEquivalent(i,j))  {
              k += n[j];
              n[j] = 0;
            }
          }
        m1 = m;
        if (k>1)  {
          if (html)  sprintf ( N,"<sub>%i</sub>",k );
               else  sprintf ( N,"(%i)",k );
          m += strlen(N);
        } else
          N[0] = char(0);
        m++;
        if (maxLen && (m-m0>maxLen) && (i<D->nDomains-1))  {
          if (html)  mmdb::CreateConcat ( S,"<br>" );
               else  mmdb::CreateConcat ( S,"\n",leadStr );
          m0 = m1;
          m += l0+1;
        }
        mmdb::CreateConcat ( S,D->domain[i]->typeId,N );
      }

    mmdb::FreeVectorMemory ( n,0 );

    return S;

  }

  mmdb::pstr Assembly::getComposition ( mmdb::pstr & S,
                                        PDomains     D,
                                        int          maxLen,
                                        mmdb::cpstr  leadStr,
                                        bool         html )  {
  mmdb::ivector n;
  mmdb::pstr    p;
  char          R[500];
  char          N[30];
  char          C;
  int           i,l,m,m0,m1,l0;

    C = char(0);  // only to keep compiler happy

    mmdb::GetVectorMemory ( n,D->nDomains,0 );
    for (i=0;i<D->nDomains;i++)
      n[i] = 0;

    for (i=0;i<asmSize;i++)
      n[M[i]->id]++;

    // check for repeat ligand names
    for (i=0;i<D->nDomains;i++)
      if (n[i]>0)  {
        p = mmdb::FirstOccurence ( D->domain[i]->range,']' );
        if (p)  {
          C    = p[1];
          p[1] = char(0);
          l    = strlen ( D->domain[i]->range );
          p[1] = C;
          for (m=i+1;m<D->nDomains;m++)
            if ((n[m]>0) &&
                (!strncmp(D->domain[i]->range,D->domain[m]->range,l))) {
              n[i] += n[m];
              n[m]  = 0;
            }
        }
      }

    mmdb::CreateCopy ( S,"" );

    if (leadStr)  l0 = strlen ( leadStr );
            else  l0 = 0;
    m  = l0;
    m0 = 0;
    for (i=0;i<D->nDomains;i++)
      if (n[i]>0)  {
        p = mmdb::FirstOccurence ( D->domain[i]->range,']' );
        if (p)  p++;
          else  p = mmdb::FirstOccurence ( D->domain[i]->range,':' );
        if (p)  {
          C  = *p;
          *p = char(0);
        }
        m1 = m;
        D->domain[i]->getDomainRange ( R,-1 );
        m += strlen ( R );
        if (n[i]>1)  {
          if (html)  sprintf ( N,"<sub>%i</sub>",n[i] );
               else  sprintf ( N,"(%i)",n[i] );
          m += strlen(N);
        } else
          N[0] = char(0);
        if (maxLen && (m-m0>maxLen) && (i<D->nDomains-1))  {
          if (html)  mmdb::CreateConcat ( S,"<br>" );
               else  mmdb::CreateConcat ( S,"\n",leadStr );
          m0 = m1;
          m += l0+1;
        }
        mmdb::CreateConcat ( S,R,N );
        if (p)  *p = C;
      }

    mmdb::FreeVectorMemory ( n,0 );

    return S;

  }

  mmdb::pstr Assembly::getDissPattern ( mmdb::pstr & S,
                                        PDomains     D,
                                        int          maxLen,
                                        mmdb::cpstr  leadStr,
                                        bool         html )  {
  mmdb::ivector  n,dl;
  mmdb::psvector d;
  mmdb::pstr     p;
  char           R[500];
  char           N[30];
  char           C;
  int            i,j,k,m,nd,m0,m1,l0;
  bool           B;

    C  = char(0);  // only to keep compiler happy
    m1 = 0;        // only to keep compiler happy

    if (leadStr)  l0 = strlen ( leadStr );
            else  l0 = 0;

    mmdb::GetVectorMemory ( n,D->nDomains,0 );

    mmdb::GetVectorMemory ( d ,asmSize,1 );
    mmdb::GetVectorMemory ( dl,asmSize,1 );
    for (i=1;i<=asmSize;i++)  {
      d [i] = NULL;
      dl[i] = 0;
    }

    k = 0;
    do {

      k++;

      for (i=0;i<D->nDomains;i++)
        n[i] = 0;

      B = false;
      for (i=0;i<asmSize;i++)
        if (M[i]->dissAsm==k)  {
          n[M[i]->id]++;
          B = true;
        }

      if (B)  {

        // check for repeated ligand names
        for (i=0;i<D->nDomains;i++)
          if (n[i]>0)  {
            p = mmdb::FirstOccurence ( D->domain[i]->range,']' );
            if (p)  {
              C    = p[1];
              p[1] = char(0);
              m    = strlen ( D->domain[i]->range );
              p[1] = C;
              for (j=i+1;j<D->nDomains;j++)
                if ((n[j]>0) &&
                    (!strncmp(D->domain[i]->range,
                              D->domain[j]->range,m))) {
                  n[i] += n[j];
                  n[j]  = 0;
                }
            }
          }

        m  = l0;
        m0 = 0;
        mmdb::CreateCopy ( d[k],"" );
        for (i=0;i<D->nDomains;i++)
          if (n[i]>0)  {
            p = mmdb::FirstOccurence ( D->domain[i]->range,']' );
            if (p)  p++;
              else  p = mmdb::FirstOccurence ( D->domain[i]->range,':' );
            if (p)  {
              C  = *p;
              *p = char(0);
            }
            m1 = m;
            D->domain[i]->getDomainRange ( R,-1 );
            m += strlen ( R );
            if (n[i]>1)  {
              if (html)  sprintf ( N,"(%i)",n[i] );
                   else  sprintf ( N,"(%i)",n[i] );
              m += strlen(N);
            } else
              N[0] = char(0);
            if (maxLen && (m-m0>maxLen) && (i<D->nDomains-1))  {
              if (html)  mmdb::CreateConcat ( d[k],"<br>" );
                   else  mmdb::CreateConcat ( d[k],"\n",leadStr );
              m0 = m1;
              m += l0+1;
            }
            mmdb::CreateConcat ( d[k],R,N );
            dl[k] = strlen ( d[k] );
            if (p)  *p = C;
          }

      } else
        k--;
    } while (B);

    mmdb::CreateCopy ( S,"" );

    m  = l0;
    m0 = 0;
    for (i=1;i<=k;i++)
      if (d[i])  {

        nd = 1;
        for (j=i+1;j<=k;j++)
          if (d[j])  {
            if (!strcmp(d[i],d[j]))  {
              delete[] d[j];
              d[j] = NULL;
              nd++;
            }
          }

        if (i>1)  {
          m += 3;
          mmdb::CreateConcat ( S," + " );
          m1 = m;
        }

        if (nd>1)  {
          sprintf ( N,"%i{",nd );
          m += strlen(N)+1;
        } else
          N[0] = char(0);
        m += dl[i];

        if (maxLen && (m-m0>maxLen) && (i<D->nDomains-1))  {
          if (html)  mmdb::CreateConcat ( S,"<br>" );
               else  mmdb::CreateConcat ( S,"\n",leadStr );
          m0 = m1;
          m += l0+1;
        }

        mmdb::CreateConcat ( S,N,d[i] );
        if (nd>1)
          mmdb::CreateConcat ( S,"}" );

      }

    for (i=1;i<=asmSize;i++)
      if (d[i])  delete[] d[i];

    mmdb::FreeVectorMemory ( d ,1 );
    mmdb::FreeVectorMemory ( dl,1 );

    mmdb::FreeVectorMemory ( n,0 );

    return S;

  }

  mmdb::pstr Assembly::getFileSafeName ( mmdb::pstr & S, PDomains D )  {
  int i,j;
    getFormula ( S,D,0,"",false );
    i = 0;
    j = 0;
    while (S[i])
      if ((S[i]=='(') || (S[i]==')'))  i++;
      else if (j<i)  S[j++] = S[i++];
      else { j++; i++; }
    S[j] = char(0);
    return S;
  }


  void Assembly::getEngagedInterfaces ( mmdb::ivector intf,
                                         int     nInterfaces )  {
  //  Returns the list of engaged interfaces; intf[i] is the
  // occurence number of ith interface in the assembly.
  int i,j,k;

    for (i=0;i<nInterfaces;i++)  {
      intf[i] = 0;
      for (j=0;j<asmSize;j++)  {
        k = 0;
        while (k<_max_n_int)
          if (M[j]->intfl[k])  {
            if (M[j]->intfl[k][i]>=0)  intf[i]++;
            k++;
          } else
            k = _max_n_int;
      }
      intf[i] /= 2;
    }

  }

  void Assembly::getDissInterfaces ( mmdb::ivector intf,
                                      int     nInterfaces )  {
  //  Returns the list of dissociating interfaces; intf[i] is the
  // occurence number of ith interface in the assembly.
  int i,j,k;

    for (i=0;i<nInterfaces;i++)  {
      intf[i] = 0;
      for (j=0;j<asmSize;j++)  {
        k = 0;
        while (k<_max_n_int)
          if (M[j]->intfl[k])  {
            if (M[j]->intfl[k][i]>=0)  {
              if (M[M[j]->intfl[k][i]]->dissAsm!=M[j]->dissAsm)
                intf[i]++;
            }
            k++;
          } else
            k = _max_n_int;
      }
      intf[i] /= 2;
    }

  }

  mmdb::PManager Assembly::getAssemblyStructure ( mmdb::PManager MMDB,
                                                  PDomains       D )  {
  mmdb::PManager AsmStructure;
  mmdb::PModel   model;
  mmdb::PPChain  chain;
  mmdb::PChain   chn;
  mmdb::ivector  selHnd;
  int            i,k,nChains;

    mmdb::GetVectorMemory ( selHnd,D->nDomains,0 );
    for (i=0;i<D->nDomains;i++)
      selHnd[i] = 0;

    model = mmdb::newModel();
    for (i=0;i<asmSize;i++)  {
      k = M[i]->id;
      if (!selHnd[k])
        selHnd[k] = D->domain[k]->SelectDomain (
                                      MMDB,mmdb::STYPE_CHAIN,1,false );
      MMDB->GetSelIndex ( selHnd[k],chain,nChains );
      if (nChains>0) {
        chn = mmdb::newChain();
        chn->Copy ( chain[0] );
        chn->SetChainID ( M[i]->visualID );
        chn->ApplyTransform ( M[i]->T );
        model->AddChain ( chn );
      }
    }

    AsmStructure = new mmdb::Manager();
    AsmStructure->AddModel ( model );

    for (i=0;i<D->nDomains;i++)
      if (selHnd[i])  MMDB->DeleteSelection ( selHnd[i] );

    mmdb::FreeVectorMemory ( selHnd,0 );

    return AsmStructure;

  }

  /*
  mmdb::PManager Assembly::getAssemblyStructure ( mmdb::PManager MMDB,
                                                  PDomains   D )  {
  mmdb::PManager AsmStructure;
  mmdb::PModel       model;
  mmdb::PPChain      chain;
  mmdb::PChain       chn;
  mmdb::ivector       selHnd;
  ChainID       chID;
  int           i,k,nChains;

    mmdb::GetVectorMemory ( selHnd,D->nDomains,0 );
    for (i=0;i<D->nDomains;i++)
      selHnd[i] = 0;

    model = newCModel();
    strcpy ( chID,"A" );
    for (i=0;i<asmSize;i++)  {
      k = M[i]->id;
      if (!selHnd[k])
        selHnd[k] = D->domain[k]->SelectDomain (
                                            MMDB,mmdb::STYPE_CHAIN,1,false );
      MMDB->GetSelIndex ( selHnd[k],chain,nChains );
      if (nChains>0) {
        chn = newCChain();
        chn->Copy ( chain[0] );
        chn->SetChainID ( chID );
        chn->ApplyTransform ( M[i]->T );
        model->AddChain ( chn );
        chID[0]++;
        if (chID[0]>'Z')  chID[0] = 'A';
      }
    }

    AsmStructure = new CMMDBManager();
    AsmStructure->AddModel ( model );

    for (i=0;i<D->nDomains;i++)
      if (selHnd[i])  MMDB->DeleteSelection ( selHnd[i] );

    mmdb::FreeVectorMemory ( selHnd,0 );

    return AsmStructure;

  }
  */


  int Assembly::getFirstMonType()  {
  // Returns type of first macromolecular monomeric unit, or first
  // ligand unit if mmSize==0.
  int monType;
  int i;

    if (asmSize<=0)  return -1;
    if (mmSize<=0)   return M[0]->type;

    monType = -2;
    for (i=0;(i<asmSize) && (monType<0);i++)
      if (M[i]->dclass!=DCLASS_Ligand)
        monType = M[0]->type;

    return monType;

  }

  int Assembly::getMonTypeOccurence ( int monType )  {
  // Returns the occurence number of specified monomer type in the
  // assembly.
  int i,n;

    n = 0;
    for (i=0;i<asmSize;i++)
      if ((M[i]->type==monType) && (!M[i]->fixed))
        n++;

    return n;

  }


  void Assembly::calcStats ( int         domainID,
                             int &       nOcc, // occurence number
                             mmdb::realtype & dom_bsa,
                             mmdb::realtype & dom_DeltaG,
                             int &       dom_nHB,
                             int &       dom_nSB,
                             int &       dom_nDS,
                             int &       dom_nAt,  // no. of buried atoms
                             int &       dom_nRes, // no. of buried res-s
                             PPInterface interface,
                             PPDomain    D,
                             int         nInterfaces )  {
  int i,j,k,n, dt,dt1,dt2;

    nOcc = 0;

    dom_bsa    = 0.0;
    dom_DeltaG = 0.0;
    dom_nHB    = 0;
    dom_nSB    = 0;
    dom_nDS    = 0;
    dom_nAt    = 0;
    dom_nRes   = 0;

    dt = D[domainID]->ncsParent;
    for (i=0;i<asmSize;i++)
      if (M[i]->id==domainID)  {
        nOcc++;
        for (j=0;j<nInterfaces;j++)  {
          k = 0;
          n = 0;
          while (n<_max_n_int)
            if (M[i]->intfl[n])  {
              if (M[i]->intfl[n][j]>=0)  k++;
              n++;
            } else
              n = _max_n_int;
          if (k)  {
            dt1 = D[interface[j]->domain1]->ncsParent;
            dt2 = D[interface[j]->domain2]->ncsParent;
            if (dt1==dt)  {
              dom_bsa += k*interface[j]->intArea;
              dom_nHB += k*interface[j]->nHBonds;
              dom_nSB += k*interface[j]->nSBridges;
              if (dt2==dt)  {
                dom_DeltaG += k*(interface[j]->intDeltaG1+
                                 interface[j]->intDeltaG2)/2.0;
                dom_nAt    += k*(interface[j]->nIntAtoms1+
                                 interface[j]->nIntAtoms2)/2;
                dom_nRes   += k*(interface[j]->nIntRes1+
                                 interface[j]->nIntRes2)/2;
              } else  {
                dom_DeltaG += k*interface[j]->intDeltaG1;
                dom_nAt    += k*interface[j]->nIntAtoms1;
                dom_nRes   += k*interface[j]->nIntRes1;
              }
            } else if (dt2==dt)  {
              dom_bsa    += k*interface[j]->intArea;
              dom_DeltaG += k*interface[j]->intDeltaG2;
              dom_nHB    += k*interface[j]->nHBonds;
              dom_nSB    += k*interface[j]->nSBridges;
              dom_nDS    += k*interface[j]->nDSBonds;
              dom_nAt    += k*interface[j]->nIntAtoms2;
              dom_nRes   += k*interface[j]->nIntRes2;
            }
          }
        }
      }

    if (nOcc>0)  {
      dom_bsa    /= nOcc;
      dom_DeltaG /= nOcc;
      dom_nHB    /= nOcc;
      dom_nSB    /= nOcc;
      dom_nDS    /= nOcc;
      dom_nAt    /= nOcc;
      dom_nRes   /= nOcc;
    }

  }

  void Assembly::makeChainMapping()  {
  //  makeChainMapping() fills AsmUnit::visualID fields, which
  //  are used as chain IDs when assembly is visualised or
  //  downloaded
  mmdb::ChainID chID;
  int           i;
    strcpy ( chID,"A" );
    for (i=0;i<asmSize;i++)  {
      strcpy ( M[i]->visualID,chID );
      chID[0]++;
      if (chID[0]>'Z')  chID[0] = 'A';
    }
  }

  bool isInteger ( mmdb::realtype x )  {
    return (fabs(mmdb::mround(x)-x) < 1.0e-8);
  }

  #define  matEps  0.00001

  void Assembly::makeOrientation ( PDomains   D,
                                   mmdb::PManager MMDB )  {
  //  Assemblies that are coming from the assembling engine
  // (CAssembler), may be placed in any unit cell and in
  // any part of a unit cell. This function brings the assembly
  // into original ASU and unit cell, if that is necessary and
  // possible.
  mmdb::SymOp    SymOp;
  mmdb::pmat44   T,symT;
  mmdb::ivector  occ;
  mmdb::mat44    rom,rfm,t1,t2;
  char           symOpTitle[300];
  int            i,j,k,n,i0,nSymOps;

    MMDB->GetROMatrix ( rom );

    mmdb::Mat4Inverse ( rom,rfm );

    //  1. Make fractional matrices for the monomers
    T = new mmdb::mat44[asmSize];
    mmdb::GetVectorMemory ( occ,asmSize,0 );
    for (i=0;i<asmSize;i++)  {
      k = M[i]->id;
      if (k>=D->nNCSParents)  {
        occ[i] = -3;  // ignore NCS mates when choosing base matrix
        mmdb::Mat4Div2 ( t2,M[i]->T,D->domain[k]->ncs_m );
        mmdb::Mat4Mult ( t1,t2,rom );
      } else  {
        occ[i] = -2;
        mmdb::Mat4Mult ( t1,M[i]->T,rom );
      }
      mmdb::Mat4Mult ( T[i],rfm,t1 );
      M[i]->ncsOpNo = D->domain[k]->ncsOpNo;
    }

    //  2. Find the most frequently met matrix, which does not
    //     include an NCS transformation
    i0 = 0;
    n  = 0;
    for (i=0;i<asmSize;i++)
      if (occ[i]==-2)  {
        k = 1;
        for (j=0;j<asmSize;j++)
          if ((j!=i) && (occ[j]==-2) &&
              mmdb::isMat4Eq(T[i],T[j],matEps,false))  {
            occ[j] = -1;
            k++;
          }
        occ[i] = -1;
        if (k>n)  {
          n  = k;
          i0 = i;
        }
      }

    //  3. Bring all monomers into orientation T[i0]^{-1}, so that
    //     the most frequent orientation is now in original ASU of
    //     the original unit cell

    mmdb::Mat4Inverse ( T[i0],t1 );

    nSymOps = MMDB->GetNumberOfSymOps();
    symT    = new mmdb::mat44[nSymOps];

    for (i=0;i<nSymOps;i++)
      MMDB->GetSymOpMatrix ( symT[i],i );

    for (i=0;i<asmSize;i++)  {

      mmdb::Mat4Copy ( T[i],t2    );  // t2 = T[i]
      mmdb::Mat4Mult ( T[i],t1,t2 );  // T[i] = t1*T[i]

      k = -1;
      for (j=0;(j<nSymOps) && (k<0);j++)
        if (mmdb::isMat4Eq(T[i],symT[j],matEps,true))  {
          if (isInteger(T[i][0][3]-symT[j][0][3]) &&
              isInteger(T[i][1][3]-symT[j][1][3]) &&
              isInteger(T[i][2][3]-symT[j][2][3]))
            k = j;
        }
      if (k>=0)  {
        M[i]->symOpNo    = k;
        M[i]->rcsb_symop = k;  // just for the moment, mapped later
        M[i]->cell_i     = mmdb::mround ( T[i][0][3]-symT[k][0][3] );
        M[i]->cell_j     = mmdb::mround ( T[i][1][3]-symT[k][1][3] );
        M[i]->cell_k     = mmdb::mround ( T[i][2][3]-symT[k][2][3] );
      }

      SymOp.CompileOpTitle ( symOpTitle,T[i],false  );
      mmdb::CreateCopy     ( M[i]->symOp,symOpTitle );

    }

    //  4. Transform all matrices back to the orthogonal space
    //     and write them into monomer structures
    for (i=0;i<asmSize;i++)  {
      k = M[i]->id;
      mmdb::Mat4Mult ( t1,T[i],rfm );
      if (k>=D->nNCSParents)  {
        mmdb::Mat4Mult ( t2,rom,t1 );
        mmdb::Mat4Mult ( M[i]->T,t2,D->domain[k]->ncs_m );
      } else
        mmdb::Mat4Mult ( M[i]->T,rom,t1 );
    }

    //  5. Dispose of temporary arrays

    if (symT)  delete[] symT;

    mmdb::FreeVectorMemory ( occ,0 );
    if (T)  delete[] T;

  }

  void Assembly::countOriginalOrientations ( PPDomain D,
                                             mmdb::ivector   icnt )  {
  int i,k;
    for (i=0;i<asmSize;i++)  {
      k = M[i]->id;
      if ((M[i]->symOpNo==0) && (M[i]->ncsOpNo==0) &&
          (D[k]->dclass!=DCLASS_Ligand))
        icnt[k]++;
    }
  }

  void Assembly::Orth2Frac ( mmdb::PManager MMDB )  {
  int i;
    for (i=0;i<asmSize;i++)
      MMDB->Orth2Frac ( M[i]->T,M[i]->TF );
  }

  void Assembly::assignRCSBSymOps ( PRCSBData rcsbData )  {
  int i;
    if (rcsbData->nSymOps>0)  {
      for (i=0;i<asmSize;i++)
        M[i]->rcsb_symop = rcsbData->rcsb_symop[M[i]->symOpNo];
    }
  }

  void Assembly::Copy ( PMultimer U, PPInterface Interface,
                        int nInterfaces, int nInterfaces0,
                        PPDomain D, mmdb::mat44 & rom )  {
  int     i;
  bool fixed;

    FreeMemory();
    if (!U)  return;

    asmSize  = U->mSize;
    mmSize   = 0;
    freeSize = 0;
    nDiss    = 0;

    if (asmSize>0)  {

      for (i=0;i<asmSize;i++)  {
        U->R[i].M->id2 = i+1;
        U->R[i].M->U   = U;
      }

      M = new PAsmUnit[asmSize];

      for (i=0;i<asmSize;i++)  {

        fixed = false;
        if (U->R[i].M->dclass!=DCLASS_Ligand)  {
          mmSize++;
          freeSize++;
        } else if (!isFixed(U->R[i].M->id,Interface,nInterfaces))  {
          fixed = true;
          freeSize++;
        }

        M[i] = new AsmUnit();
        M[i]->makeUnit ( U->R[i],Interface,nInterfaces,nInterfaces0,D,rom );
        M[i]->fixed = fixed;

      }

      if (U->dissMult)  {
        for (i=0;i<asmSize;i++)  {
          M[i]->dissAsm = U->dissMult[i];
          if (M[i]->dissAsm>nDiss)  nDiss = M[i]->dissAsm;
        }
      }

      for (i=0;i<asmSize;i++)
        U->R[i].M->id2 = 0;

    }

    dissEn      = U->dissEn;       // maximal enthalpy of dissociation
    entropy     = U->entropy;      // entropy change at dissociation
    freeEn      = -dissEn-entropy; // ground-level free energy of diss-n
    dissIntArea = U->dissIntArea;  // dissociation interface area
    dissEn0     = U->dissEn0;      // ground-level enthalpy of diss-n
    entropy0    = U->entropy0;     // ground-level entropy change at diss-n
    freeEn0     = -dissEn0-entropy0; // ground-level free energy of diss-n
    asa         = U->asa;          // accessible surface area
    bsa         = U->bsa;          // buried surface area
    seGain      = U->seGain;       // solvation energy gain
    nUC         = U->nUC;          // number of this assemblies in u.cells
    type        = U->type;         // assembly type id in assembly set
    stable      = U->stable;       // true if assembly is stable

  }

  void Assembly::averageData ( PAssembly A, int & nSample )  {
  int i;
  int n1 = nSample + 1;

    if (nSample<=0)  {

      FreeMemory();

      serNo     = A->serNo;     // assembly serial number
      asmSize   = A->asmSize;   // assembly size (number of AsmUnits)
      mmSize    = A->mmSize;    // macromolecular size (number of macromols)
      freeSize  = A->freeSize;  // number of free (not fixed) monomers

      nUC       = A->nUC;       // number of this assemblies in unit cells
      type      = A->type;      // assembly type id in assembly set
      nDiss     = A->nDiss;     // number of dissociating parts
      symNumber = A->symNumber; // assembly symmetry number
      stable    = A->stable;    // true if assembly is stable

      M = new PAsmUnit[asmSize];
      for (i=0;i<asmSize;i++)  {
        M[i] = new AsmUnit();
        M[i]->Copy ( A->M[i] );
      }

  /*
    dissEn      = A->dissEn;
    entropy     = A->entropy;
    freeEn      = A->freeEn;
    dissIntArea = A->dissIntArea;
    dissEn0     = A->dissEn0;
    entropy0    = A->entropy0;
    freeEn0     = A->freeEn0;
    asa         = A->asa;
    bsa         = A->bsa;
    seGain      = A->seGain;
  */

    }

    dissEn      = (nSample*dissEn      + A->dissEn     ) / n1;
    entropy     = (nSample*entropy     + A->entropy    ) / n1;
    freeEn      = (nSample*freeEn      + A->freeEn     ) / n1;
    dissIntArea = (nSample*dissIntArea + A->dissIntArea) / n1;
    dissEn0     = (nSample*dissEn0     + A->dissEn0    ) / n1;
    entropy0    = (nSample*entropy0    + A->entropy0   ) / n1;
    freeEn0     = (nSample*freeEn0     + A->freeEn0    ) / n1;
    asa         = (nSample*asa         + A->asa        ) / n1;
    bsa         = (nSample*bsa         + A->bsa        ) / n1;
    seGain      = (nSample*seGain      + A->seGain     ) / n1;

    nSample = n1;

  }


  mmdb::xml::PXMLObject Assembly::getAssemblyXML ( PDomains    D,
                                                   PInterfaces PI,
                                                   int   nCellOut,
                                                   int      Score )  {
  mmdb::xml::PXMLObject xml;
  mmdb::xml::PXMLObject xmlInterfaces;
  mmdb::xml::PXMLObject xmlInterface;
  PInterface            interface;
  mmdb::ivector         intf,dintf;
  mmdb::pstr            F;
  int                   i,nInterfaces,n;

    F     = NULL;
    intf  = NULL;
    dintf = NULL;

    xml = new mmdb::xml::XMLObject ( xml_assembly );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_ser_no  ,serNo    ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_id      ,type+1   ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_size    ,asmSize  ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_mmsize  ,mmSize   ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_freesize,freeSize ) );

    if (Score>=0)  {
      if (asmSize<=1)
        xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_score,
          "Analysis of crystal interfaces has not revealed any "
          "strong indications that this assembly may form stable "
          "complexes in solution." ) );
      else if (Score<=2)
        xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_score,
          "This assembly appears to be stable in solution." ) );
      else if (Score<=5)
        xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_score,
          "This assembly falls into a grey region of "
          "complexation criteria. It may or may not be "
          "stable in solution." ) );
      else
        xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_score,
          "This assembly may be formed from crystallographic "
          "considerations, however it does not appear to be "
          "stable in solution." ) );
    }

    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_diss_energy,freeEn    ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_diss_energy_0,freeEn0 ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_asa        ,asa       ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_bsa        ,bsa       ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_entropy    ,entropy   ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_entropy_0  ,entropy0  ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_diss_area,dissIntArea ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_int_dg     ,seGain    ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_nuc        ,nUC       ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_ndiss      ,nDiss     ) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_sym_num    ,symNumber ) );


    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_formula,
                                          getFormula(F,D,0,"",false)) );

    xml->AddObject ( new mmdb::xml::XMLObject(xml_asm_composition,
                                      getComposition(F,D,0,"",false)) );


    nInterfaces = PI->getNofInterfaces();

    mmdb::GetVectorMemory ( intf,nInterfaces,0  );
    getEngagedInterfaces  ( intf,nInterfaces    );

    mmdb::GetVectorMemory ( dintf,nInterfaces,0 );
    getDissInterfaces     ( dintf,nInterfaces   );

    xmlInterfaces = new mmdb::xml::XMLObject ( xml_asm_interfaces );

    n = 0;
    for (i=0;i<nInterfaces;i++)
      if (intf[i]>0)  {
        interface = PI->getInterface ( i );
        if (interface)  {
          n++;
          xmlInterface = new mmdb::xml::XMLObject ( xml_asm_interface );
          xmlInterface->AddObject ( new mmdb::xml::XMLObject (
                                    xml_asm_intf_id,interface->id ) );
          xmlInterface->AddObject ( new mmdb::xml::XMLObject (
                                    xml_asm_intf_diss,(dintf[i]>0) ) );
          xmlInterfaces->AddObject ( xmlInterface );
        }
      }

    xmlInterfaces->InsertObject (
                         new mmdb::xml::XMLObject(xml_asm_nints,n),0 );
    xml->AddObject ( xmlInterfaces );


    for (i=0;i<asmSize;i++)
      xml->AddObject ( M[i]->getAsmUnitXML(D->domain,nCellOut) );

    mmdb::FreeVectorMemory ( intf ,0 );
    mmdb::FreeVectorMemory ( dintf,0 );
    if (F)  delete[] F;

    return xml;

  }

  mmdb::cpstr stateNames[20] = {
     "MONOMERIC",
     "DIMERIC",
     "TRIMERIC",
     "TETRAMERIC",
     "PENTAMERIC",
     "HEXAMERIC",
     "HEPTAMERIC",
     "OCTAMERIC",
     "NONAMERIC",
     "DECAMERIC",
     "UNDECAMERIC",
     "DODECAMERIC",
     "TRIDECAMERIC",
     "TETRADECAMERIC",
     "PENTADECAMERIC",
     "HEXADECAMERIC",
     "HEPTADECAMERIC",
     "OCTADECAMERIC",
     "NONADECAMERIC",
     "EICOSAMERIC"
  };

  mmdb::pstr Assembly::sizeName ( mmdb::pstr S )  {
  // returns S filled
    if (mmSize<=20)
          strcpy  ( S,stateNames[mmSize-1] );
    else  sprintf ( S,"%iMERIC",mmSize     );
    return S;
  }

  void Assembly::write ( mmdb::io::RFile f, int nInterfaces )  {
  int  i;
  mmdb::byte Version=3;

    f.WriteByte ( &Version );

    if (!M)  asmSize = -asmSize;

    f.WriteInt ( &serNo    );
    f.WriteInt ( &asmSize  );
    f.WriteInt ( &mmSize   );
    f.WriteInt ( &freeSize );

    for (i=0;i<asmSize;i++)
      M[i]->write ( f,nInterfaces );

    if (!M)  asmSize = -asmSize;

    f.WriteReal ( &dissEn      );
    f.WriteReal ( &entropy     );
    f.WriteReal ( &freeEn      );
    f.WriteReal ( &dissIntArea );
    f.WriteReal ( &dissEn0     );
    f.WriteReal ( &entropy0    );
    f.WriteReal ( &freeEn0     );
    f.WriteReal ( &asa         );
    f.WriteReal ( &bsa         );
    f.WriteReal ( &seGain      );
    f.WriteInt  ( &nUC         );
    f.WriteInt  ( &type        );
    f.WriteInt  ( &nDiss       );
    f.WriteInt  ( &symNumber   );
    f.WriteBool ( &stable      );

  }

  void Assembly::read  ( mmdb::io::RFile f, int nInterfaces )  {
  int  i;
  mmdb::byte Version;

    FreeMemory();

    f.ReadByte ( &Version );

    if (Version>1)
      f.ReadInt ( &serNo );
    f.ReadInt ( &asmSize );
    f.ReadInt ( &mmSize  );
    if (Version>2)
          f.ReadInt ( &freeSize  );
    else  freeSize = asmSize;  // just to keep it positive

    if (asmSize>0)  {
      M = new PAsmUnit[asmSize];
      for (i=0;i<asmSize;i++)  {
        M[i] = new AsmUnit();
        M[i]->read ( f,nInterfaces );
      }
    } else
      asmSize = -asmSize;

    f.ReadReal ( &dissEn      );
    f.ReadReal ( &entropy     );
    if (Version>2)
      f.ReadReal ( &freeEn    );
    f.ReadReal ( &dissIntArea );
    if (Version>2)  {
      f.ReadReal ( &dissEn0   );
      f.ReadReal ( &entropy0  );
      f.ReadReal ( &freeEn0   );
    } else  {
      dissEn0  = 0.0;
      entropy0 = 0.0;
      freeEn   = -dissEn  - entropy;
      freeEn0  = -dissEn0 - entropy0;
    }
    f.ReadReal ( &asa         );
    f.ReadReal ( &bsa         );
    f.ReadReal ( &seGain      );
    f.ReadInt  ( &nUC         );
    f.ReadInt  ( &type        );
    f.ReadInt  ( &nDiss       );
    f.ReadInt  ( &symNumber   );
    f.ReadBool ( &stable      );

  }


}  // namespace pisa
