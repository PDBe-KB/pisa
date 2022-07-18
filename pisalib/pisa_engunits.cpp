// $Id: pisa_asmunits.cpp $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  AsmUnits <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Monomer
//       ~~~~~~~~~  pisa::Multimer
//                  pisa::IntfRecStack
//                  pisa::MultimerSet
//
//  (C) E. Krissinel 2007-2014
//
// =================================================================
//

#include <math.h>

#include "pisa_engunits.h"
#include "pisa_types.h"

#include "mmdb2/mmdb_math_linalg.h"

namespace pisa  {

  //  ================================================================

  static mmdb::realtype entropy_const = -6.8249;
  static mmdb::realtype entropy_rbody = 0.912179;
  static mmdb::realtype entropy_surf  = 0.104834;

  void writeEntropyFactors ( mmdb::io::RFile f )  {
  mmdb::byte Version=1;
    f.WriteByte ( &Version       );
    f.WriteReal ( &entropy_const );
    f.WriteReal ( &entropy_rbody );
    f.WriteReal ( &entropy_surf  );
  }

  void readEntropyFactors ( mmdb::io::RFile f )  {
  mmdb::byte Version;
    f.ReadByte ( &Version       );
    f.ReadReal ( &entropy_const );
    f.ReadReal ( &entropy_rbody );
    f.ReadReal ( &entropy_surf  );
  }

  void makeEntropyFactorsCIF ( mmdb::mmcif::PStruct mmCIFStruct )  {
    mmCIFStruct->PutReal ( entropy_const,"entropy_const",10 );
    mmCIFStruct->PutReal ( entropy_rbody,"entropy_rbody",10 );
    mmCIFStruct->PutReal ( entropy_surf ,"entropy_surf" ,10 );
  }

  void getEntropyFactorsCIF ( mmdb::mmcif::PStruct mmCIFStruct )  {
    mmCIFStruct->GetReal ( entropy_const,"entropy_const",false );
    mmCIFStruct->GetReal ( entropy_rbody,"entropy_rbody",false );
    mmCIFStruct->GetReal ( entropy_surf ,"entropy_surf" ,false );
  }

  void SetEntropyFactors ( mmdb::realtype ent_const,
                           mmdb::realtype ent_rbody,
                           mmdb::realtype ent_surf ) {
    entropy_const = ent_const;
    entropy_rbody = ent_rbody;
    entropy_surf  = ent_surf;
  }

  void GetEntropyFactors ( mmdb::realtype & ent_const,
                           mmdb::realtype & ent_rbody,
                           mmdb::realtype & ent_surf ) {
    ent_const = entropy_const;
    ent_rbody = entropy_rbody;
    ent_surf  = entropy_surf;
  }


  //  =========================  MonRef  ============================

  void MonRef::Init()  {
    M = NULL;
    i = 0;
    j = 0;
    k = 0;
    c = false;  // used for specifying internal interfaces
  }

  bool MonRef::isInternal ( RMonRef L )  {
  //   This function returns true if interface L is internal
  // for multimer M->U
  PMultimer U;
  int       di,dj,dk;
    L.c = false;
    if (L.M)  {
      U = M->U;
      if (U==L.M->U)  {
        // monomers M and L.M are found in the same multimer, now
        // check that the unit cell placement of M and L.M corresponds
        // to the interface
        di = U->R[L.M->mx].i - U->R[M->mx].i;
        dj = U->R[L.M->mx].j - U->R[M->mx].j;
        dk = U->R[L.M->mx].k - U->R[M->mx].k;
        L.c = ((L.i==di) && (L.j==dj) && (L.k==dk));
      }
    }
    return L.c;
  }

  void MonRef::GetTMatrix ( mmdb::mat44 & T, mmdb::mat44 & rom ) {
  int ii,jj;
    for (ii=0;ii<4;ii++)
      for (jj=0;jj<4;jj++)
        T[ii][jj] = M->uct[ii][jj];
    T[0][3] += i*rom[0][0] + j*rom[0][1] + k*rom[0][2];
    T[1][3] += i*rom[1][0] + j*rom[1][1] + k*rom[1][2];
    T[2][3] += i*rom[2][0] + j*rom[2][1] + k*rom[2][2];
  }


  //  =========================  Monomer  ==========================

  Monomer::Monomer ( int monId, int monNCSParent, int monDClass,
                       mmdb::mat44 & TMatrix, int nInterfaces )  {
  UNUSED_ARGUMENT(nInterfaces);
  int i,j;

    id        = monId;        // monomer serial number in ASU
    ncsParent = monNCSParent; // monomer ncs parent == domain->ncsParent
    ix        = 0;            // monomer serial number in unit cell
    mx        = 0;            // monomer serial number in multimer
    id2       = 0;            // all functions should leave it 0
    symOpNo   = 0;            // symmetry operation (ASU) number
    dclass    = monDClass;    // protein, dna/rna or ligand

    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        uct[i][j] = TMatrix[i][j];

    for (i=0;i<_max_n_int;i++)
      L[i] = NULL;

    U = NULL;

  }

  Monomer::~Monomer()  {
  int i;
    for (i=0;i<_max_n_int;i++)
      if (L[i])  delete[] L[i];
  }

  bool repeated_assignment = false;

  void Monomer::AddInterface ( PMonomer Mon, int interfaceNo,
                                int i, int j, int k, int nInterfaces ) {
  int n,m;
    n = 0;
    while (n<_max_n_int)  {
      if (!L[n])  {
        L[n] = new MonRef[nInterfaces];
        for (m=0;m<nInterfaces;m++)
          L[n][m].Init();
      }
      if (!L[n][interfaceNo].M)  {
        L[n][interfaceNo].M = Mon;
        L[n][interfaceNo].i = i;
        L[n][interfaceNo].j = j;
        L[n][interfaceNo].k = k;
        n = _max_n_int;
      } else if ((n>=1) && (dclass!=DCLASS_Ligand) &&
                           (Mon->dclass!=DCLASS_Ligand))  {
        repeated_assignment = true;
        n = _max_n_int;
      } else
        n++;
    }
  }

  bool Monomer::isParallel ( PMonomer Mon )  {
  int     i,j;
  bool parallel;
    parallel = true;
    for (i=0;(i<3) && parallel;i++)
      for (j=0;(j<3) && parallel;j++)
        parallel = fabs(uct[i][j]-Mon->uct[i][j])<0.00001;
    return parallel;
  }

  void  Monomer::getMassCenter ( PDomains D,
                                 mmdb::realtype & x, mmdb::realtype & y,
                                 mmdb::realtype & z )  {
  PDomain domain;
    domain = D->domain[ncsParent];
    x = uct[0][0]*domain->mx + uct[0][1]*domain->my +
        uct[0][2]*domain->mz + uct[0][3];
    y = uct[1][0]*domain->mx + uct[1][1]*domain->my +
        uct[1][2]*domain->mz + uct[1][3];
    z = uct[2][0]*domain->mx + uct[2][1]*domain->my +
        uct[2][2]*domain->mz + uct[2][3];
  }


  // ==========================  Multimer  ==========================

  Multimer::Multimer ( int maxSize )  {
  int i;

    if (maxSize>0)  {
      R = new MonRef[maxSize];  // list of monomers in the multimer
      for (i=0;i<maxSize;i++)
        R[i].Init();
    } else
      R = NULL;
    mSize = 0;    // multimer size

    dissEn      = 0.0;   // maximal free energy of dissociation
    entropy     = 0.0;   // entropy change at dissociation
    dissIntArea = 0.0;   // dissociation interface area
    dissEn0     = 0.0;   // ground-level free energy of dissociation
    entropy0    = 0.0;   // ground-level entropy change at dissociation
    asa         = 0.0;   // accessible surface area
    bsa         = 0.0;   // buried surface area
    seGain      = 0.0;   // solvation energy gain
    dissMult    = NULL;  // dissociation subunits ids
    nUC         = 0;     // number of these multimers in unit cell
    type        = 0;     // multimer type id in multimer set
    stable      = false;

  }


  Multimer::~Multimer()  {
    FreeMemory();
  }

  void Multimer::FreeMemory()  {

    if (R)  {
      delete[] R;
      R = NULL;
    }
    mSize = 0;

    mmdb::FreeVectorMemory ( dissMult,0 );

  }

  void Multimer::MakeMultimer ( PMonomer Monomer )  {
    R[0].M     = Monomer;
    R[0].M->U  = this;
    R[0].M->mx = 0;
    mSize = 1;
  }

  void Multimer::MakeMultimer ( PPMonomer Monomer, int nMonomers )  {
  int i;
    if (!R)  {
      R = new MonRef[nMonomers];  // list of monomers in the multimer
      for (i=0;i<nMonomers;i++)
        R[i].Init();
    }
    mSize = nMonomers;
    for (i=0;i<mSize;i++)  {
      R[i].M     = Monomer[i];
      R[i].M->U  = this;
      R[i].M->mx = i;
    }
  }

  void Multimer::MakeBoundMultimer ( PPMonomer Monomer,
                                     int       nMonomers )  {
  int i;
    if (!R)  {
      R = new MonRef[nMonomers];  // list of monomers in the multimer
      for (i=0;i<nMonomers;i++)
        R[i].Init();
    }
    mSize = 0;
    for (i=0;i<nMonomers;i++)
      if (Monomer[i]->L[0])  {
        R[mSize].M     = Monomer[i];
        R[mSize].M->U  = this;
        R[mSize].M->mx = mSize;
        mSize++;
      }
  }

  void Multimer::SetMultimerReferences()  {
  int i;
    for (i=0;i<mSize;i++)  {
      R[i].M->U  = this;
      R[i].M->mx = i;
    }
  }

  bool Multimer::isMacroMolecule()  {
  int     i;
  bool isMM;
    isMM = false;
    for (i=0;(i<mSize) && (!isMM);i++)
      isMM = (R[i].M->dclass!=DCLASS_Ligand);
    return isMM;
  }

  int Multimer::getMMSize()  {
  int i,mmSize;
    mmSize = 0;
    for (i=0;i<mSize;i++)
      if (R[i].M->dclass!=DCLASS_Ligand)
        mmSize++;
    return mmSize;
  }


  int Multimer::getMaxEType ( mmdb::realtype _entropy,
                              mmdb::realtype dissThreshold,
                              int            r,
                              mmdb::ivector  eType,
                              mmdb::ivector  iflag,
                              PPInterface    interface,
                              int            nInterfaces )  {
  mmdb::realtype intfEn,ien;
  int      i,j,itp,k,n;

    if (stable)  return mmdb::MaxInt;

    intfEn = 0.0;
    k      = -1;

    for (i=0;i<mSize;i++)
      R[i].M->id2 = 1;

    for (i=nInterfaces-1;(i>=0) && (k<0);i--)  {
      itp = eType[i];
      if ((itp>=r) && (iflag[itp]<0) &&
          (iflag[itp]!=INTF_Void))  {
        ien = interface[i]->stabEn +
              2.0*entropy_surf*interface[i]->intArea/1000.0;
        if (ien<0.0)  {
          //  summate only on affine interfaces
          for (j=0;j<mSize;j++)  {
            n = 0;
            while (n<_max_n_int)
              if (R[j].M->L[n])  {
                if (R[j].M->L[n][i].M)  {
                  if (R[j].M->L[n][i].M->id2!=1)  intfEn += ien;
                }
                n++;
              } else
                n = _max_n_int;
          }
        }
      }
      if (intfEn+_entropy<dissThreshold)  k = itp;
    }

    for (i=0;i<mSize;i++)
      R[i].M->id2 = 0;

    return k;

  }

  void Multimer::CalcInternalInterfaces ( int nInterfaces )  {
  //  Multimer references and index mx MUST be set
  PMonomer M1;
  int      i,j,n;

    for (i=0;i<mSize;i++)
      for (j=0;j<nInterfaces;j++)  {
        n = 0;
        while (n<_max_n_int)
          if (R[i].M->L[n])  {
            R[i].M->L[n][j].c = false;
            M1 = R[i].M->L[n][j].M;
            if (M1)  {
              if (M1->U==this)  {
                // check that L[n][j] links to the right unit cell
                R[i].M->L[n][j].c  =
                         ((R[i].i+R[i].M->L[n][j].i-R[M1->mx].i)==0) &&
                         ((R[i].j+R[i].M->L[n][j].j-R[M1->mx].j)==0) &&
                         ((R[i].k+R[i].M->L[n][j].k-R[M1->mx].k)==0);
              }
            }
            n++;
          } else
            n = _max_n_int;

      }

  }

  void Multimer::SetMultimer ( PMultimer     U,
                               int           nInterfaces,
                               PPInterface   Interface,
                               mmdb::ivector eType )  {
  UNUSED_ARGUMENT(nInterfaces);
  UNUSED_ARGUMENT(Interface);
  UNUSED_ARGUMENT(eType);
  //   Makes a copy of multimer U. This function does not require
  // that multimer references in monomers are set up.
  int  i;

    FreeMemory();

    mSize = U->mSize;  // multimer size
    nUC   = 0;

    R = new MonRef[mSize];  // list of monomers in multimer
    for (i=0;i<mSize;i++)  {
      R[i].M = U->R[i].M;
      R[i].i = U->R[i].i;
      R[i].j = U->R[i].j;
      R[i].k = U->R[i].k;
    }

  }

  #define _frame_dif_threshold 3.5

  bool Multimer::isIdentical ( PMultimer U, PDomains domains,
                               mmdb::mat44 & rom, bool withLigands )  {
  //    When withLigands is set false, the function compares
  // the number, composition and spatial arrangement of
  // macromolecular units only.
  PPDomain      D;
  PDFrame     * df2;
  DFrame        df1;
  mmdb::ivector dt1,dt2;
  mmdb::ovector c1,c2;
  mmdb::mat44   Ti,Tj,Ti1,Tj1;
  int           i,j,t0,i1,j1,ii,jj;
  bool          B;

    D = domains->domain;

    if (withLigands)
      if (U->mSize!=mSize) return false;
    if ((mSize==1) && (U->mSize==1))
      return (D[U->R[0].M->ncsParent]->type==D[R[0].M->ncsParent]->type);

    B = true;

    mmdb::GetVectorMemory ( c1,mSize   ,0 );
    mmdb::GetVectorMemory ( c2,U->mSize,0 );

    if (withLigands)  {
      for (i=0;i<mSize;i++)  {
        c1[i] = true;
        c2[i] = true;
      }
    } else  {
      i1 = 0;
      for (i=0;i<mSize;i++) {
        c1[i] = (R[i].M->dclass!=DCLASS_Ligand);
        if (c1[i]) i1++;
      }
      j1 = 0;
      for (i=0;i<U->mSize;i++)  {
        c2[i] = (U->R[i].M->dclass!=DCLASS_Ligand);
        if (c2[i]) j1++;
      }
      if (i1!=j1) return false;
    }


    mmdb::GetVectorMemory ( dt1,domains->nDTypes,1 );
    mmdb::GetVectorMemory ( dt2,domains->nDTypes,1 );

    for (i=1;i<=domains->nDTypes;i++)  {
      dt1[i] = 0;
      dt2[i] = 0;
    }

    for (i=0;i<mSize;i++)
      if (c1[i])
        dt1[D[R[i].M->ncsParent]->type]++;

    for (i=0;i<U->mSize;i++)
      if (c2[i])
        dt2[D[U->R[i].M->ncsParent]->type]++;

    t0 = 1;
    j  = mmdb::MaxInt;
    i1 = 0;
    for (i=1;(i<=domains->nDTypes) && B;i++)  {
      B = (dt1[i]==dt2[i]);
      if (B && (dt1[i]>0))  {
        i1 += dt1[i];
        if (dt1[i]<j)  {
          j  = dt1[i];
          t0 = i;
        }
      }
    }

    mmdb::FreeVectorMemory ( dt1,1 );
    mmdb::FreeVectorMemory ( dt2,1 );

    if (B && (i1>1))  {

      df2 = new PDFrame[U->mSize];
      for (i=0;i<U->mSize;i++)
        df2[i] = NULL;

      B = false;
      for (i=0;(i<mSize) && (!B);i++)
        if (D[R[i].M->ncsParent]->type==t0)  {
          R[i].GetTMatrix ( Ti,rom );
          for (i1=0;(i1<U->mSize) && (!B);i1++)
            if (D[U->R[i1].M->ncsParent]->type==t0)  {
              U->R[i1].GetTMatrix ( Ti1,rom );
              if (!df2[i1])  {
                df2[i1] = new DFrame[U->mSize];
                for (j1=0;j1<U->mSize;j1++)
                  if (c2[j1] && (j1!=i1))  {
                    U->R[j1].GetTMatrix ( Tj1,rom );
                    MakeDFrame ( df2[i1][j1],
                                 D[U->R[i1].M->ncsParent],Ti1,
                                 D[U->R[j1].M->ncsParent],Tj1 );
                  }
              }
              B = true;
              for (j=0;(j<mSize) && B;j++)
                if (c1[j] && (j!=i))  {
                  R[j].GetTMatrix ( Tj,rom );
                  MakeDFrame ( df1,D[R[i].M->ncsParent],Ti,
                                   D[R[j].M->ncsParent],Tj );
                  B = false;
                  for (j1=0;(j1<U->mSize) && (!B);j1++)
                    if (c2[j1] && (j1!=i1))  {
                      B = true;
                      for (ii=0;(ii<frameLen) && B;ii++)
                        for (jj=0;(jj<frameLen) && B;jj++)
                          B = fabs(df1[ii][jj]-df2[i1][j1][ii][jj]) <
                              _frame_dif_threshold;
                    }
                }
            }
        }

      for (i=0;i<U->mSize;i++)
        if (df2[i])  delete[] df2[i];
      delete[] df2;

    }

    mmdb::FreeVectorMemory ( c2,0 );
    mmdb::FreeVectorMemory ( c1,0 );

    return B;

  }

  bool Multimer::isComplemenatry ( PMultimer U )  {
  //   This function identifies whether multimers U and 'this' may
  // potentially "complement" each other in ASU. The multimers may
  // "complement" each other if (a) they have different size or
  // (b) they include different monomers in original positions
  int     i,j;
  bool complement;
    if (U->mSize!=mSize) return true;
    complement = true;
    for (i=0;(i<mSize) && complement;i++)
      for (j=0;(j<mSize) && complement;j++)
        if (R[i].M->id==U->R[j].M->id) complement = false;
    return complement;
  }

  void  Multimer::CalcProperties ( PDomains      D,
                                   PPInterface   interface,
                                   mmdb::ivector eType,
                                   int           nInterfaces,
                                   int           nETypes,
                                   PSymNumber    SNum,
                                   mmdb::mat44 & rom )  {
  // The function calculates asa, bsa (accessible surface area and
  // buried surface area) and solvation energy gain (seGain).
  int  i,j,k,n;

    asa    = 0.0;
    bsa    = 0.0;
    seGain = 0.0;

    for (i=0;i<mSize;i++)
      asa += D->domain[R[i].M->id]->surfArea;

    for (i=0;i<nInterfaces;i++)  {
      k = 0;
      for (j=0;j<mSize;j++)  {
        n = 0;
        while (n<_max_n_int)
          if (R[j].M->L[n])  {
            if (R[j].M->L[n][i].c)  k++;
            n++;
          } else
            n = _max_n_int;
      }
      asa    -= k*interface[i]->intArea;
      bsa    += k*interface[i]->intArea;
      seGain += k*interface[i]->intDeltaG;
    }

    seGain /= 2.0;

    Dissociate ( D,interface,eType,nInterfaces,nETypes,SNum,rom );

    if (dissIntArea>0.0)  stable  = ((dissEn+entropy)<=0.0);
                    else  stable  = true;

  }


  void Multimer::calcRBETerms ( int              n,
                                mmdb::realtype & logW,
                                mmdb::realtype & logJ,
                                RDIRecursion     P )  {
  //   This function calculates the rigid-body entropy terms for
  // dissociating part number n (n=1,2,...) of the multimer:
  //     weight    - molecular weight
  //     pJ        - product of moments of inertia in terms of
  //                 molecular weight
  //   If n equals 0, then entropy terms of the multimer itself
  // are calculated.
  //   Other parameters:
  //     D   pointer to the array of domains
  //     mc  array of 3-vectors with mass centers of the monomeric
  //         units
  mmdb::rmatrix  A,T;
  mmdb::rvector  Eigen,Aik;
  mmdb::mat33    itm;
  mmdb::realtype weight,mx,my,mz, w;
  int            i,j;

    //  1.  Calculate the multimer's weight and mass center

    mx     = 0.0;  // mass center
    my     = 0.0;  //  of the dissociating
    mz     = 0.0;  //    part #n
    weight = 0.0;
    for (i=0;i<mSize;i++)
      if ((!n) || (R[i].M->id2==n))  {
        w       = P.D[R[i].M->id]->weight;
        weight += w;
        mx     += w*P.mc[i][0];
        my     += w*P.mc[i][1];
        mz     += w*P.mc[i][2];
      }

    if (weight<=0.0)  {
      logW = -mmdb::MaxReal;
      logJ = 0.0;
      return;
    }

    mx /= weight;
    my /= weight;
    mz /= weight;

    //  2.  Calculate the multimer's tensor of inertia

    P.SNum->getStorage ( A,T,Eigen,Aik );

    for (i=1;i<=3;i++)
      for (j=i;j<=3;j++)
        A[i][j] = 0.0;

    for (i=0;i<mSize;i++)
      if ((!n) || (R[i].M->id2==n))  {
        P.D[R[i].M->id]->getInertiaTensor ( itm,P.TMatrix[i],mx,my,mz );
        A[1][1] += itm[0][0];
        A[1][2] += itm[0][1];
        A[1][3] += itm[0][2];
        A[2][2] += itm[1][1];
        A[2][3] += itm[1][2];
        A[3][3] += itm[2][2];
      }

    A[2][1] = A[1][2];
    A[3][1] = A[1][3];
    A[3][2] = A[2][3];

    //  3.  Diagonalise the inertia tensor and get the product
    //      of principal moments of inertia

    for (i=1;i<=3;i++)
      for (j=1;j<=3;j++)
        A[i][j] /= weight;

    mmdb::math::Jacobi ( 3,A,T,Eigen,Aik,i );

    logW = log(weight);
    logJ = 3.0*logW + log(Eigen[1]) + log(Eigen[2]) + log(Eigen[3]);

  }


  void Multimer::calcRBEntropy ( int nDiss, mmdb::realtype & entropy_change,
                                  RDIRecursion P )  {
  //  This function calculates the rigid-body entropy *change* upon
  // dissociation of the multimer into nDiss parts. Each part is
  // identified by index id2 of the monomeric units, 1<=id2<=nDiss.
  mmdb::realtype logW,logJ, eTrans,eRot, symNo;
  int      i;

    if (nDiss<=1)  {
      entropy = 0;
      return;
    }

    eTrans = 0.0;  // translation entropy contribution
    eRot   = 0.0;  // rotation entropy contribution

    //  1.  Calculate entropy of dissociated parts less constant terms

    for (i=1;i<=nDiss;i++)  {
      calcRBETerms ( i,logW,logJ,P );
      if (logW>-mmdb::MaxReal)  {
        eTrans += logW;
        symNo   = P.SNum->getSymNumber ( this,i,P.domains,P.TMatrix );
        eRot   += logJ - 2.0*log(symNo);
      }
    }

    //  2.  Substract entropy of the multimer less constant term

    calcRBETerms ( 0,logW,logJ,P );
    if (logW>-mmdb::MaxReal)  {
      eTrans -= logW;
      symNo   = P.SNum->getSymNumber ( this,0,P.domains,P.TMatrix );
      eRot   -= logJ - 2.0*log(symNo);
    }

    //  3.  Sum the terms up

    entropy_change = entropy_const*(nDiss-1) +
                     entropy_rbody*TRconst*(3.0*eTrans+eRot)/2.0;

  //  coincides with the old one - CHECKED:
  //  entropy_change = entropy_const*(nDiss-1) + entropy_rbody*eTrans;

  }


  // interface type flag for fixed ligand interfaces
  #define  _itp_flig    20000

  void Multimer::dissInterface (
                int           typeNo,  // interface type to disengage
                mmdb::ivector       itp,     // list of all interface types
                RDIRecursion P,
                mmdb::realtype &    dissEnX, // minimal dissociation free energy
                mmdb::realtype &    entropyX,
                mmdb::realtype &    dissIntAreaX,
                mmdb::ivector       dissMultX )  {
  PMultimer U;
  mmdb::ivector    id2,itp1,dissMult1;
  mmdb::realtype   dissEn1,entropy1,dissIntArea1;
  mmdb::realtype   dEn,dEntropy,dSAS;
  int        i,j,k,k0,n,m,ti,t,tNo, engaged;
  bool    st;

    if (typeNo>=0)  {

      tNo = typeNo;
      while (tNo<P.nit)
        if ((itp[tNo]>=0) && (itp[tNo]<_itp_flig))  break;
                                              else  tNo++;

      if (tNo>=P.nit)  return;

      // disengage interfaces of type tNo: negative interface type
      // means that the interface is disengaged
      itp[tNo] = -(itp[tNo]+1);

    } else if (typeNo==-2)  {
      // disengage all interfaces. This is used for the measurement
      // of the "absolute" free energy of complex

      for (i=0;i<P.nit;i++)
        if ((itp[i]>=0) && (itp[i]<_itp_flig))
          itp[i] = -(itp[i]+1);

      tNo = -2;

    } else
      tNo = -1;

    //   At this point, all monomers are marked 0. We shall see
    // whether disengaging the interface splits the multimer
    // into two or more parts. Start with first monomer and
    // mark it with 1. Then mark with 1 all monomers linked by
    // engaged (closed) interfaces and count them. If number of
    // linked monomers is less than multimer size, the multimer
    // splits.

    R[0].M->id2 = 1;  // dissociating unit id
    k = 1;            // k will count linked monomers
    do {
      k0 = k;
      for (i=0;i<P.nInterfaces;i++)  {
        ti = P.eType[i];
        t  = -(ti+1);
        // check whether ith interface is internal and engaged
        engaged = 0;  // if remains 0 then external interface
        for (j=0;(j<P.nit) && (!engaged);j++)
          if ((ti==itp[j]) || (ti==itp[j]-_itp_flig))
                               engaged =  1; // internal and engaged
          else if (t==itp[j])  engaged = -1; // int-l and disengaged
        if (engaged==1)
          for (j=0;j<mSize;j++)
            if (R[j].M->id2==0)  {  // consider only unmarked monomers
              m = 0;
              while (m<_max_n_int)
                if (R[j].M->L[m])  {
                  if (R[j].M->L[m][i].c)  {
                    if (R[j].M->L[m][i].M->id2==1)  {
                      R[j].M->id2 = 1;
                      m = _max_n_int;
                    }
                  }
                  m++;
                } else
                  m = _max_n_int;
              if (R[j].M->id2==1)  {
                // jth monomer is linked, mark and count it
                k++;
              }
            }
      }
    } while (k>k0);

    if (k<mSize)  {
      // assembly splits into at least 2 parts - check energetics
      // and stop the recursion if all parts are stable

      // identify all dissociated parts
      n = 1;  // this will be mark number for the split parts
      do  {
        n++;
        // mark first yet unmarked monomer with number n
        k = -1;
        for (i=0;(i<mSize) && (k<0);i++)
          if (R[i].M->id2==0)  {
            R[i].M->id2 = n;
            k = i;
          }
        if (k>=0)  {
          k = 1;
          do {
            k0 = k;
            for (i=0;i<P.nInterfaces;i++)  {
              ti = P.eType[i];
              t  = -(ti+1);
              engaged = 0;  // if remains 0 then external interface
              for (j=0;(j<P.nit) && (!engaged);j++)
                if ((ti==itp[j]) || (ti==itp[j]-_itp_flig))
                                     engaged =  1; // engaged
                else if (t==itp[j])  engaged = -1; // disengaged
              if (engaged==1)
                for (j=0;j<mSize;j++)
                  if (R[j].M->id2==0)  {
                    m = 0;
                    while (m<_max_n_int)
                      if (R[j].M->L[m])  {
                        if (R[j].M->L[m][i].c)  {
                          if (R[j].M->L[m][i].M->id2==n)  {
                            R[j].M->id2 = n;
                            m = _max_n_int;
                          }
                        }
                        m++;
                      } else
                        m = _max_n_int;
                    if (R[j].M->id2==n)
                      k++;
                  }
            }
          } while (k>k0);  // until no new counts are recorded

        }

      } while (k>=0);

      n--;  // number of dissociated parts

      // check whether all dissociated parts are stable

      mmdb::GetVectorMemory ( id2      ,mSize,0 );
      mmdb::GetVectorMemory ( dissMult1,mSize,0 );
      mmdb::GetVectorMemory ( itp1     ,P.nit,0 );
      U = new Multimer ( mSize );

      // save id2; R[].M->id2 are reused in recursive calls
      for (i=0;i<mSize;i++)
        id2[i] = R[i].M->id2;

      for (i=0;i<P.nit;i++)
        itp1[i] = itp[i];

      st = true;
      for (k=1;(k<=n) && st;k++)  {

        U->mSize = 0;
        for (i=0;i<mSize;i++)
          if (id2[i]==k)  U->R[U->mSize++].M = R[i].M;

        if (U->mSize>1)  {

          for (i=0;i<mSize;i++)  {
            R[i].M->id2  = 0;
            dissMult1[i] = 1;
          }

          dissEn1      = -mmdb::MaxReal/2.0;
          entropy1     = -mmdb::MaxReal/2.0;
          dissIntArea1 = 0.0;

          U->dissInterface ( tNo+1,itp1,P,dissEn1,entropy1,
                             dissIntArea1,dissMult1 );

          st = (dissEn1+entropy1<0.0);

        }

      }

      mmdb::FreeVectorMemory ( dissMult1,0 );
      mmdb::FreeVectorMemory ( itp1     ,0 );
      delete U;

      if (st)  {

        // calculate entropy and free energy of dissociation

        for (i=0;i<mSize;i++)
          R[i].M->id2 = id2[i];

        mmdb::FreeVectorMemory ( id2,0 );

        dEn  = 0.0;
        dSAS = 0.0;
        for (i=0;i<P.nInterfaces;i++)  {
          k = 0;
          for (j=0;j<mSize;j++)  {
            m = 0;
            while (m<_max_n_int)
              if (R[j].M->L[m])  {
                if (R[j].M->L[m][i].c)  {
                  if ((R[j].M->L[m][i].M->id2>0) &&
                      (R[j].M->L[m][i].M->id2!=R[j].M->id2)) k++;
                }
                m++;
              } else
                m = _max_n_int;
          }
          dEn  += k*P.interface[i]->stabEn;
          dSAS += k*P.interface[i]->intArea;
        }

        dEn  /= 2.0;
        dSAS /= 2.0;

        calcRBEntropy ( n,dEntropy,P );
        dEntropy += 2.0*entropy_surf*dSAS/1000.0;

  //      printf ( " dSAS=%15.8g\n",dSAS );

        // Check for the minimal free energy of dissociation.
        // Note that dEn+dEntropy gives free energy of association,
        // so that the inequality sign in the following statement
        // is correct
        if (dEn+dEntropy>dissEnX+entropyX)  {
          dissEnX      = dEn;
          entropyX     = dEntropy;
          dissIntAreaX = dSAS;
          for (i=0;i<mSize;i++)
            dissMultX[i] = R[i].M->id2;
        }

        for (i=0;i<mSize;i++)
          R[i].M->id2 = 0;

      } else if (tNo<P.nit-1)  {

        mmdb::FreeVectorMemory ( id2,0 );
        dissInterface ( tNo+1,itp,P,dissEnX,entropyX,
                        dissIntAreaX,dissMultX );

      }

    } else  {

      // assembly does not split - try opening other interfaces

      for (i=0;i<mSize;i++)
        R[i].M->id2 = 0;

      if (tNo<P.nit-1)
        dissInterface ( tNo+1,itp,P,dissEnX,entropyX,
                        dissIntAreaX,dissMultX );
    }

    // restore the interface state before moving up the recursion
    if (tNo>=0)
      itp[tNo] = -(itp[tNo]+1);

  }


  #ifdef __debug
  extern void out4 ( mmdb::mat44 & T, pstr name );
  #endif

  void Multimer::CalcTMatrices ( mmdb::rpmat44 TMatrix,
                                 mmdb::mat44 & rom )  {
  //   This function calculates transformation matrices for all
  // monomeric units in the multimer. Array TMatrix is allocated
  // with dimension [0..mSize-1].
  int i,j,k;

    TMatrix = new mmdb::mat44[mSize];

    for (i=0;i<mSize;i++)  {

      for (j=0;j<4;j++)
        for (k=0;k<4;k++)
          TMatrix[i][j][k] = R[i].M->uct[j][k];
      TMatrix[i][0][3] += R[i].i*rom[0][0] + R[i].j*rom[0][1] +
                                             R[i].k*rom[0][2];
      TMatrix[i][1][3] += R[i].i*rom[1][0] + R[i].j*rom[1][1] +
                                             R[i].k*rom[1][2];
      TMatrix[i][2][3] += R[i].i*rom[2][0] + R[i].j*rom[2][1] +
                                             R[i].k*rom[2][2];
    }

  }


  void DIRecursion::Init ( PDomains         Doms,
                           PPInterface      Ints,
                           PSymNumber       symNum,
                           mmdb::ivector    intType,
                           int              nInts,
                           PMultimer        U,
                           int              nETypes,
                           mmdb::mat44    & rom,
                           mmdb::ivector  & itp )  {
  //   Makes the list of engaged interface types in multimer U
  // (itp) and of their stabilisation energies (iSE)
  mmdb::rvector iSE;
  mmdb::ovector lig;
  int           i,j,k,n;

    domains     = Doms;
    D           = Doms->domain;
    interface   = Ints;
    SNum        = symNum;
    eType       = intType;
    nInterfaces = nInts;

    mmdb::GetVectorMemory ( itp,nETypes,0 );
    mmdb::GetVectorMemory ( lig,nETypes,0 );
    mmdb::GetVectorMemory ( iSE,nETypes,0 );
    for (i=0;i<nETypes;i++)  {
      itp[i] = -1;    // initially disengaged
      lig[i] = false;
      iSE[i] = 0.0;
    }

    //  mark monomers belonging to multimer U
    for (i=0;i<U->mSize;i++)
      U->R[i].M->id2 = 1;

    for (i=0;i<nInterfaces;i++) {
      k = 0;
      for (j=0;j<U->mSize;j++)  {
        n = 0;
        while (n<_max_n_int)
          if (U->R[j].M->L[n])  {
            if (U->R[j].M->L[n][i].c)  {
              if (U->R[j].M->L[n][i].M->id2==1)  k++;
            }
            n++;
          } else
            n = _max_n_int;
      }
      if (k>0)  {  // interface internal to U, engage it
        j = eType[i];
        if (interface[i]->fixedLigand)  itp[j] = _itp_flig+j;
                                  else  itp[j] = j;
        lig[j] = (interface[i]->dclass1==DCLASS_Ligand) ||
                 (interface[i]->dclass2==DCLASS_Ligand);
        iSE[j] += k*(interface[i]->stabEn +
                  2.0*entropy_surf*interface[i]->intArea/1000.0)/2.0;
      }
    }

    // calculate the number of engaged interfaces
    nit = 0;
    for (i=0;i<nETypes;i++)
      if (itp[i]>=0)  nit++;

    for (i=nETypes-1;(i>=0) && (nit>1);i--)  {
      k = 0;
      for (j=0;(j<nInterfaces) && (k<=0);j++)
        if (eType[j]==itp[i])  {
          if (interface[j]->casual)  k = -1;
                               else  k =  1;
        }
      if (k<0)  {
        itp[i] = -1;
        nit--;
      }
    }

    //  move disengaged and insignificant interfaces to the end
    // of list, 'nit'  will be the number of significant engaged
    // interfaces
    nit = nETypes; // nit will be the number of significant intf types
    i   = 0;
    while (i<nit)
      if (itp[i]<0)  {
        nit--;
        if (i<nit)  {
          mmdb::ISwap ( itp[i],itp[nit] );
          mmdb::OSwap ( lig[i],lig[nit] );
          mmdb::RSwap ( iSE[i],iSE[nit] );
        }
      } else
        i++;

    // Move ligand interfaces to the end of list of significant
    // interfaces. This allows one to keep focus on macromolecular
    // dissociation.

    nMMInt = nit; // nMMInt will be the number of m'molecular interfaces
    i      = 0;
    while (i<nMMInt)
      if (lig[i])  {
        nMMInt--;
        if (i<nMMInt)  {
          mmdb::ISwap ( itp[i],itp[nMMInt] );
          mmdb::OSwap ( lig[i],lig[nMMInt] );
          mmdb::RSwap ( iSE[i],iSE[nMMInt] );
        }
      } else
        i++;

    // Sort significant interfaces in order of decreasing stabilization
    // (binding) energy so that the weakest interfaces come first.

    for (i=0;i<nMMInt;i++)
      for (j=i+1;j<nMMInt;j++)
        if (iSE[i]<iSE[j])  {
          mmdb::ISwap ( itp[i],itp[j] );
          mmdb::RSwap ( iSE[i],iSE[j] );
        }

    for (i=nMMInt;i<nit;i++)
      for (j=i+1;j<nit;j++)
        if (iSE[i]<iSE[j])  {
          mmdb::ISwap ( itp[i],itp[j] );
          mmdb::RSwap ( iSE[i],iSE[j] );
        }

    //  allow ligand dissociation if there is no macromolecular
    //  interfaces:
    if (nMMInt<=0)  nMMInt = nit;

    //  unmark monomers belonging to multimer U
    for (i=0;i<U->mSize;i++)
      U->R[i].M->id2 = 0;

    mmdb::FreeVectorMemory ( lig,0 );
    mmdb::FreeVectorMemory ( iSE,0 );


    U->CalcTMatrices ( TMatrix,rom );

    mc = new mmdb::vect3[U->mSize];

    for (i=0;i<U->mSize;i++)  {

      k = U->R[i].M->id;
      mc[i][0] = D[k]->mx;
      mc[i][1] = D[k]->my;
      mc[i][2] = D[k]->mz;
      mmdb::TransformXYZ ( TMatrix[i],mc[i][0],mc[i][1],mc[i][2] );

    }

  }


  void DIRecursion::FreeMemory()  {
    if (mc)      delete[] mc;
    if (TMatrix) delete[] TMatrix;
  }


  void Multimer::Dissociate ( PDomains      D,
                              PPInterface   interface,
                              mmdb::ivector eType,
                              int           nInterfaces,
                              int           nETypes,
                              PSymNumber    SNum,
                              mmdb::mat44 & rom )  {
  //
  //  eType[i] is the engagement type of interface i
  //
  DIRecursion   P;
  mmdb::ivector itp;
  mmdb::ivector itp0;
  int            i,i1;

    // Re-initialise subunit ids
    mmdb::FreeVectorMemory ( dissMult,0 );
    mmdb::GetVectorMemory  ( dissMult,mSize,0 );
    for (i=0;i<mSize;i++)
      dissMult[i] = 0;

    if (mSize<=1)  {
      dissEn      = 0.0;
      entropy     = 0.0;
      dissIntArea = 0.0;
      dissEn0     = 0.0;
      entropy0    = 0.0;
      dissMult[0] = 1;
      return;
    }

    P.Init ( D,interface,SNum,eType,nInterfaces,this,nETypes,rom,itp );

    // Calculate the ground-level free energy first
    mmdb::GetVectorMemory ( itp0,nETypes,0 );
    for (i=0;i<nETypes;i++)
      itp0[i] = itp[i];

    dissEn0     = -mmdb::MaxReal/2.0;  // make it infinitely low initial
    entropy0    = -mmdb::MaxReal/2.0;  //   free energy of dissociation
    dissIntArea = 0.0;

    dissInterface ( -2,itp,P,dissEn0,entropy0,dissIntArea,dissMult );

  //  printf ( " %5i dia=%15.8g  bsa=%15.8g ",mSize,dissIntArea,bsa );

    if (dissIntArea<=0.0)  {
      dissEn0  = 0.0;
      entropy0 = 0.0;
    }

    // Now calculate dissociation pattern and maximal free energy
    // of dissociation

    for (i=0;i<nETypes;i++)
      itp[i] = itp0[i];

    for (i=0;i<mSize;i++)
      dissMult[i] = 0;

    dissEn      = -mmdb::MaxReal/2.0;  // make it infinitely low initial
    entropy     = -mmdb::MaxReal/2.0;  //   free energy of dissociation
    dissIntArea = 0.0;

    //   If P.nMMInt<nETypes, it means that there are insignificant
    // interfaces that are considered as non-engageable. In this
    // case, set i1=-1, which will cause dissInterface to check
    // whether these insignificant interfaces split the assembly
    // apart. Otherwise, set i1=0 and dissInterface will try to
    // split the assembly by disengaging only the real interfaces.
    if (P.nMMInt<nETypes)  i1 = -1;
                     else  i1 = 0;

    for (i=i1;i<P.nMMInt;i++)
      dissInterface ( i,itp,P,dissEn,entropy,dissIntArea,dissMult );

  //  printf ( "  dia1=%15.8g\n",dissIntArea );

    if (dissIntArea<=0.0)  {
      dissEn      = 0.0;
      entropy     = 0.0;
      dissIntArea = 0.0;
      for (i=0;i<mSize;i++)
        dissMult[i] = 1;
    }

    mmdb::FreeVectorMemory ( itp ,0 );
    mmdb::FreeVectorMemory ( itp0,0 );
    P.FreeMemory();

  }

  bool Multimer::isInternal ( int intfType, mmdb::ivector eType,
                                  int nInterfaces )  {
  int     i,j,n;
  bool internal;
    internal = false;
    for (i=0;(i<nInterfaces) && (!internal);i++)
      if (eType[i]==intfType)  {
        for (j=0;(j<mSize) && (!internal);j++)  {
          n = 0;
          while ((n<_max_n_int) && (!internal))
            if (R[j].M->L[n])  {
              if (R[j].M->L[n][i].M)
                internal = R[j].isInternal ( R[j].M->L[n][i] );
              n++;
            } else
              n = _max_n_int;
        }
      }
    return internal;
  }

  int Multimer::EngageInterface ( int            interfaceNo,
                                  mmdb::omatrix  parMon,
                                  mmdb::realtype intArea )  {
  UNUSED_ARGUMENT(intArea);
  //    Engages interface number interfaceNo for the multimer and
  //  merges the interfacing multimers found. Returns the number
  //  of merged interfaces or -1 if the interface cannot be engaged.
  //    This version of multimer engaging function checks on
  //  parallel monomers, which might include collaterality of
  //  physically equivalent, although not crystallographically
  //  identical, monomers.
  PMonRef   L1;
  PMultimer U1;
  int        i,j,k,m,n, dic,djc,dkc, rc;
    n  = 0;
    i  = 0;
    rc = 0;
    while ((i<mSize) && (!rc))  {
      m = 0;
      while ((m<_max_n_int) && (!rc))
        if (R[i].M->L[m])  {
          L1 = &(R[i].M->L[m][interfaceNo]);
          if (L1->M)  {
            n++;
            U1 = L1->M->U;
            if (U1->mSize>0)  {
              if (U1==this)  {
                // check that interfacing monomer is found in the
                // correct unit cell. If the unit cell is wrong then
                // engaging the interface makes an infinite multimer
                dic = R[i].i + L1->i - R[L1->M->mx].i;
                djc = R[i].j + L1->j - R[L1->M->mx].j;
                dkc = R[i].k + L1->k - R[L1->M->mx].k;
                if (dic || djc || dkc)  rc = -2;
                                  else  L1->c = true; // internal intfce
              } else  {
                // check for parallel molecules
                for (j=0;(j<mSize) && (!rc);j++)
                  for (k=0;(k<U1->mSize) && (!rc);k++)
                    if (parMon[R[j].M->ix][U1->R[k].M->ix]) rc = -1;
                if (!rc)  {
                  // assemblies can be merged, so merge them
                  dic = R[i].i + L1->i - U1->R[L1->M->mx].i;
                  djc = R[i].j + L1->j - U1->R[L1->M->mx].j;
                  dkc = R[i].k + L1->k - U1->R[L1->M->mx].k;
                  for (j=0;j<U1->mSize;j++)  {
                    R[mSize].M     = U1->R[j].M;
                    R[mSize].M->U  = this; // monomer belongs to
                                           //   'this' now
                    R[mSize].M->mx = mSize;
                    R[mSize].i     = U1->R[j].i + dic;
                    R[mSize].j     = U1->R[j].j + djc;
                    R[mSize].k     = U1->R[j].k + dkc;
                    mSize++;
                  }
                  L1->c     = true; // internal interface
                  U1->mSize = 0;    // monomer U1 is now empty
                }
              }
            }
          }
          m++;
        } else
          m = _max_n_int;
      i++;
    }

    if (rc)  return rc;

    return n;

  }


  void  Multimer::getEntropyTerms ( mmdb::realtype & const_ent,
                                     mmdb::realtype & rbody_ent,
                                     mmdb::realtype & surf_ent )  {
  int  i,n;

    surf_ent = 2.0*dissIntArea/1000.0;

    if ((!dissMult) || (mSize<=1))  {
      const_ent = 0.0;
      rbody_ent = 0.0;
      return;
    }

    n = 0;
    for (i=0;i<mSize;i++)
      n = mmdb::IMax ( n,dissMult[i] );

    const_ent = n - 1;

    if (entropy_rbody!=0.0)  {
      rbody_ent  = entropy - entropy_const*const_ent -
                             entropy_surf*surf_ent;
      rbody_ent /= entropy_rbody;
    } else
      rbody_ent = 0.0;

  }


  void Multimer::GetDissProperties ( PPInterface      interface,
                                     int              nInterfaces,
                                     mmdb::rvector    DeltaSAS,
                                     mmdb::realtype & DeltaG,
                                     int            & nHBonds,
                                     int            & nSBridges )  {
  //  This function returns the total solvation energy gain,
  // number of hydrogen bonds and solvent bridges for the interfaces
  // between dissociated parts of the multimer. The function does
  // require the multimer references to be set prior calling, and
  // internal interfaces to be calculated:
  //    SetMultimerReferences  ();
  //    CalcInternalInterfaces ( nInterfaces );
  int i,j,k,ii,n, a1,a2;

    for (ii=0;ii<nASPs;ii++)
      DeltaSAS[ii] = 0.0;
    DeltaG    = 0.0;
    nHBonds   = 0;
    nSBridges = 0;

    if (dissMult && (mSize>1))  {

      for (i=0;i<mSize;i++)
        R[i].M->id2 = dissMult[i];

      for (i=0;i<mSize;i++)  {
        a1 = R[i].M->id2;
        for (j=0;j<nInterfaces;j++)  {
          k = 0;
          n = 0;
          while (n<_max_n_int)
            if (R[i].M->L[n])  {
              if (R[i].M->L[n][j].c)  {
                a2 = R[i].M->L[n][j].M->id2;
                if ((a2>0) && (a2!=a1))  k++;
              }
              n++;
            } else
              n = _max_n_int;
          if (k)  {
            for (ii=0;ii<nASPs;ii++)
              DeltaSAS[ii] += k*(interface[j]->DeltaSAS1[ii] +
                                 interface[j]->DeltaSAS2[ii]);
            DeltaG    += k*interface[j]->intDeltaG;
            nHBonds   += k*interface[j]->nHBonds;
            nSBridges += k*interface[j]->nSBridges;
          }
        }
      }

      for (i=0;i<mSize;i++)
        R[i].M->id2 = 0;

    }

    for (ii=0;ii<nASPs;ii++)
      DeltaSAS[ii] /= 2.0;
    DeltaG    /= 2.0;
    nHBonds   /= 2;
    nSBridges /= 2;

  }


  //  ======================  IntfRecStack  =========================

  IntfRecStack::IntfRecStack ( int nETypes, int maxNofMultimers )  {
    mmdb::GetVectorMemory ( iflag,nETypes        ,0 );
    mmdb::GetVectorMemory ( mSize,maxNofMultimers,0 );
    U = new PMultimer[maxNofMultimers];
    nMultimers = 0;
  }

  IntfRecStack::~IntfRecStack()  {
    mmdb::FreeVectorMemory ( iflag,0 );
    mmdb::FreeVectorMemory ( mSize,0 );
    if (U)  delete[] U;
  }


  //  ========================  MultimerSet  ========================

  MultimerSet::MultimerSet()  {
    intf          = NULL;
    U             = NULL;
    nMultimers    = 0;      // number of multimers obtained
    maxStableSize = 0;      // maximal size of stable multimer
    equiv_all     = false;  // true if all multimers are equivalent
    stable_all    = false;
  }

  MultimerSet::~MultimerSet()  {
    FreeMemory();
  }

  void MultimerSet::FreeMemory()  {
  int i;
    mmdb::FreeVectorMemory ( intf,0 );
    if (U)  {
      for (i=0;i<nMultimers;i++)
        if (U[i])  delete U[i];
      delete[] U;
      U = NULL;
    }
    nMultimers = 0;
  }

  int MultimerSet::getMaxEType ( int          r,
                                  mmdb::realtype     dissThreshold,
                                  mmdb::ivector      eType,
                                  PPInterface Interface,
                                  int          nInterfaces )  {
  mmdb::realtype entropy;
  int      i,k;

    k = -1;
    for (i=0;i<nMultimers;i++)  {
      entropy = U[i]->entropy;
      k = mmdb::IMax ( k,U[i]->getMaxEType(entropy,dissThreshold,r,eType,intf,
                   Interface,nInterfaces) );
    }

    return k;

  }

  bool MultimerSet::isEqual ( PIntfRecStack recStack,
                                  int nETypes )  {
  //   This function returns true if assembly set in recStack
  // has been obtained by engaging the same set of interfaces
  // as 'this' results.
  int     i;
  bool same;

    same = (intf!=NULL);
    for (i=0;(i<nETypes) && same;i++)
      same = ((intf[i]>INTF_Open)  && (recStack->iflag[i]>INTF_Open)) ||
             ((intf[i]<=INTF_Open) && (recStack->iflag[i]<=INTF_Open));

    return same;

  }


  bool multimer_overflow = false;

  void MultimerSet::makeSet ( PIntfRecStack recStack,
                              int           nETypes,
                              int           nInterfaces,
                              int           nMonomers,
                              PPInterface   interface,
                              PDomains      D,
                              mmdb::mat44 & rom,
                              mmdb::ivector eType )  {
  //   Collects assembly types found in recStack.  nInterfaces is total
  // number of interfaces in the system, nMonomers is number of
  // monomers in ASU.
  PMultimer U1;
  int       i,j,k;

    FreeMemory();

    mmdb::GetVectorMemory ( intf,nETypes,0 );
    for (i=0;i<nETypes;i++)
      intf[i] = recStack->iflag[i];

    //   As a consequence of crystal symmetry, the number of different
    // multimers cannot exceed the number of different monomers in ASU
    U = new PMultimer[nMonomers];
    for (i=0;i<nMonomers;i++)
      U[i] = NULL;

    //   Identify all different multimers

    for (i=0;i<recStack->nMultimers;i++)  {
      U1 = new Multimer ( 0 );
      U1->SetMultimer ( recStack->U[i],nInterfaces,interface,eType );
      U1->nUC = U1->mSize;
      k = -1;
      for (j=0;(j<nMultimers) && (k<0);j++)
        if (!U[j]->isComplemenatry(U1))  {
          if (U[j]->isIdentical(U1,D,rom,true))  k = j;
        }
      if (k>=0) {
        U[k]->nUC += U1->nUC;
        delete U1;
      } else if (nMultimers<nMonomers)  {
        U[nMultimers++] = U1;
      } else
        multimer_overflow = true;
    }

    for (i=0;i<nMultimers;i++)
      U[i]->nUC /= U[i]->mSize;

  }


  void  MultimerSet::CalcProperties ( PDomains      D,
                                      PPInterface   interface,
                                      int           nInterfaces,
                                      mmdb::ivector eType,
                                      int           nETypes,
                                      PSymNumber    SNum,
                                      mmdb::mat44 & rom )  {
  int i,j,k,nmm;

    if (nMultimers>1)  {

      for (i=0;i<nMultimers;i++)
        U[i]->type = -1;
      k   = 0;
      nmm = 0;
      for (i=0;i<nMultimers;i++)
        if (U[i]->type<0)  {
          U[i]->type = k;
          for (j=i+1;j<nMultimers;j++)
            if (U[j]->type<0)  {
              if (U[j]->isIdentical(U[i],D,rom,false))
                U[j]->type = k;
            }
          k++;
          if (U[i]->isMacroMolecule())  nmm++;
        }

  //    equiv_all = (k<=1);
      equiv_all = (nmm==1);

    } else  {

      U[0]->type = 0;
      equiv_all  = true;

    }

    maxStableSize = 1;
    stable_all    = true;
    for (i=0;i<nMultimers;i++)  {
      U[i]->CalcProperties ( D,interface,eType,nInterfaces,nETypes,
                             SNum,rom );
      if (U[i]->stable)  maxStableSize = U[i]->mSize;
                   else  stable_all    = false;
    }

  }


  int MultimerSet::getMaxMultimerSize()  {
  int i,k;
    k = 0;
    for (i=0;i<nMultimers;i++)
      if (U[i]->mSize>k)  k = U[i]->mSize;
    return k;
  }

}  // namespace pisa
