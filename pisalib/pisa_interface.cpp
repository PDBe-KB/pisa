// $Id: pisa_interface.cpp $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_interface <implemetation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Interface
//       ~~~~~~~~~  pisa::Interfaces
//
//  (C) E. Krissinel, 2007-2014
//
// =================================================================
//

#include <string.h>
#include <math.h>
#include<iostream>
using namespace std;

#include "mmdb2/mmdb_math_linalg.h"
#include "mmdb2/mmdb_tables.h"
#include "mmdb2/mmdb_utils.h"

#include "pisa_interface.h"
#include "pisa_types.h"
#include "pisa_defs.h"

#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end

namespace pisa  {

  // ====================  Calibrated parameters  ====================

  static mmdb::realtype  hb_en  = 0.444037;  // hydrogen bond factor
  // The following is because of overlooking, on early stage of
  // development, that HB is detected in about 50% of instances where SB
  // is formed. From version 2, HBs are removed from the list where
  // they are identical to SB, therefore the correction.
  static mmdb::realtype  sb_en  = hb_en/2.0+0.150028; // sald bridge factor
  static mmdb::realtype  ds_en  = 4.0;       // disulfide bond factor
  static mmdb::realtype  epsinv = 0.0/80.0;  // inverse dielectric
                                             //     susceptibility

  //  interfaces with interface area smaller than casualIntArea are
  //  considered as artifacts unless they make more than
  //  casualIntThresh fraction of the interfacing domains area
  static mmdb::realtype casualIntArea   = 40.0;
  static mmdb::realtype casualIntThresh = 0.03;
  static mmdb::realtype DSBondThresh    = 2.3;  // angstrom. This is an
               // average, which is not good everywhere. E.g.
               // a salt-bridge bound hexamer in 2cgy may be picked
               // only if 2.7 A S-S distance in disulphide bond is
               // allowed.
  static mmdb::realtype CovBondThresh   = 1.9;  // angstrom
  static mmdb::realtype OthBondThresh   = 5.0;  // angstrom

  void SetBondFactors ( mmdb::realtype HBEn, mmdb::realtype SBEn,
                        mmdb::realtype DSEn, mmdb::realtype epsf )  {
    hb_en  = HBEn;
    sb_en  = SBEn;
    ds_en  = DSEn;
    epsinv = epsf;
  }

  void GetBondFactors ( mmdb::realtype & HBEn, mmdb::realtype & SBEn,
                        mmdb::realtype & DSEn, mmdb::realtype & epsf )  {
    HBEn = hb_en;
    SBEn = sb_en;
    DSEn = ds_en;
    epsf = epsinv;
  }

  void writeBondFactors ( mmdb::io::RFile f )  {
  mmdb::byte Version=3;
    f.WriteByte ( &Version );
    f.WriteReal ( &hb_en   );
    f.WriteReal ( &sb_en   );
    f.WriteReal ( &ds_en   );
    f.WriteReal ( &epsinv  );
  }

  void readBondFactors ( mmdb::io::RFile f )  {
  mmdb::byte Version;
    f.ReadByte ( &Version );
    f.ReadReal ( &hb_en   );
    f.ReadReal ( &sb_en   );
    // The following is because of overlooking, on early stage of
    // development, that HB is detected in about 50% of instances where
    // SB is formed. From version 2, HBs are removed from the list where
    // they are identical to SB, therefore the correction.
    if (Version<3)
      sb_en += hb_en/2.0;
    if (Version>1)
      f.ReadReal ( &ds_en );
    f.ReadReal ( &epsinv  );

//    printf ( " hb_en=%6.2f  sb_en=%6.2f  ds_en=%6.2f\n",
//             hb_en,sb_en,ds_en );

  }

  void makeBondFactorsCIF (  mmdb::mmcif::PStruct mmCIFStruct )  {
    mmCIFStruct->PutReal ( hb_en,"hb_energy",10 );
    mmCIFStruct->PutReal ( sb_en,"sb_energy",10 );
    mmCIFStruct->PutReal ( ds_en,"ds_energy",10 );
    mmCIFStruct->PutReal ( epsinv, "eps_inv",10 );
  }

  void getBondFactorsCIF (  mmdb::mmcif::PStruct mmCIFStruct )  {
    mmCIFStruct->GetReal ( hb_en,"hb_energy",false );
    mmCIFStruct->GetReal ( sb_en,"sb_energy",false );
    mmCIFStruct->GetReal ( ds_en,"ds_energy",false );
    mmCIFStruct->GetReal ( epsinv, "eps_inv",false );
  }


  // ==========================  IBond  =============================

  void IBond::SetAtomPair ( ccp4srs::RAtomPair AP )  {
    serNum1 = AP.a1->serNum;
    serNum2 = AP.a2->serNum;
    res1    = AP.a1->GetResidueNo();
    res2    = AP.a2->GetResidueNo();
    dist    = sqrt(AP.a1->GetDist2(AP.a2));
  }

  void IBond::SetContact ( mmdb::RContact C )  {
    serNum1 = C.id1;
    serNum2 = C.id2;
    res1    = C.id1;
    res2    = C.id2;
    dist    = C.dist;
  }

  void IBond::Copy ( RIBond IB )  {
    serNum1 = IB.serNum1;
    serNum2 = IB.serNum2;
    res1    = IB.res1;
    res2    = IB.res2;
    dist    = IB.dist;
  }

  void IBond::write ( mmdb::io::RFile f )  {
    f.WriteInt  ( &serNum1 );
    f.WriteInt  ( &serNum2 );
    f.WriteInt  ( &res1    );
    f.WriteInt  ( &res2    );
    f.WriteReal ( &dist    );
  }

  void IBond::read ( mmdb::io::RFile f )  {
    f.ReadInt  ( &serNum1 );
    f.ReadInt  ( &serNum2 );
    f.ReadInt  ( &res1    );
    f.ReadInt  ( &res2    );
    f.ReadReal ( &dist    );
  }

  void SREffect::SetEffect ( mmdb::realtype eff, int resN )  {
    effect = eff;
    resNo  = resN;
  }

  void SREffect::Copy ( RSREffect sre )  {
    effect = sre.effect;
    resNo  = sre.resNo;
  }

  void SREffect::write ( mmdb::io::RFile f )  {
    f.WriteReal ( &effect );
    f.WriteInt  ( &resNo  );
  }

  void SREffect::read ( mmdb::io::RFile f )  {
    f.ReadReal ( &effect );
    f.ReadInt  ( &resNo  );
  }

  void MakeDFrame ( RDFrame dframe, PDomain D1, PDomain D2,
                    mmdb::mat44 & TMatrix )  {
  TFrame         frame2;
  mmdb::realtype dx,dy,dz,mx2,my2,mz2;
  int            i,j,k;

    if (D2->nAtoms>1)  {
      for (i=0;i<frameLen;i++)
        for (j=0;j<3;j++)  {
          frame2[i][j] = TMatrix[j][3];
          for (k=0;k<3;k++)
            frame2[i][j] += TMatrix[j][k]*D2->frame[i][k];
        }
    } else  {
      mx2 = D2->mx;
      my2 = D2->my;
      mz2 = D2->mz;
      mmdb::TransformXYZ ( TMatrix,mx2,my2,mz2 );
    }

    if ((D1->nAtoms>1) && (D2->nAtoms>1))  {
      // molecule against molecule
      for (i=0;i<frameLen;i++)
        for (j=0;j<frameLen;j++)  {
          dx = D1->frame[i][0] - frame2[j][0];
          dy = D1->frame[i][1] - frame2[j][1];
          dz = D1->frame[i][2] - frame2[j][2];
          dframe[i][j] = sqrt ( dx*dx + dy*dy + dz*dz );
        }
    } else if (D1->nAtoms==D2->nAtoms)  {
      // atom against atom
      dx = D1->mx - mx2;
      dy = D1->my - my2;
      dz = D1->mz - mz2;
      dx = sqrt ( dx*dx + dy*dy + dz*dz );
      for (i=0;i<frameLen;i++)
        for (j=0;j<frameLen;j++)
          dframe[i][j] = dx;
    } else if (D1->nAtoms==1)  {
      // atom against molecule
      for (j=0;j<frameLen;j++)  {
        dx = D1->mx - frame2[j][0];
        dy = D1->my - frame2[j][1];
        dz = D1->mz - frame2[j][2];
        dx = sqrt ( dx*dx + dy*dy + dz*dz );
        for (i=0;i<frameLen;i++)
          dframe[i][j] = dx;
      }
    } else  {
      // molecule against atom
      for (i=0;i<frameLen;i++)  {
        dx = D1->frame[i][0] - mx2;
        dy = D1->frame[i][1] - my2;
        dz = D1->frame[i][2] - mz2;
        dx = sqrt ( dx*dx + dy*dy + dz*dz );
        for (j=0;j<frameLen;j++)
          dframe[i][j] = dx;
      }
    }

  }

  void MakeDFrame ( RDFrame dframe, PDomain D1, mmdb::mat44 & TM1,
                                    PDomain D2, mmdb::mat44 & TM2 )  {
  TFrame   frame1;
  TFrame   frame2;
  mmdb::realtype mx1,my1,mz1,mx2,my2,mz2;
  mmdb::realtype dx,dy,dz;
  int      i,j,k;

    if (D1->nAtoms>1)  {
      for (i=0;i<frameLen;i++)
        for (j=0;j<3;j++)  {
          frame1[i][j] = TM1[j][3];
          for (k=0;k<3;k++)
            frame1[i][j] += TM1[j][k]*D1->frame[i][k];
        }
    } else  {
      mx1 = D1->mx;
      my1 = D1->my;
      mz1 = D1->mz;
      mmdb::TransformXYZ ( TM1,mx1,my1,mz1 );
    }

    if (D2->nAtoms>1)  {
      for (i=0;i<frameLen;i++)
        for (j=0;j<3;j++)  {
          frame2[i][j] = TM2[j][3];
          for (k=0;k<3;k++)
            frame2[i][j] += TM2[j][k]*D2->frame[i][k];
        }
    } else  {
      mx2 = D2->mx;
      my2 = D2->my;
      mz2 = D2->mz;
      mmdb::TransformXYZ ( TM2,mx2,my2,mz2 );
    }

    if ((D1->nAtoms>1) && (D2->nAtoms>1))  {
      for (i=0;i<frameLen;i++)
        for (j=0;j<frameLen;j++)  {
          dx = frame1[i][0] - frame2[j][0];
          dy = frame1[i][1] - frame2[j][1];
          dz = frame1[i][2] - frame2[j][2];
          dframe[i][j] = sqrt ( dx*dx + dy*dy + dz*dz );
        }
    } else if (D1->nAtoms==D2->nAtoms)  {
      // atom against atom
      dx = mx1 - mx2;
      dy = my1 - my2;
      dz = mz1 - mz2;
      dx = sqrt ( dx*dx + dy*dy + dz*dz );
      for (i=0;i<frameLen;i++)
        for (j=0;j<frameLen;j++)
          dframe[i][j] = dx;
    } else if (D1->nAtoms==1)  {
      // atom against molecule
      for (j=0;j<frameLen;j++)  {
        dx = mx1 - frame2[j][0];
        dy = my1 - frame2[j][1];
        dz = mz1 - frame2[j][2];
        dx = sqrt ( dx*dx + dy*dy + dz*dz );
        for (i=0;i<frameLen;i++)
          dframe[i][j] = dx;
      }
    } else  {
      // molecule against atom
      for (i=0;i<frameLen;i++)  {
        dx = frame1[i][0] - mx2;
        dy = frame1[i][1] - my2;
        dz = frame1[i][2] - mz2;
        dx = sqrt ( dx*dx + dy*dy + dz*dz );
        for (j=0;j<frameLen;j++)
          dframe[i][j] = dx;
      }
    }

  }


  // ========================  Interface  ===========================

  Interface::Interface() : mmdb::io::Stream()  {
    InitInterface();
  }

  Interface::Interface (  mmdb::io::RPStream Object )
            : mmdb::io::Stream ( Object )  {
    InitInterface();
  }

  Interface::~Interface()  {
    FreeMemory();
  }

  void Interface::set_asisKey_xml_int ( AS_IS_KEY as_is_Proc_int )  {
    asisKey = as_is_Proc_int;
  }
  
  void Interface::FreeMemory()  {
    if (symOp)  {
      delete[] symOp;
      symOp = NULL;
    }
    deleteHSDBonds ();
    deleteCovBonds ();
    deleteSREffects();
    mmdb::FreeVectorMemory ( bsa1   ,0 );
    mmdb::FreeVectorMemory ( bsa2   ,0 );
    mmdb::FreeVectorMemory ( solvEn1,0 );
    mmdb::FreeVectorMemory ( solvEn2,0 );
  }

  void Interface::InitInterface()  {
  int i,j;
    id            = 0;
    type          = 0;       // equivalence type 1,2,...
    domain1       = -1;
    domain2       = -1;
    dclass1       = DCLASS_None;
    dclass2       = DCLASS_None;
    symOpNo       = 0;
    rcsb_symop    = 0;       // rcsb symop serial number
    cell_i        = 0;
    cell_j        = 0;
    cell_k        = 0;
    nOcc          = 0;       // occurence number in ASU
    mmdb::Mat4Init ( TMatrix );
    for (i=0;i<frameLen;i++)  {
      for (j=0;j<frameLen;j++)
        dframe[i][j] = 0.0;
      dframe[i][i] = 1.0;
    }
    selHndInt1    = 0;
    selHndInt2    = 0;
    nIntAtoms1    = 0;
    nIntAtoms2    = 0;
    nRes1         = 0;
    nRes2         = 0;
    nIntRes1      = 0;
    nIntRes2      = 0;
    intArea1      = 0.0;
    intArea2      = 0.0;
    intArea       = 0.0;
    intDeltaG1    = 0.0;
    intDeltaG2    = 0.0;
    intDeltaG     = 0.0;
    aveDeltaG1    = 0.0;
    aveDeltaG2    = 0.0;
    aveDeltaG     = 0.0;
    probDeltaG1   = 0.0;
    probDeltaG2   = 0.0;
    probDeltaG    = 0.0;
    PValue1       = 0.0;
    PValue2       = 0.0;
    PValue        = 0.0;
    for (i=0;i<nASPs;i++)  {
      DeltaSAS1[i] = 0.0;
      DeltaSAS2[i] = 0.0;
    }
    stabEn        = 0.0;
    css           = 0.0;      // complexation significance score
    symOp         = NULL;
    HBond         = NULL;
    nHBonds       = 0;
    SBridge       = NULL;
    nSBridges     = 0;
    DSBond        = NULL;
    nDSBonds      = 0;
    CovBond       = NULL;
    nCovBonds     = 0;
    OthBond       = NULL;     //GDL:'other contacts' that are NOT h-bonds,s-bridges,cov bonds and disulfide bonds
    nOthBonds     = 0;        //GDL: No of 'other contacts' that are NOT h-bonds,s-bridges,cov bonds and disulfide bonds
    bsa1          = NULL;
    bsa2          = NULL;     // buried surface area
    solvEn1       = NULL;
    solvEn2       = NULL;     // solvation energy changes
    sre_stab1     = NULL;
    sre_stab2     = NULL;
    sre_destab1   = NULL;
    sre_destab2   = NULL;
    nStab1        = 0;
    nStab2        = 0;
    nDestab1      = 0;
    nDestab2      = 0;
    overlap       = false;    // flag "overlap suspicion"
    casual        = false;    // flag "casual interface"
    Xrel          = false;    // flag "crystallographically related"
    fixedLigand   = false;    // flag "fixed ligand interface"
  }

  void Interface::deleteHSDBonds()  {

    if (HBond)  {
      delete[] HBond;
      HBond = NULL;
    }
    nHBonds = 0;

    if (SBridge)  {
      delete[] SBridge;
      SBridge = NULL;
    }
    nSBridges = 0;

    if (DSBond)  {
      delete[] DSBond;
      DSBond = NULL;
    }
    nDSBonds = 0;

  }

  void Interface::deleteCovBonds()  {

    if (CovBond)  {
      delete[] CovBond;
      CovBond = NULL;
    }
    nCovBonds = 0;

  }
  //GDL: delete 'other contacts'
  void Interface::deleteOthBonds()  {

    if (OthBond)  {
      delete[] OthBond;
      OthBond = NULL;
    }
    nOthBonds = 0;

  }

  void Interface::deleteSREffects()  {

    if (sre_stab1)  {
      delete[] sre_stab1;
      sre_stab1 = NULL;
    }
    nStab1 = 0;

    if (sre_stab2)  {
      delete[] sre_stab2;
      sre_stab2 = NULL;
    }
    nStab2 = 0;

    if (sre_destab1)  {
      delete[] sre_destab1;
      sre_destab1 = NULL;
    }
    nDestab1 = 0;

    if (sre_destab2)  {
      delete[] sre_destab2;
      sre_destab2 = NULL;
    }
    nDestab2 = 0;

  }


  void Interface::Reset()  {
    FreeMemory   ();
    InitInterface();
  }

  void Interface::SetContact ( int   dom1, int  dom2, PDomains D,
                               int nSymOp, mmdb::mat44 & TM,
                               int  icell, int  jcell, int kcell )  {
    domain1    = dom1;
    domain2    = dom2;
    dclass1    = D->domain[dom1]->dclass;
    dclass2    = D->domain[dom2]->dclass;
    symOpNo    = nSymOp;
    rcsb_symop = nSymOp;   // just for now
    mmdb::Mat4Copy ( TM,TMatrix );
    MakeDFrame     ( dframe,D->domain[dom1],D->domain[dom2],TMatrix );
    cell_i     = icell;
    cell_j     = jcell;
    cell_k     = kcell;
    nOcc       = 1;
  }


  bool ASPfit = false;

  void SetASPFitSE ( bool fitMode )  {
    ASPfit = fitMode;
  }

  mmdb::realtype Interface::getStabEn() {
  mmdb::realtype DG;
  int            i;

    if (ASPfit)  {
      DG = 0.0;
      for (i=0;i<nASPs;i++)
        DG += ASP[i]*(DeltaSAS1[i]+DeltaSAS2[i]);
    } else
      DG = intDeltaG;

    return DG - hb_en*nHBonds - sb_en*nSBridges - ds_en*nDSBonds;

  }


  void  Interface::calcResBSA ( mmdb::PManager  MMDB,
                                int             selHndA,
                                mmdb::rvector & bsa,
                                mmdb::rvector & solvEn,
                                int           & nRes,
                                mmdb::rvector   intSAS,
                                mmdb::rvector   atomSE )  {
  mmdb::PPResidue res;
  mmdb::PPAtom    atom;
  mmdb::ivector   index;
  int             selHndR,nAtoms, i;

    selHndR = MMDB->NewSelection();
    MMDB->Select ( selHndR,mmdb::STYPE_RESIDUE,selHndA,mmdb::SKEY_NEW );
    MMDB->GetSelIndex ( selHndR,res,nRes );
    if (nRes>0)  {
      mmdb::GetVectorMemory ( bsa   ,nRes,0 );
      mmdb::GetVectorMemory ( solvEn,nRes,0 );
      if (intSAS && atomSE)  {
        mmdb::GetVectorMemory ( index ,nRes,0 );
        for (i=0;i<nRes;i++)  {
          bsa   [i] = 0.0;
          solvEn[i] = 0.0;
          index [i] = res[i]->index;
          res   [i]->index = i;
        }
        MMDB->GetSelIndex ( selHndA,atom,nAtoms );
        for (i=0;i<nAtoms;i++)  {
          bsa   [atom[i]->GetResidue()->index] += intSAS[i];
          solvEn[atom[i]->GetResidue()->index] += atomSE[i];
        }
        for (i=0;i<nRes;i++)
          res[i]->index = index[i];
        mmdb::FreeVectorMemory ( index,0 );
      } else  {
        for (i=0;i<nRes;i++)  {
          bsa   [i] = 0.0;
          solvEn[i] = 0.0;
        }
      }
    }
    MMDB->DeleteSelection ( selHndR );

  }

  DefineClass(SortSRE);

  class SortSRE : public mmdb::QuickSort  {

    public :

      SortSRE () : QuickSort()  {}
      ~SortSRE() {}

      virtual int  Compare ( int i, int j )  {
      // sort by decreasing data[i]
        if (((mmdb::rvector)data)[i]<((mmdb::rvector)data)[j])  return  1;
        if (((mmdb::rvector)data)[i]>((mmdb::rvector)data)[j])  return -1;
        return 0;
      }

      virtual void Swap ( int i, int j )  {
      mmdb::realtype sre = ((mmdb::rvector)data)[i];
      int            rn  = resN[i];
        ((mmdb::rvector)data)[i] = ((mmdb::rvector)data)[j];
        ((mmdb::rvector)data)[j] = sre;
        resN[i] = resN[j];
        resN[j] = rn;
      }

      void Sort ( mmdb::rvector sre, mmdb::ivector resNo, int nRes )  {
        resN = resNo;
        QuickSort::Sort ( sre,nRes );
      }

    protected :
      mmdb::ivector resN;

  };


  void Interface::calcSREffs ( mmdb::rvector solvEn,
                               mmdb::rvector asa,
                               mmdb::rvector bsa,
                               int           nRes,
                               int           domNo,
                               RPSREffect    sre_stab,
                               int         & nStab,
                               RPSREffect    sre_destab,
                               int         & nDestab )  {
  // Calculation of maximal and minimal single-residue effects.
  // Chemical bonds must be calculated prior calling this function.
  SortSRE        sort;
  mmdb::rvector  sreff;
  mmdb::ivector  resNo;
  mmdb::realtype fsa;
  int            i,j;

    mmdb::GetVectorMemory ( sreff,nRes,0 );
    mmdb::GetVectorMemory ( resNo,nRes,0 );

    for (i=0;i<nRes;i++)  {

      sreff[i] = solvEn[i];

      for (j=0;j<nHBonds;j++)
        if (((domNo<=1) && (HBond[j].res1==i)) ||
            ((domNo>=2) && (HBond[j].res2==i)))
          sreff[i] += hb_en;
      for (j=0;j<nSBridges;j++)
        if (((domNo<=1) && (SBridge[j].res1==i)) ||
            ((domNo>=2) && (SBridge[j].res2==i)))
          sreff[i] += sb_en;
      for (j=0;j<nDSBonds;j++)
        if (((domNo<=1) && (DSBond[j].res1==i)) ||
            ((domNo>=2) && (DSBond[j].res2==i)))
          sreff[i] += ds_en;

      resNo[i] = i;

    }

    sort.Sort ( sreff,resNo,nRes );

    nStab = 0;
    while ((nStab<nRes) && (nStab<6) && (sreff[0]>0.0) &&
           (sreff[0]-sreff[nStab]<=0.5))
      nStab++;
    if (nStab>0)  {
      sre_stab = new SREffect[nStab];
      for (i=0;i<nStab;i++)
        sre_stab[i].SetEffect ( sreff[i],resNo[i] );
    }

    nDestab = 0;
    i       = nRes-1;
    j       = nRes-1;
    fsa     = bsa[resNo[j]]/asa[resNo[j]]/2.5;
    while ((i>=0) && (nDestab<6) && (sreff[i]<0.0) &&
           (bsa[resNo[i]]/asa[resNo[i]]>fsa) &&
           (sreff[i]-sreff[j]<=0.5))  {
      nDestab++;
      i--;
    }
    if (nDestab>0)  {
      sre_destab = new SREffect[nDestab];
      for (i=0;i<nDestab;i++)  {
        sre_destab[i].SetEffect ( sreff[j],resNo[j] );
        j--;
      }
    }

    mmdb::FreeVectorMemory ( resNo,0 );
    mmdb::FreeVectorMemory ( sreff,0 );

  }


  PROSURF_RC Interface::calcInterface ( mmdb::PManager    MMDB,
                                        PDomain           D1,
                                        PDomain           D2,
                                        PProSurf          proSurf,
                                        PMolRefIndex      molRef,
                                        ccp4srs::PManager SRS )  {
  //   The function assumes that there are two identical models
  // in the structure. If second model is not found, it is
  // generated automatically and is not deleted on return.
  mmdb::SymOp   SymOp;
  SolvEnergy    SolvEnergy;
  mmdb::PModel  model;
  mmdb::PPAtom  atom1,atom2;
  mmdb::rvector x1,y1,z1, x2,y2,z2, SAS1,intSAS1,SAS2,intSAS2, atomSE;
  mmdb::mat44   ST,tm;
  char          symOpTitle[300];
  int           nAtoms1,nAtoms2;
  int           selHnd1,selHnd2,selHndSurf2,selHndR;
  int           i;//,mrc;
  PROSURF_RC    rc;


    deleteHSDBonds  ();
    deleteCovBonds  ();
    deleteSREffects ();
    mmdb::FreeVectorMemory ( bsa1   ,0 );
    mmdb::FreeVectorMemory ( bsa2   ,0 );
    mmdb::FreeVectorMemory ( solvEn1,0 );
    mmdb::FreeVectorMemory ( solvEn2,0 );
    nRes1 = 0;
    nRes2 = 0;

    

    //  1. Check that molecule has been complemented with hydrogens,
    //     which are needed for electrostatics calculations.
/* --- candidate for removal
    uddHnd = MMDB->GetUDDHandle ( mmdb::UDR_HIERARCHY,hydrogen_udd );
    if (uddHnd>0)
      MMDB->GetUDData ( uddHnd,mrc );
    else  {
      uddHnd = MMDB->RegisterUDInteger ( mmdb::UDR_HIERARCHY,hydrogen_udd );
      mrc = 0;
    }

    if (!mrc)  {
      // hydrogens were not added, add them and store a flag in MMDB
      SRS->addHydrogens     ( MMDB,NULL );
      MMDB->FinishStructEdit();
      // serial numbers should never be touched here as they are
      // used for carrying the atom selections between the CPUs
      MMDB->PDBCleanup      ( mmdb::PDBCLEAN_INDEX );
      MMDB->PutUDData       ( uddHnd,1  );
      MMDB->GetUDData       ( uddHnd,mrc );
    }
*/
    //  2. Check that 2nd model is there and generate it if it is not.

    if (!MMDB->GetModel(2))  {
      model = mmdb::newModel();
      model->Copy ( MMDB->GetModel(1) );
      MMDB->AddModel ( model );
    }


    //  3. Select first domain in 1st model.

    selHnd1 = D1->SelectDomain ( MMDB,mmdb::STYPE_ATOM,1,false );
    MMDB->Select ( selHnd1,mmdb::STYPE_ATOM,0,"*",
                           mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                           "*","*","H","*",mmdb::SKEY_CLR );

    //    and set up its surface selection for faster calculations
    proSurf->selHndSurf1 = D1->selHndSurf;


    //  4. In anticipation of coordinate transformations, select 2nd
    //     domain in 2nd model

    selHnd2 = D2->SelectDomain ( MMDB,mmdb::STYPE_ATOM,2,false );
    MMDB->Select ( selHnd2,mmdb::STYPE_ATOM,0,"*",
                           mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                           "*","*","H","*",mmdb::SKEY_CLR );

    //    now make a new selection of its surface in 2nd model
    MMDB->GetSelIndex ( D2->selHndSurf,atom2,nAtoms2 );
    if (nAtoms2>0)  {
      selHndSurf2 = MMDB->NewSelection();
      MMDB->SelectNeighbours ( selHndSurf2,mmdb::STYPE_ATOM,
                               atom2,nAtoms2,0.0,0.001,mmdb::SKEY_NEW );
      proSurf->selHndSurf2 = selHndSurf2;
    } else  {
      //  We do not have surfaces selected. Ok, then calculations
      // will take somewhat longer.
      proSurf->selHndSurf1 = 0;
      proSurf->selHndSurf2 = 0;
      selHndSurf2 = 0;
    }

    // 5. Apply symmetry operations to domains

        MMDB->GetSelIndex ( selHnd1,atom1,nAtoms1 );
    if (D1->ncsOpNo>0)  {
      mmdb::GetVectorMemory ( x1,nAtoms1,0 );
      mmdb::GetVectorMemory ( y1,nAtoms1,0 );
      mmdb::GetVectorMemory ( z1,nAtoms1,0 );
      for (i=0;i<nAtoms1;i++)
        if (atom1[i])  {
          x1[i] = atom1[i]->x;
          y1[i] = atom1[i]->y;
          z1[i] = atom1[i]->z;
          atom1[i]->Transform ( D1->ncs_m );
        }
    }

    MMDB->GetSelIndex ( selHnd2,atom2,nAtoms2 );
    mmdb::GetVectorMemory ( x2,nAtoms2,0 );
    mmdb::GetVectorMemory ( y2,nAtoms2,0 );
    mmdb::GetVectorMemory ( z2,nAtoms2,0 );
    if (D2->ncsOpNo>0)  mmdb::Mat4Mult ( tm,TMatrix,D2->ncs_m );
                  else  mmdb::Mat4Copy ( TMatrix,tm );
    for (i=0;i<nAtoms2;i++)
      if (atom2[i])  {
        x2[i] = atom2[i]->x;
        y2[i] = atom2[i]->y;
        z2[i] = atom2[i]->z;
        atom2[i]->Transform ( tm );
      }

    
    // 6. Calculate the interface

    mmdb::GetVectorMemory ( SAS1   ,nAtoms1,0 );
    mmdb::GetVectorMemory ( intSAS1,nAtoms1,0 );
    mmdb::GetVectorMemory ( SAS2   ,nAtoms2,0 );
    mmdb::GetVectorMemory ( intSAS2,nAtoms2,0 );
    rc = proSurf->calcInterface ( atom1,nAtoms1,atom2,nAtoms2,
                                  molRef,SURFF_Selected,
                                  SAS1,intSAS1,
                                  SAS2,intSAS2 );
    //Grisell test print interface residues
    /*
    if (rc==PROSURF_Ok)  {
    for (i=0;i<nAtoms1;i++)
      {
	cout << atom1[i]->GetResName()<<"\t"<<atom1[i]->GetSeqNum()<<"\t"<<atom1[i]->element<<"\n";
      }
    for (i=0;i<nAtoms2;i++)
      {
        cout << atom2[i]->GetResName()<<"\t"<<atom2[i]->GetSeqNum()<<"\t"<<atom2[i]->element<<"\n";
      }
      }*/
    
    //cout << atom1->GetResName()<<"\t"<<atom1->GetSeqNum ()<<"\t"<<atom1->element<<"\t"<< atom2->GetResName() <<"\t"<<atom2->GetSeqNum()<<"\t"<<atom2->element<< "\n";
    
    // 7. Get the interface selections and parameters

    if (rc==PROSURF_Ok)  {

      selHndInt1 = proSurf->selHndInt1;
      selHndInt2 = proSurf->selHndInt2;

      nIntAtoms1 = MMDB->GetSelLength ( selHndInt1 );
      nIntAtoms2 = MMDB->GetSelLength ( selHndInt2 );
      
      checkOverlap ( MMDB,atom1,nAtoms1,atom2,nAtoms2 );
      
      selHndR    = MMDB->NewSelection();
      MMDB->Select ( selHndR,mmdb::STYPE_RESIDUE,selHndInt1,mmdb::SKEY_NEW );
      nIntRes1   = MMDB->GetSelLength ( selHndR );
      MMDB->Select ( selHndR,mmdb::STYPE_RESIDUE,selHndInt2,mmdb::SKEY_NEW );
      nIntRes2   = MMDB->GetSelLength ( selHndR );
      MMDB->DeleteSelection ( selHndR );

      intArea1   = proSurf->intArea1;
      intArea2   = proSurf->intArea2;
      intArea    = (intArea1+intArea2)/2.0;

      if (intArea>casualIntArea)
            casual = false;
      else  casual = ((intArea1<=casualIntThresh*D1->surfArea) &&
                      (intArea2<=casualIntThresh*D2->surfArea));

    } else  {

      selHndInt1 = 0;
      selHndInt2 = 0;

      nIntAtoms1 = 0;
      nIntAtoms2 = 0;
      nIntRes1   = 0;
      nIntRes2   = 0;

      overlap    = false;
      casual     = true;

      intArea1   = 0.0;
      intArea2   = 0.0;
      intArea    = 0.0;

    }

    if (symOpNo>=0)  {
      SymOp.SetSymOp   ( MMDB->GetSymOp(symOpNo) );
      SymOp.GetTMatrix ( ST );
    } else
      mmdb::Mat4Init ( ST );
    ST[0][3] += cell_i;
    ST[1][3] += cell_j;
    ST[2][3] += cell_k;
    SymOp.CompileOpTitle ( symOpTitle,ST,false );
    mmdb::CreateCopy ( symOp,symOpTitle );

    MMDB->DeleteSelection ( selHndSurf2 );


    // 8. Calculate interface properties

    if (intArea>0.0)  {

      mmdb::GetVectorMemory ( atomSE,mmdb::IMax(nAtoms1,nAtoms2),0 );

      intDeltaG1 = SolvEnergy.calcIntSolvEnergy ( atom1,SAS1,intSAS1,
                                                  nAtoms1,atomSE,
                                                  molRef );
      for (i=0;i<nASPs;i++)
        DeltaSAS1[i] = SolvEnergy.aspArea[i];

      calcResBSA ( MMDB,selHnd1,bsa1,solvEn1,nRes1,intSAS1,atomSE );

      intDeltaG2 = SolvEnergy.calcIntSolvEnergy ( atom2,SAS2,intSAS2,
                                                  nAtoms2,atomSE,
                                                  molRef );
      for (i=0;i<nASPs;i++)
        DeltaSAS2[i] = SolvEnergy.aspArea[i];

      calcResBSA ( MMDB,selHnd2,bsa2,solvEn2,nRes2,intSAS2,atomSE );

      mmdb::FreeVectorMemory ( atomSE,0 );

      intDeltaG    = intDeltaG1 + intDeltaG2;

      PValue1      = SolvEnergy.calcIntPValue ( atom1,SAS1,nAtoms1,
                                                nIntAtoms1,intDeltaG1,
                                                molRef );
      aveDeltaG1   = SolvEnergy.aveDeltaG;
      probDeltaG1  = SolvEnergy.probDeltaG;

      PValue2      = SolvEnergy.calcIntPValue ( atom2,SAS2,nAtoms2,
                                                nIntAtoms2,intDeltaG2,
                                                molRef );
      aveDeltaG2   = SolvEnergy.aveDeltaG;
      probDeltaG2  = SolvEnergy.probDeltaG;

      if ((PValue1>=0.0) && (PValue2>=0.0))
            PValue = sqrt(PValue1*PValue2);
      else  PValue = -1.0;

      aveDeltaG  = aveDeltaG1  + aveDeltaG2;
      probDeltaG = probDeltaG1 + probDeltaG2;

      calcChemBonds ( MMDB,SRS);

      //GDL: calculate other contacts - start
      
      CalcOtherContacts (MMDB,atom1,nAtoms1,atom2,nAtoms2);
      
      //calculate othercontacts - end 
      
      calcSREffs ( solvEn1,D1->asa,bsa1,nRes1,1,
                   sre_stab1,nStab1,sre_destab1,nDestab1 );
      calcSREffs ( solvEn2,D2->asa,bsa2,nRes2,2,
                   sre_stab2,nStab2,sre_destab2,nDestab2 );

      // stabilization energy
      stabEn = getStabEn(); // intDeltaG - 0.8*nHBonds - 0.4*nSBridges;

      css = 0.0;  // to be defined from assemblies

    } else  {

      calcResBSA ( MMDB,selHnd1,bsa1,solvEn1,nRes1,NULL,NULL );
      calcResBSA ( MMDB,selHnd2,bsa2,solvEn2,nRes2,NULL,NULL );

      intDeltaG1  = 0.0;
      intDeltaG2  = 0.0;
      intDeltaG   = 0.0;

      aveDeltaG1  = 0.0;
      aveDeltaG2  = 0.0;
      aveDeltaG   = 0.0;

      probDeltaG1 = 0.0;
      probDeltaG2 = 0.0;
      probDeltaG  = 0.0;

      PValue1     = 0.0;
      PValue2     = 0.0;
      PValue      = 0.0;

    }

     mmdb::FreeVectorMemory ( intSAS2,0 );
     mmdb::FreeVectorMemory ( SAS2   ,0 );
     mmdb::FreeVectorMemory ( intSAS1,0 );
     mmdb::FreeVectorMemory ( SAS1   ,0 );


    // 9. Restore domains

    if (D1->ncsOpNo>0)  {
      for (i=0;i<nAtoms1;i++)
        if (atom1[i])  {
          atom1[i]->x = x1[i];
          atom1[i]->y = y1[i];
          atom1[i]->z = z1[i];
        }
       mmdb::FreeVectorMemory ( x1,0 );
       mmdb::FreeVectorMemory ( y1,0 );
       mmdb::FreeVectorMemory ( z1,0 );
    }

    for (i=0;i<nAtoms2;i++)
      if (atom2[i])  {
        atom2[i]->x = x2[i];
        atom2[i]->y = y2[i];
        atom2[i]->z = z2[i];
      }
     mmdb::FreeVectorMemory ( x2,0 );
     mmdb::FreeVectorMemory ( y2,0 );
     mmdb::FreeVectorMemory ( z2,0 );

    MMDB->DeleteSelection ( selHnd1  );


    // 10. Reflect selection of the interface in 2nd structure
    //     onto 1st model

    if ((nIntAtoms2>0) || (selHndInt2>0))  {
      MMDB->GetSelIndex ( selHndInt2,atom2,nAtoms2 );
      MMDB->SelectNeighbours ( selHnd2,mmdb::STYPE_ATOM,
                               atom2,nAtoms2,0.0,0.001,mmdb::SKEY_NEW );
      MMDB->DeleteSelection ( selHndInt2 );
      selHndInt2 = selHnd2;
    } else
      MMDB->DeleteSelection ( selHnd2 );

    return rc;

  }

  void Interface::CalcOtherContacts ( mmdb::PManager MMDB,mmdb::PPAtom atom1, int nat1,
				       mmdb::PPAtom atom2, int nat2 )  {
    mmdb::PContact   cont;
    mmdb::PAtom      a1,a2;
    int              i,j,ncont, sec1, sec2, s1, s2;
    mmdb::ivector    SBsec1,SBsec2,HBsec1,HBsec2,CovBsec1,CovBsec2,DSsec1,DSsec2;
    bool             isBond;           
    mmdb::ivector    bondsec1,bondsec2;
    mmdb::PContact   othcont;
    mmdb::PAtomName  atnames1,atnames2;
    mmdb::AtomName   at1name,at2name;
    int              nothcont;                                                                                                                                    
    mmdb::PAtom      at1,at2;

    //delete contacts in array OthBond and reset number of other contacts to nOthBonds = 0 
    deleteOthBonds();

    mmdb::GetVectorMemory ( SBsec1  ,nSBridges,0 );
    mmdb::GetVectorMemory ( SBsec2  ,nSBridges,0 );
    mmdb::GetVectorMemory ( HBsec1  ,nHBonds  ,0 );
    mmdb::GetVectorMemory ( HBsec2  ,nHBonds  ,0 );
    mmdb::GetVectorMemory ( CovBsec1,nCovBonds,0 );
    mmdb::GetVectorMemory ( CovBsec2,nCovBonds,0 );
    mmdb::GetVectorMemory ( DSsec1  ,nDSBonds ,0 );
    mmdb::GetVectorMemory ( DSsec2  ,nDSBonds ,0 );

    ncont = nSBridges + nHBonds + nCovBonds + nDSBonds;
    mmdb::GetVectorMemory ( bondsec1,ncont,0 );
    mmdb::GetVectorMemory ( bondsec2,ncont,0 );

    
    atnames1 = new mmdb::AtomName[ncont];
    atnames2 = new mmdb::AtomName[ncont];

    deleteOthBonds();

    othcont = NULL;
    nothcont = 0 ;

    // Seek contacts within distance threshold OthBondThresh
    MMDB->SeekContacts(atom1,nat1,atom2,nat2,0.0,OthBondThresh,0,othcont,nothcont,0,NULL,0,0);

    // Get list of atoms from interface residues 1 and 2 that form salt-bridges        
    for (i=0;i<nSBridges;i++)  {
      a1 = GetAtom ( MMDB,SBridge[i].serNum1 );
      a2 = GetAtom ( MMDB,SBridge[i].serNum2 );
      SBsec1[i]   = a1->GetSeqNum ();
      SBsec2[i]   = a2->GetSeqNum ();
      bondsec1[i] = a1->GetSeqNum ();
      bondsec2[i] = a2->GetSeqNum ();
      strcpy ( atnames1[i],a1->name );
      strcpy ( atnames2[i],a2->name );
   
    }
    
    // Get list of atoms from interface residues 1 and 2 that form H-bonds
    for (i=0;i<nHBonds;i++)  {
      a1 = GetAtom ( MMDB,HBond[i].serNum1 );
      a2 = GetAtom ( MMDB,HBond[i].serNum2 );
      HBsec1[i] = a1->GetSeqNum ();
      HBsec2[i] = a2->GetSeqNum ();
      bondsec1[nSBridges+i]=a1->GetSeqNum ();
      bondsec2[nSBridges+i]=a2->GetSeqNum ();
      strcpy ( atnames1[nSBridges+i],a1->name );
      strcpy ( atnames2[nSBridges+i],a2->name );

    }
    // Get list of atoms from interface residues 1 and 2 that form covalent bonds  
    for (i=0;i<nCovBonds;i++)  {
      a1 = GetAtom ( MMDB,CovBond[i].serNum1 );
      a2 = GetAtom ( MMDB,CovBond[i].serNum2 );
      CovBsec1[i] = a1->GetSeqNum ();
      CovBsec2[i] = a2->GetSeqNum ();
      bondsec1[nSBridges+nHBonds+i]=a1->GetSeqNum ();
      bondsec2[nSBridges+nHBonds+i]=a2->GetSeqNum ();
      strcpy ( atnames1[nSBridges+nHBonds+i],a1->name );
      strcpy ( atnames2[nSBridges+nHBonds+i],a2->name );
      
    }

    // Get list of atoms from interface residues 1 and 2 that form disulfide bonds
    for (i=0;i<nDSBonds;i++)  {
      a1 = GetAtom ( MMDB,DSBond[i].serNum1 );
      a2 = GetAtom ( MMDB,DSBond[i].serNum2 );
      DSsec1[i] = a1->GetSeqNum ();
      DSsec2[i] = a2->GetSeqNum ();
      bondsec1[nSBridges+nHBonds+nCovBonds+i] = a1->GetSeqNum ();
      bondsec2[nSBridges+nHBonds+nCovBonds+i] = a2->GetSeqNum ();
      strcpy ( atnames1[nSBridges+nHBonds+nCovBonds+i],a1->name );
      strcpy ( atnames2[nSBridges+nHBonds+nCovBonds+i],a2->name );
      
    }
      
    nOthBonds = 0;
    OthBond   = new IBond[nothcont];

    // Seek for contacts that are not already labeled as h-bonds, s-bridges, cov-bonds or disulfide bonds, and save them in the list OthBond and counter nOthBonds
    for (i=0;i<nothcont;i++)  {

      isBond = false;
	  
      if (othcont[i].dist<OthBondThresh)  {
        at1     = atom1[othcont[i].id1] ;
        at2     = atom2[othcont[i].id2] ;
        sec1    = at1->GetSeqNum ();
        sec2    = at2->GetSeqNum ();
        strcpy ( at1name,at1->name );
        strcpy ( at2name,at2->name );

	      for (j=0;j<ncont;j++) {
    		  if (sec1 == bondsec1[j] && sec2 == bondsec2[j] && 
		      (strcmp(at1name,atnames1[j])==0) && (strcmp(at2name,atnames2[j])==0))  {
		        isBond = true ; 
		      }
		    }
	      
	      if(!isBond)  {		
		OthBond[nOthBonds].serNum1 = at1->serNum;
		OthBond[nOthBonds].serNum2 = at2->serNum;
		OthBond[nOthBonds].dist    = othcont[i].dist;
		nOthBonds++;  
    		}
	      
	    }
	  
	  }
    
    mmdb::FreeVectorMemory ( SBsec1  ,0 );
    mmdb::FreeVectorMemory ( SBsec2  ,0 );
    mmdb::FreeVectorMemory ( HBsec1  ,0 );
    mmdb::FreeVectorMemory ( HBsec2  ,0 );
    mmdb::FreeVectorMemory ( CovBsec1,0 );
    mmdb::FreeVectorMemory ( CovBsec2,0 );
    mmdb::FreeVectorMemory ( DSsec1  ,0 );
    mmdb::FreeVectorMemory ( DSsec2  ,0 );

    mmdb::FreeVectorMemory ( bondsec1,0 );
    mmdb::FreeVectorMemory ( bondsec2,0 );

    delete[] atnames1;
    delete[] atnames2;

  }


  void Interface::checkOverlap ( mmdb::PManager MMDB,
                                 mmdb::PPAtom atom1, int nat1,
                                 mmdb::PPAtom atom2, int nat2 )  {

  mmdb::PContact cont;
  PIBond         cov;
  mmdb::PAtom    a1,a2;
  int            type1[4],type2[4];
  int            i,j,ncont,maxncont,el1,el2,nCovAlloc;
  bool           B1,B2;
  
    deleteCovBonds();

    if ((nat1<=1) && (nat2<=1))
      overlap = false;
    else  {

      cont  = NULL;
      ncont = 0;

      MMDB->SeekContacts ( atom1,nat1,atom2,nat2,0.0,2.4,0,
                           cont,ncont,0,NULL,0,0 );
      
      maxncont = mmdb::mround ( 1.1*mmdb::IMax(nIntAtoms1,nIntAtoms2) );

      overlap = (ncont>maxncont);

      if (!overlap)  {

        type1[0] = mmdb::getElementNo ( "BE" );
        type2[0] = mmdb::getElementNo ( "F"  );
        type1[1] = mmdb::getElementNo ( "MG" );
        type2[1] = mmdb::getElementNo ( "CL" );
        type1[2] = mmdb::getElementNo ( "CA" );
        type2[2] = mmdb::getElementNo ( "BR" );
        type1[3] = mmdb::getElementNo ( "SR" );
        type2[3] = mmdb::getElementNo ( "I" );

        nCovAlloc = 0;
        for (i=0;i<ncont;i++)
          if (cont[i].dist<CovBondThresh)  {
            // assume a covalent link if atoms are right
            if ((dclass1==DCLASS_Ligand) || (dclass2==DCLASS_Ligand))  {
              a1  = atom1[cont[i].id1];
              a2  = atom2[cont[i].id2];
              el1 = mmdb::getElementNo ( a1->element );
              el2 = mmdb::getElementNo ( a2->element );
              B1  = false;
              B2  = false;
              for (j=0;j<4;j++)  {
                if (!B1)  B1 = (el1>type1[j]) && (el1<type2[j]);
                if (!B2)  B2 = (el2>type1[j]) && (el2<type2[j]);
              }
              if (B1 && B2)  {
                if (nCovBonds<=nCovAlloc)  {
                  nCovAlloc = nCovBonds + 10;
                  cov = new IBond[nCovAlloc];
                  for (j=0;j<nCovBonds;j++)
                    cov[j].Copy ( CovBond[j] );
                  if (CovBond)  delete[] CovBond;
                  CovBond = cov;
                }
                CovBond[nCovBonds].serNum1 = a1->serNum;
                CovBond[nCovBonds].serNum2 = a2->serNum;
                CovBond[nCovBonds].dist    = cont[i].dist;
                nCovBonds++;
              }
            }
          }

      }

      if (cont)  delete[] cont;

    }

  }

  void  Interface::calcChemBonds ( mmdb::PManager    MMDB,
                                  ccp4srs::PManager SRS )  {
  ccp4srs::PAtomPair HB,SB;
  mmdb::PPResidue    Res1,Res2;
  mmdb::PPAtom       S1,S2;
  mmdb::realtype     dist;
  int                selRes1,selRes2, selHnd1,selHnd2, selSulph1,selSulph2;
  int                nIRes1,nIRes2, nS1,nS2;
  int                i,j,k,r;

    deleteHSDBonds();

    selRes1 = MMDB->NewSelection();
    MMDB->Select      ( selRes1   ,mmdb::STYPE_RESIDUE,
                        selHndInt1,mmdb::SKEY_NEW );
    MMDB->GetSelIndex ( selRes1,Res1,nIRes1 );
    selRes2 = MMDB->NewSelection();
    MMDB->Select      ( selRes2   ,mmdb::STYPE_RESIDUE,
                        selHndInt2,mmdb::SKEY_NEW );
    MMDB->GetSelIndex ( selRes2,Res2,nIRes2 );

    HB = NULL;
    SB = NULL;
    
    k  = SRS->CalcHBonds ( Res1,nIRes1,Res2,nIRes2,HB,nHBonds,
                           SB,nSBridges,NULL,NULL,true );
        
    if (k>=0)  {

      if (nHBonds>0)  {
        HBond = new IBond[nHBonds];
        for (i=0;i<nHBonds;i++)
          HBond[i].SetAtomPair ( HB[i] );
      }
      if (nSBridges>0)  {
        SBridge = new IBond[nSBridges];
        for (i=0;i<nSBridges;i++)
          SBridge[i].SetAtomPair ( SB[i] );
      }

    } else  {
      nHBonds   = 0;
      nSBridges = 0;
    }
    

    if (HB)  delete[] HB;
    if (SB)  delete[] SB;


    selHnd1   = MMDB->NewSelection();
    selHnd2   = MMDB->NewSelection();
    selSulph1 = MMDB->NewSelection();
    selSulph2 = MMDB->NewSelection();

    MMDB->Select ( selHnd1,mmdb::STYPE_RESIDUE,0,
                           "*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                           "CYS","*","*","*",mmdb::SKEY_AND );
    MMDB->Select ( selSulph1,mmdb::STYPE_ATOM,selHnd1,mmdb::SKEY_NEW );
    MMDB->Select ( selSulph1,mmdb::STYPE_ATOM,0,
                             "*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                             "*","*","S","*",mmdb::SKEY_AND );
    MMDB->GetSelIndex ( selSulph1,S1,nS1 );

    if (nS1>0)  {
      MMDB->Select ( selHnd2,mmdb::STYPE_RESIDUE,0,
                             "*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                             "CYS","*","*","*",mmdb::SKEY_AND );
      MMDB->Select ( selSulph2,mmdb::STYPE_ATOM,selHnd2,mmdb::SKEY_NEW );
      MMDB->Select ( selSulph2,mmdb::STYPE_ATOM,0,
                               "*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                               "*","*","S","*",mmdb::SKEY_AND );
      MMDB->GetSelIndex ( selSulph2,S2,nS2 );
      if (nS2>0)  {
        // recon it's just a few sulphurs in the interface, so do not
        // mess with the SeekContact
        nDSBonds = 0;
        DSBond   = new IBond[nS1*nS2];
        for (i=0;i<nS1;i++)
          for (j=0;j<nS2;j++)  {
            dist = sqrt ( S1[i]->GetDist2(S2[j]) );
            if (dist<=DSBondThresh)  {
              DSBond[nDSBonds].serNum1 = S1[i]->serNum;
              DSBond[nDSBonds].serNum2 = S2[j]->serNum;
              DSBond[nDSBonds].dist    = dist;
              nDSBonds++;
            }
          }
      }
    }

    MMDB->DeleteSelection ( selSulph1 );
    MMDB->DeleteSelection ( selSulph2 );
    MMDB->DeleteSelection ( selHnd1   );
    MMDB->DeleteSelection ( selHnd2   );
    MMDB->DeleteSelection ( selRes1   );
    MMDB->DeleteSelection ( selRes2   );

    
  }



  #define rot_threshold    0.001
  #define trans_threshold  1.0

  bool Interface::checkTMatrix ( mmdb::mat44  & TM )  {
  // Returns true if TMatrix equals to TM
  mmdb::realtype dT,dR, d;
  int      i;

    dT = 0.0;
    dR = 0.0;

    for (i=0;i<3;i++)  {
      d   = TMatrix[i][3] - TM[i][3];
      dT += d*d;
      d   = TMatrix[i][0] - TM[i][0];
      dR += d*d;
      d   = TMatrix[i][1] - TM[i][1];
      dR += d*d;
      d   = TMatrix[i][2] - TM[i][2];
      dR += d*d;
    }

    return ((dR<rot_threshold) && (dT<trans_threshold));

  }


  bool Interface::checkDFrame ( RDFrame df, bool direct )  {
  //  Returns true if frames compare
  int      i,j;
  bool  ok;

    ok = true;
    if (direct)  {
      for (i=0;(i<frameLen) && ok;i++)
        for (j=0;(j<frameLen) && ok;j++)
          ok = fabs(df[i][j]-dframe[i][j])<2.25;
    } else  {
      for (i=0;(i<frameLen) && (ok);i++)
        for (j=0;(j<frameLen) && (ok);j++)
          ok = fabs(df[j][i]-dframe[i][j])<2.25;
    }

    return ok;

  }

  void Interface::writeMaxFrameDiff ( mmdb::mat44 & df )  {
  mmdb::realtype d1,d2;
  int      i,j;

    d1 = 0.0;
    d2 = 0.0;
    for (i=0;i<frameLen;i++)
      for (j=0;j<frameLen;j++)
        d1 = mmdb::RMax ( d1,fabs(df[i][j]-dframe[i][j]) );
    for (i=0;i<frameLen;i++)
      for (j=0;j<frameLen;j++)
        d2 = mmdb::RMax ( d2,fabs(df[j][i]-dframe[i][j]) );

  }


  #ifdef __debug
  void out4 ( mmdb::mat44 & T, mmdb::pstr name )  {
    printf ( " -- %s\n",name );
    printf ( " %10.3f %10.3f %10.3f    %10.3f\n"
             " %10.3f %10.3f %10.3f    %10.3f\n"
             " %10.3f %10.3f %10.3f    %10.3f\n"
             " %10.3f %10.3f %10.3f    %10.3f\n",
          T[0][0],T[0][1],T[0][2],T[0][3],
          T[1][0],T[1][1],T[1][2],T[1][3],
          T[2][0],T[2][1],T[2][2],T[2][3],
          T[3][0],T[3][1],T[3][2],T[3][3] );
  }
  #endif

  mmdb::realtype vdistance ( mmdb::mat44  & t, mmdb::realtype x,
                             mmdb::realtype y, mmdb::realtype z )  {
  mmdb::realtype vx,vy,vz;
    vx = t[0][0]*x + t[0][1]*y + t[0][2]*z + t[0][3];
    vy = t[1][0]*x + t[1][1]*y + t[1][2]*z + t[1][3];
    vz = t[2][0]*x + t[2][1]*y + t[2][2]*z + t[2][3];
    return  sqrt ( vx*vx + vy*vy + vz*vz );
  }



  bool Interface::isSimilar ( PDomain A,
                                  mmdb::mat44 &    Ta,
                                  mmdb::realtype   rmsdA,
                                  PDomain B,
                                  mmdb::mat44 &    Tb,
                                  mmdb::realtype   rmsdB,
                                  mmdb::mat44 &    Tm )  {
  //
  //   Suppose there is an interface made of structures A and B,
  // such that structure A is superposable with [domain1] and
  // structure B - with [domain2] through the transformation
  // matrices Ta and Tb :
  //
  //    [domain1] = Ta*[A] + O(rmsdA)
  //    [domain2] = Tb*[B] + O(rmsdB)
  //
  // Tm is the transformation matrix of that interface, i.e. being
  // applied to B, brings it into interface with A:
  //
  //    Tm*[B]             makes an interface with [A]
  //    TMatrix*[domain2]  makes an interface with [domain1]
  //
  // This function returns true if interface between A and B is
  // structurally similar to that between [domain1] and [domain2].
  //
  mmdb::mat44     t1,t2;
  mmdb::realtype  R,dR, r1,r2,r3;
  int       i,j;

    //  t1 brings [B] into contact with [domain1] through superposing it
    //  on [domain2]
    mmdb::Mat4Mult ( t1,TMatrix,Tb );
    //  t2 brings [B] into contact with [domain1] through superposing
    //  A on [domain1]
    mmdb::Mat4Mult ( t2,Ta,Tm );

    //  If interfaces are structurally similar, t1 should be equal to t2.
    //  Check this by multiplying their difference with three
    //  representative vectors.

    R = (A->R+B->R)/2.0;
    if (R<=0.0)  R = 10.0;

    for (i=0;i<3;i++)
      for (j=0;j<4;j++)
        t1[i][j] -= t2[i][j];

    r1 = vdistance ( t1, B->mx+R,B->my  ,B->mz   );
    r2 = vdistance ( t1, B->mx  ,B->my+R,B->mz   );
    r3 = vdistance ( t1, B->mx  ,B->my  ,B->mz+R );

    dR = 0.075*R + 3.0*(rmsdA+rmsdB);

    return (r1<dR) && (r2<dR) && (r3<dR);

  }


  int Interface::isEquivalent ( PInterface Interface,
                                 PDomains Domains )  {
  //   This function returns non-zero if interface 'Interface' is
  // equivalent to 'this' one. Both interfaces belong to the same
  // crystal here. The interfaces are identified as equivalent if
  //   a) they are formed by highly similar or identical domains
  //   b) the domains are found in the same relative positions
  //      in both interfaces.
  //   'Domains' provides domain description for both interfaces,
  // which must include the structure superposition data.
  //   The function returns
  //      +1    - if interfaces match domain-to-domain
  //      -1    - if interfaces match inversely, i.e. if
  //              Interface->domain1 should be equivalenced with
  //              this->domain2 and Interface->domain2 - with
  //              this->domain1
  //       0    - if interface equivalence has not been identified
  //

    //   Interfaces between single-atom structures without symmetry
    // are not equivalenced here. The following statement makes
    // a trivial check, which should be good enough for the way
    // of ligand treatment in PISA (partial fixing to macromolecules).
    // A rigorous approach should include looking for equivalence
    // between single-atom ligands together with macromolecules they
    // attach to.

    if ((Domains->getNAtoms(domain1)<=1) &&
        (Domains->getNAtoms(domain2)<=1))  return 0;

    if (Domains->isEquivalent(domain1,Interface->domain1) &&
        Domains->isEquivalent(domain2,Interface->domain2))  {

      if (checkDFrame(Interface->dframe,true))  return 1;

    }

    if (Domains->isEquivalent(domain1,Interface->domain2) &&
        Domains->isEquivalent(domain2,Interface->domain1))  {

      if (checkDFrame(Interface->dframe,false))  return -1;

    }

    return 0;

  }


  void addXML0 ( mmdb::xml::PXMLObject xml, mmdb::cpstr tag,
                 mmdb::realtype V )  {
    if (fabs(V)>1.0e-15)
         xml->AddObject ( new mmdb::xml::XMLObject(tag,V)   );
    else xml->AddObject ( new mmdb::xml::XMLObject(tag,0.0) );
  }


  mmdb::PAtom GetAtom ( mmdb::PManager MMDB, int serNum )  {
  mmdb::PPAtom A;
  mmdb::PAtom  atom;
  int          nAtoms,i;

    if (!MMDB)  return NULL;

    MMDB->GetAtomTable ( A,nAtoms );

    if ((serNum>0) && (serNum<=nAtoms))  {
      // by design, this should be a direct hit (all serial numbers
      // get renumbered when the interface is calculated)
      atom = A[serNum-1];
      if (atom)  {
        if (atom->serNum==serNum)  return atom;
      }
    }

    // this part of code should never be engaged, but we do this
    // for safety
    atom = NULL;
    for (i=0;(i<nAtoms) && (!atom);i++)
      if (A[i])  {
        if (A[i]->serNum==serNum)  atom = A[i];
      }

    return atom;  // so that it can return NULL

  }


  mmdb::xml::PXMLObject Interface::getStructXML ( mmdb::PManager MMDB,
                                         PPDomain     Domain,
                                         int           structNo,
                                         int           dclass,
                                         int           symop,
                                         int           icell,
                                         int           jcell,
                                         int           kcell,
					 int           asis_param,
                                         mmdb::mat44 &       t,
                                         int           nIntAtoms,
                                         int           nIntRes,
                                         mmdb::realtype      intArea,
                                         mmdb::realtype      intDeltaG,
                                         mmdb::realtype      PValue,
                                         mmdb::cpstr         symOpn,
                                         mmdb::rvector       bsa,
                                         mmdb::rvector       solvEn )  {
  mmdb::xml::PXMLObject xml,xmlRes,xmli;
  PDomain    D;
  mmdb::PPResidue  Res;
  mmdb::PAtom      atom;
  char        S[500];
  int         selHndRes,selHndResS,selHndResI;
  int         selHndResH,selHndResB,selHndResD;
  int         selHndResC;
  int         selHndResOC; //GDL: add Other Contacts  
  int         i,nRes;


  
    
    xml = new mmdb::xml::XMLObject ( xml_asmunit );

    if (structNo==1)  D = Domain[domain1];
                else  D = Domain[domain2];

    xml->AddObject ( new mmdb::xml::XMLObject ( xml_asmunit_id,structNo ) );

    xml->AddObject ( new mmdb::xml::XMLObject ( xml_asmunit_name,
                                      D->getDomainID(S) ) );

    switch (dclass)  {
      default :
      case DCLASS_Protein :  strcpy ( S,"Protein" ); break;
      case DCLASS_DNA     :  strcpy ( S,"DNA"     ); break;
      case DCLASS_RNA     :  strcpy ( S,"RNA"     ); break;
      case DCLASS_Ligand  :  strcpy ( S,"Ligand"  );
    }
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_asmunit_class,S ) );
    
    //GDL: if as-is flag is used, don't write the following information
    if(asis_param==0){
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_symop_no,symop  ) );
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_symop   ,symOpn ) );
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_cell_i  ,icell  ) );
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_cell_j  ,jcell  ) );
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_cell_k  ,kcell  ) );
    
    
      addXML0 ( xml,xml_asmunit_rxx,t[0][0] );
      addXML0 ( xml,xml_asmunit_rxy,t[0][1] );
      addXML0 ( xml,xml_asmunit_rxz,t[0][2] );
      addXML0 ( xml,xml_asmunit_tx ,t[0][3] );
      addXML0 ( xml,xml_asmunit_ryx,t[1][0] );
      addXML0 ( xml,xml_asmunit_ryy,t[1][1] );
      addXML0 ( xml,xml_asmunit_ryz,t[1][2] );
      addXML0 ( xml,xml_asmunit_ty ,t[1][3] );
      addXML0 ( xml,xml_asmunit_rzx,t[2][0] );
      addXML0 ( xml,xml_asmunit_rzy,t[2][1] );
      addXML0 ( xml,xml_asmunit_rzz,t[2][2] );
      addXML0 ( xml,xml_asmunit_tz ,t[2][3] );
    }
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_natoms,nIntAtoms ) );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_nres ,nIntRes    ) );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_area ,intArea    ) );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_se   ,intDeltaG  ) );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_pval ,PValue     ) );


    selHndRes  = MMDB->NewSelection();
    selHndResS = MMDB->NewSelection();
    selHndResI = MMDB->NewSelection();
    selHndResH = MMDB->NewSelection();
    selHndResB = MMDB->NewSelection();
    selHndResD = MMDB->NewSelection();
    selHndResC = MMDB->NewSelection();

    D->SelectDomain ( selHndRes,MMDB,mmdb::STYPE_RESIDUE,mmdb::SKEY_NEW,1,false );
    MMDB->GetSelIndex ( selHndRes,Res,nRes );

    MMDB->Select ( selHndResS,mmdb::STYPE_RESIDUE,
                   D->selHndSurf,mmdb::SKEY_NEW );

    if (structNo<=1)  i = selHndInt1;
                else  i = selHndInt2;
    MMDB->Select ( selHndResI,mmdb::STYPE_RESIDUE,i,mmdb::SKEY_NEW );

    for (i=0;i<nHBonds;i++)  {
      if (structNo<=1)  atom = GetAtom(MMDB,HBond[i].serNum1);
                  else  atom = GetAtom(MMDB,HBond[i].serNum2);
      if (atom)
        MMDB->SelectResidue ( selHndResH,atom->GetResidue(),
                              mmdb::STYPE_RESIDUE,mmdb::SKEY_OR,false );
    }
    MMDB->MakeSelIndex ( selHndResH );

    for (i=0;i<nSBridges;i++)  {
      if (structNo<=1)
            atom = GetAtom ( MMDB,SBridge[i].serNum1 );
      else  atom = GetAtom ( MMDB,SBridge[i].serNum2 );
      if (atom)
        MMDB->SelectResidue ( selHndResB,atom->GetResidue(),
                              mmdb::STYPE_RESIDUE,mmdb::SKEY_OR,false );
    }
    MMDB->MakeSelIndex ( selHndResB );
    //GDL add other contacts: start
    for (i=0;i<nOthBonds;i++)  {
      if (structNo<=1)
            atom = GetAtom ( MMDB,OthBond[i].serNum1 );
      else  atom = GetAtom ( MMDB,OthBond[i].serNum2 );
      if (atom)
        MMDB->SelectResidue ( selHndResOC,atom->GetResidue(),
                              mmdb::STYPE_RESIDUE,mmdb::SKEY_OR,false );
    }
    MMDB->MakeSelIndex ( selHndResOC );
    //GDL add other contacts: end 

    for (i=0;i<nDSBonds;i++)  {
      if (structNo<=1)
            atom = GetAtom ( MMDB,DSBond[i].serNum1 );
      else  atom = GetAtom ( MMDB,DSBond[i].serNum2 );
      if (atom)
        MMDB->SelectResidue ( selHndResD,atom->GetResidue(),
                              mmdb::STYPE_RESIDUE,mmdb::SKEY_OR,false );
    }
    MMDB->MakeSelIndex ( selHndResD );

    for (i=0;i<nCovBonds;i++)  {
      if (structNo<=1)
            atom = GetAtom ( MMDB,CovBond[i].serNum1 );
      else  atom = GetAtom ( MMDB,CovBond[i].serNum2 );
      if (atom)
        MMDB->SelectResidue ( selHndResC,atom->GetResidue(),
                              mmdb::STYPE_RESIDUE,mmdb::SKEY_OR,false );
    }
    MMDB->MakeSelIndex ( selHndResC );

    xmlRes = new mmdb::xml::XMLObject ( xml_interface_res );

    for (i=0;i<nRes;i++)  {

      xmli = new mmdb::xml::XMLObject ( xml_residue );
      xmli->AddObject ( new mmdb::xml::XMLObject ( xml_residue_serno,i+1 ) );
      xmli->AddObject ( new mmdb::xml::XMLObject ( xml_residue_name,
                                         Res[i]->GetResName() ) );
      xmli->AddObject ( new mmdb::xml::XMLObject ( xml_residue_seqnum,
                                         Res[i]->GetSeqNum() ) );
      xmli->AddObject ( new mmdb::xml::XMLObject ( xml_residue_label_seqnum,
                                         Res[i]->GetLabelSeqID() ) );
      xmli->AddObject ( new mmdb::xml::XMLObject ( xml_residue_inscode,
                                         Res[i]->GetInsCode() ) );

      S[0] = char(0);
      if (Res[i]->isInSelection(selHndResH))  strcat ( S,"H" );
      if (Res[i]->isInSelection(selHndResB))  strcat ( S,"S" );
      if (Res[i]->isInSelection(selHndResD))  strcat ( S,"D" );
      if (Res[i]->isInSelection(selHndResC))  strcat ( S,"C" );
      xmli->AddObject ( new mmdb::xml::XMLObject ( xml_residue_bonds,S ) );

      xmli->AddObject ( new mmdb::xml::XMLObject ( xml_residue_asa,D->asa[i] ) );
      xmli->AddObject ( new mmdb::xml::XMLObject ( xml_residue_bsa,bsa   [i] ) );
      xmli->AddObject ( new mmdb::xml::XMLObject ( xml_residue_se ,solvEn[i] ) );

      xmlRes->AddObject ( xmli );

    }

    xml->AddObject ( xmlRes );

    MMDB->DeleteSelection ( selHndRes  );
    MMDB->DeleteSelection ( selHndResS );
    MMDB->DeleteSelection ( selHndResI );
    MMDB->DeleteSelection ( selHndResH );
    MMDB->DeleteSelection ( selHndResB );
    MMDB->DeleteSelection ( selHndResD );
    MMDB->DeleteSelection ( selHndResC );

    return xml;

  }

  mmdb::xml::PXMLObject Interface::makeBondTXML ( PIBond Bond, int nBonds,
                                         mmdb::cpstr Tag,
                                         mmdb::PManager MMDB ) {
  mmdb::xml::PXMLObject xml,xmlb;
  mmdb::PAtom      atom;
  int         i;
  
    xml = new mmdb::xml::XMLObject ( Tag );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_bond_nbonds,nBonds ) );

    for (i=0;i<nBonds;i++)  {

      xmlb = new mmdb::xml::XMLObject ( xml_bond_bond );
      
      atom = GetAtom ( MMDB,Bond[i].serNum1 );
      //cout<<atom->serNum<<"\n";
      
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_chain1,
                                         atom->GetChainID() ) );
      //xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_atom_site_id1,
      //                                   atom->GetAtomSiteID() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_label_asym_id1,
                                         atom->GetLabelAsymID() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_orig_label_asym_id1,
                                         atom->GetOrigLabelAsymID() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_pdbx_sifts_xref_db_num_1,
                                         atom->GetUniprotnum() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_pdbx_sifts_xref_db_name_1,
                                         atom->GetUniprotName() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_pdbx_sifts_xref_db_acc_1,
                                         atom->GetUniprotAcc() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_res1,
                                         atom->GetResName() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_seqnum1,
                                         atom->GetSeqNum() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_label_seqnum1,
                                         atom->GetLabelSeqID() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_inscode1,
                                         atom->GetInsCode() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_atname1,atom->name ) );

      atom = GetAtom ( MMDB,Bond[i].serNum2 );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_chain2,
                                         atom->GetChainID() ) );
      //xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_atom_site_id2,
      //                                   atom->GetAtomSiteID() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_label_asym_id2,
                                         atom->GetLabelAsymID() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_orig_label_asym_id2,
                                         atom->GetOrigLabelAsymID() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_pdbx_sifts_xref_db_acc_2,
                                         atom->GetUniprotAcc() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_pdbx_sifts_xref_db_num_2,
                                         atom->GetUniprotnum() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_pdbx_sifts_xref_db_name_2,
					 atom->GetUniprotName() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_res2,
                                         atom->GetResName() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_seqnum2,
                                         atom->GetSeqNum() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_label_seqnum2,
                                         atom->GetLabelSeqID() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_inscode2,
                                         atom->GetInsCode() ) );
      xmlb->AddObject ( new mmdb::xml::XMLObject ( xml_bond_atname2,atom->name ) );

      xmlb->AddObject ( new mmdb::xml::XMLObject (xml_bond_dist,Bond[i].dist ) );

      xml->AddObject ( xmlb );

    }

    return xml;

  }

  mmdb::xml::PXMLObject Interface::getXML ( mmdb::PManager MMDB,
					    PPDomain Domain, int as_is_param )  {
  mmdb::xml::PXMLObject xml;
  mmdb::mat44       t;
    
    xml = new mmdb::xml::XMLObject ( xml_interface );

    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_id  ,id        ) );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_type,type      ) );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_nocc,nOcc      ) );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_area,intArea   ) );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_se  ,intDeltaG ) );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_pval,PValue    ) );
    xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_sten,stabEn    ) );
    
    //GDL start: if status 'as-is' ommit this info in the xml file  
    if(as_is_param==0){
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_css ,css       ) );
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_overlap,overlap) );
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_xrel ,Xrel     ) );
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_interface_fixed,fixedLigand) );
    }
    //GDL end: if status 'as-is' ommit this info in the xml file
    
    xml->AddObject ( makeBondTXML ( HBond,nHBonds,
                                    xml_bond_hbonds,MMDB) );
    xml->AddObject ( makeBondTXML ( SBridge,nSBridges,
                                    xml_bond_sbridges,MMDB) );
    xml->AddObject ( makeBondTXML ( DSBond,nDSBonds,
                                    xml_bond_dsbonds,MMDB) );
    xml->AddObject ( makeBondTXML ( CovBond,nCovBonds,
                                    xml_bond_covbonds,MMDB) );
    xml->AddObject ( makeBondTXML ( OthBond,nOthBonds,
                                    xml_bond_otherbonds,MMDB) ); //GDL:add other contacts to xml 

    mmdb::Mat4Init ( t );
    xml->AddObject ( getStructXML ( MMDB,Domain,1,dclass1,1,
			      0,0,0,as_is_param,t,nIntAtoms1,nIntRes1,
                              intArea1,intDeltaG1,PValue1,"x,y,z",
                              bsa1,solvEn1 ) );

    xml->AddObject ( getStructXML ( MMDB,Domain,2,dclass2,rcsb_symop,
			      cell_i,cell_j,cell_k,as_is_param,TMatrix,
                              nIntAtoms2,nIntRes2,
                              intArea2,intDeltaG2,PValue2,symOp,
                              bsa2,solvEn2 ) );

    return xml;

  }

  void copySREffects ( RPSREffect  sre,
                       int       & nSRE,
                       PSREffect   sre_src,
                       int         nSRE_src )  {
    nSRE = nSRE_src;
    if (nSRE>0)  {
      sre = new SREffect[nSRE];
      for (int i=0;i<nSRE;i++)
        sre[i].Copy ( sre_src[i] );
    }
  }


  void Interface::Copy  ( PInterface interface )  {
  int i,j;

    type          = interface->type;
    domain1       = interface->domain1;
    domain2       = interface->domain2;
    dclass1       = interface->dclass1;
    dclass2       = interface->dclass2;
    symOpNo       = interface->symOpNo;
    rcsb_symop    = interface->rcsb_symop;
    cell_i        = interface->cell_i;
    cell_j        = interface->cell_j;
    cell_k        = interface->cell_k;
    nOcc          = interface->nOcc;
    mmdb::Mat4Copy ( interface->TMatrix,TMatrix );
    for (i=0;i<frameLen;i++)
      for (j=0;j<frameLen;j++)
        dframe[i][j] = interface->dframe[i][j];
    selHndInt1    = interface->selHndInt1;
    selHndInt2    = interface->selHndInt2;
    nIntAtoms1    = interface->nIntAtoms1;
    nIntAtoms2    = interface->nIntAtoms2;
    nIntRes1      = interface->nIntRes1;
    nIntRes2      = interface->nIntRes2;
    intArea1      = interface->intArea1;
    intArea2      = interface->intArea2;
    intArea       = interface->intArea;
    intDeltaG1    = interface->intDeltaG1;
    intDeltaG2    = interface->intDeltaG2;
    intDeltaG     = interface->intDeltaG;
    aveDeltaG1    = interface->aveDeltaG1;
    aveDeltaG2    = interface->aveDeltaG2;
    aveDeltaG     = interface->aveDeltaG;
    probDeltaG1   = interface->probDeltaG1;
    probDeltaG2   = interface->probDeltaG2;
    probDeltaG    = interface->probDeltaG;
    PValue1       = interface->PValue1;
    PValue2       = interface->PValue2;
    PValue        = interface->PValue;
    stabEn        = interface->stabEn;
    css           = interface->css;
    overlap       = interface->overlap;
    casual        = interface->casual;
    Xrel          = interface->Xrel;
    fixedLigand   = interface->fixedLigand;
    mmdb::CreateCopy ( symOp,interface->symOp );

    deleteHSDBonds();

    nHBonds = interface->nHBonds;
    if (nHBonds>0)  {
      HBond = new IBond[nHBonds];
      for (i=0;i<nHBonds;i++)  {
        HBond[i].Copy ( interface->HBond[i] );
      }
    }

    nSBridges = interface->nSBridges;
    if (nSBridges>0)  {
      SBridge = new IBond[nSBridges];
      for (i=0;i<nSBridges;i++)  {
        SBridge[i].Copy ( interface->SBridge[i] );
      }
    }
    //GDL add copy other contacts : start
    nOthBonds = interface->nOthBonds;
    if (nOthBonds>0)  {
      OthBond = new IBond[nOthBonds];
      for (i=0;i<nOthBonds;i++)  {
        OthBond[i].Copy ( interface->OthBond[i] );
      }
    }
    //GDL add copy other contacts : end
    
    
    nDSBonds = interface->nDSBonds;
    if (nDSBonds>0)  {
      DSBond = new IBond[nDSBonds];
      for (i=0;i<nDSBonds;i++)  {
        DSBond[i].Copy ( interface->DSBond[i] );
      }
    }

    deleteCovBonds();

    nCovBonds = interface->nCovBonds;
    if (nCovBonds>0)  {
      CovBond = new IBond[nCovBonds];
      for (i=0;i<nCovBonds;i++)  {
        CovBond[i].Copy ( interface->CovBond[i] );
      }
    }

    deleteSREffects();

    copySREffects ( sre_stab1,nStab1,
                    interface->sre_stab1,interface->nStab1 );
    copySREffects ( sre_stab2,nStab2,
                    interface->sre_stab2,interface->nStab2 );
    copySREffects ( sre_destab1,nDestab1,
                    interface->sre_destab1,interface->nDestab1 );
    copySREffects ( sre_destab2,nDestab2,
                    interface->sre_destab2,interface->nDestab2 );

    mmdb::FreeVectorMemory ( bsa1   ,0 );
    mmdb::FreeVectorMemory ( bsa2   ,0 );
    mmdb::FreeVectorMemory ( solvEn1,0 );
    mmdb::FreeVectorMemory ( solvEn2,0 );
    nRes1 = interface->nRes1;
    nRes2 = interface->nRes2;
    if (nRes1>0)  {
      mmdb::GetVectorMemory ( bsa1   ,nRes1,0 );
      mmdb::GetVectorMemory ( solvEn1,nRes1,0 );
      for (i=0;i<nRes1;i++)  {
        bsa1   [i] = interface->bsa1   [i];
        solvEn1[i] = interface->solvEn1[i];
      }
    }
    if (nRes2>0)  {
      mmdb::GetVectorMemory ( bsa2   ,nRes2,0 );
      mmdb::GetVectorMemory ( solvEn2,nRes2,0 );
      for (i=0;i<nRes2;i++)  {
        bsa2   [i] = interface->bsa2   [i];
        solvEn2[i] = interface->solvEn2[i];
      }
    }

  }

  void Interface::write ( mmdb::io::RFile f )  {
  int i,j,Version, dcl1,dcl2;

    Version = 3;
    f.WriteInt ( &Version    );

    f.WriteInt ( &id         );
    f.WriteInt ( &type       );
    f.WriteInt ( &domain1    );
    f.WriteInt ( &domain2    );
    dcl1 = dclass1;
    dcl2 = dclass2;
    f.WriteInt ( &dcl1       );
    f.WriteInt ( &dcl2       );
    f.WriteInt ( &symOpNo    );
    f.WriteInt ( &rcsb_symop );
    f.WriteInt ( &cell_i     );
    f.WriteInt ( &cell_j     );
    f.WriteInt ( &cell_k     );
    f.WriteInt ( &nOcc       );
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        f.WriteReal ( &(TMatrix[i][j]) );
    for (i=0;i<frameLen;i++)
      for (j=0;j<frameLen;j++)
        f.WriteReal ( &(dframe[i][j]) );
    f.WriteInt  ( &selHndInt1    );
    f.WriteInt  ( &selHndInt2    );
    f.WriteInt  ( &nIntAtoms1    );
    f.WriteInt  ( &nIntAtoms2    );
    f.WriteInt  ( &nIntRes1      );
    f.WriteInt  ( &nIntRes2      );
    f.WriteReal ( &intArea1      );
    f.WriteReal ( &intArea2      );
    f.WriteReal ( &intArea       );
    f.WriteReal ( &intDeltaG1    );
    f.WriteReal ( &intDeltaG2    );
    f.WriteReal ( &intDeltaG     );
    f.WriteReal ( &aveDeltaG1    );
    f.WriteReal ( &aveDeltaG2    );
    f.WriteReal ( &aveDeltaG     );
    f.WriteReal ( &probDeltaG1   );
    f.WriteReal ( &probDeltaG2   );
    f.WriteReal ( &probDeltaG    );
    f.WriteReal ( &PValue1       );
    f.WriteReal ( &PValue2       );
    f.WriteReal ( &PValue        );
    for (i=0;i<nASPs;i++)  {
      f.WriteReal ( &(DeltaSAS1[i]) );
      f.WriteReal ( &(DeltaSAS2[i]) );
    }
    f.WriteReal ( &stabEn        );
    f.WriteReal ( &css           );

    f.CreateWrite ( symOp        );

    f.WriteInt ( &nHBonds );
    for (i=0;i<nHBonds;i++)
      HBond[i].write ( f );

    f.WriteInt ( &nSBridges );
    for (i=0;i<nSBridges;i++)
      SBridge[i].write ( f );
    
    //GDL add write other contacts: start
    f.WriteInt ( &nOthBonds );
    for (i=0;i<nOthBonds;i++)
      OthBond[i].write ( f );
    //GDL add write other contacts: end
    
    f.WriteInt ( &nDSBonds );
    for (i=0;i<nDSBonds;i++)
      DSBond[i].write ( f );

    f.WriteInt ( &nCovBonds );
    for (i=0;i<nCovBonds;i++)
      CovBond[i].write ( f );

    f.WriteInt ( &nRes1 );
    for (i=0;i<nRes1;i++)  {
      f.WriteReal ( &(bsa1   [i]) );
      f.WriteReal ( &(solvEn1[i]) );
    }

    f.WriteInt ( &nRes2 );
    for (i=0;i<nRes2;i++)  {
      f.WriteReal ( &(bsa2   [i]) );
      f.WriteReal ( &(solvEn2[i]) );
    }

    f.WriteInt ( &nStab1 );
    for (i=0;i<nStab1;i++)
      sre_stab1[i].write ( f );

    f.WriteInt ( &nStab2 );
    for (i=0;i<nStab2;i++)
      sre_stab2[i].write ( f );

    f.WriteInt ( &nDestab1 );
    for (i=0;i<nDestab1;i++)
      sre_destab1[i].write ( f );

    f.WriteInt ( &nDestab2 );
    for (i=0;i<nDestab2;i++)
      sre_destab2[i].write ( f );

    f.WriteBool ( &overlap     );
    f.WriteBool ( &casual      );
    f.WriteBool ( &Xrel        );
    f.WriteBool ( &fixedLigand );

  }

  void Interface::read  ( mmdb::io::RFile f )  {
  int i,j,Version, dcl1,dcl2;

    deleteHSDBonds ();
    deleteCovBonds ();
    deleteSREffects();
    mmdb::FreeVectorMemory ( bsa1   ,0 );
    mmdb::FreeVectorMemory ( bsa2   ,0 );
    mmdb::FreeVectorMemory ( solvEn1,0 );
    mmdb::FreeVectorMemory ( solvEn2,0 );

    f.ReadInt ( &Version   );

    f.ReadInt ( &id         );
    f.ReadInt ( &type       );
    f.ReadInt ( &domain1    );
    f.ReadInt ( &domain2    );
    f.ReadInt ( &dcl1       );
    f.ReadInt ( &dcl2       );
    dclass1 = (DOMAIN_CLASS)dcl1;
    dclass2 = (DOMAIN_CLASS)dcl2;
    f.ReadInt ( &symOpNo    );
    f.ReadInt ( &rcsb_symop );
    f.ReadInt ( &cell_i     );
    f.ReadInt ( &cell_j     );
    f.ReadInt ( &cell_k     );
    f.ReadInt ( &nOcc       );
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        f.ReadReal ( &(TMatrix[i][j]) );
    for (i=0;i<frameLen;i++)
      for (j=0;j<frameLen;j++)
        f.ReadReal ( &(dframe[i][j]) );
    f.ReadInt  ( &selHndInt1    );
    f.ReadInt  ( &selHndInt2    );
    f.ReadInt  ( &nIntAtoms1    );
    f.ReadInt  ( &nIntAtoms2    );
    f.ReadInt  ( &nIntRes1      );
    f.ReadInt  ( &nIntRes2      );
    f.ReadReal ( &intArea1      );
    f.ReadReal ( &intArea2      );
    f.ReadReal ( &intArea       );
    f.ReadReal ( &intDeltaG1    );
    f.ReadReal ( &intDeltaG2    );
    f.ReadReal ( &intDeltaG     );
    f.ReadReal ( &aveDeltaG1    );
    f.ReadReal ( &aveDeltaG2    );
    f.ReadReal ( &aveDeltaG     );
    f.ReadReal ( &probDeltaG1   );
    f.ReadReal ( &probDeltaG2   );
    f.ReadReal ( &probDeltaG    );
    f.ReadReal ( &PValue1       );
    f.ReadReal ( &PValue2       );
    f.ReadReal ( &PValue        );

    for (i=0;i<nASPs;i++)  {
      f.ReadReal ( &(DeltaSAS1[i]) );
      f.ReadReal ( &(DeltaSAS2[i]) );
    }

    f.ReadReal ( &stabEn        );
    f.ReadReal ( &css           );

    f.CreateRead ( symOp         );

    f.ReadInt  ( &nHBonds );
    if (nHBonds>0)  {
      HBond = new IBond[nHBonds];
      for (i=0;i<nHBonds;i++)
        HBond[i].read ( f );
    }

    f.ReadInt  ( &nSBridges );
    if (nSBridges>0)  {
      SBridge = new IBond[nSBridges];
      for (i=0;i<nSBridges;i++)
        SBridge[i].read ( f );
    }
    //GDL add read Other Contacts: start
    f.ReadInt  ( &nOthBonds );
    if (nOthBonds>0)  {
      OthBond = new IBond[nOthBonds];
      for (i=0;i<nOthBonds;i++)
        OthBond[i].read ( f );
    }
    // GDL add read Other Contacts: end
    f.ReadInt  ( &nDSBonds );
    if (nDSBonds>0)  {
      DSBond = new IBond[nDSBonds];
      for (i=0;i<nDSBonds;i++)
        DSBond[i].read ( f );
    }

    f.ReadInt  ( &nCovBonds );
    if (nCovBonds>0)  {
      CovBond = new IBond[nCovBonds];
      for (i=0;i<nCovBonds;i++)
        CovBond[i].read ( f );
    }

    f.ReadInt ( &nRes1 );
    if (nRes1>0)  {
      mmdb::GetVectorMemory ( bsa1   ,nRes1,0 );
      mmdb::GetVectorMemory ( solvEn1,nRes1,0 );
      for (i=0;i<nRes1;i++)  {
        f.ReadReal ( &(bsa1   [i]) );
        f.ReadReal ( &(solvEn1[i]) );
      }
    }

    f.ReadInt ( &nRes2 );
    if (nRes2>0)  {
      mmdb::GetVectorMemory ( bsa2   ,nRes2,0 );
      mmdb::GetVectorMemory ( solvEn2,nRes2,0 );
      for (i=0;i<nRes2;i++) {
        f.ReadReal ( &(bsa2   [i]) );
        f.ReadReal ( &(solvEn2[i]) );
      }
    }

    if (Version>2)  {

      f.ReadInt ( &nStab1 );
      if (nStab1>0)  {
        sre_stab1 = new SREffect[nStab1];
        for (i=0;i<nStab1;i++)
          sre_stab1[i].read ( f );
      }

      f.ReadInt ( &nStab2 );
      if (nStab2>0)  {
        sre_stab2 = new SREffect[nStab2];
        for (i=0;i<nStab2;i++)
          sre_stab2[i].read ( f );
      }

      f.ReadInt ( &nDestab1 );
      if (nDestab1>0)  {
        sre_destab1 = new SREffect[nDestab1];
        for (i=0;i<nDestab1;i++)
          sre_destab1[i].read ( f );
      }

      f.ReadInt ( &nDestab2 );
      if (nDestab2>0)  {
        sre_destab2 = new SREffect[nDestab2];
        for (i=0;i<nDestab2;i++)
          sre_destab2[i].read ( f );
      }

    }

    f.ReadBool ( &overlap     );
    f.ReadBool ( &casual      );
    f.ReadBool ( &Xrel        );
    f.ReadBool ( &fixedLigand );

  }

  MakeStreamFunctions(Interface)



  // ========================  Interfaces  ==========================

  Interfaces::Interfaces() : mmdb::QuickSort()  {
    InitInterfaces();
  }

  Interfaces::Interfaces (  mmdb::io::RPStream Object )
             : mmdb::QuickSort ( Object )  {
    InitInterfaces();
  }

  Interfaces::~Interfaces()  {
    FreeMemory();
  }

  void Interfaces::InitInterfaces()  {
    rcsb_symops = false; // true if rcsb symops have been assigned
    PI          = NULL;
    intStatus   = INTS_notCalculated; // interface status
    nInterfaces = 0;
    nITypes     = 0;     // number of equivalent interface types
    nAlloc      = 0;
  }

  void Interfaces::FreeMemory()  {
  int i;
    if (PI)  {
      for (i=0;i<nAlloc;i++)
        if (PI[i])  delete PI[i];
      delete[] PI;
    }
    PI          = NULL;
    nInterfaces = 0;
    nITypes     = 0;   // number of equivalent interface types
    nAlloc      = 0;
  }

  void Interfaces::Reset()  {
    nInterfaces = 0;
    nITypes     = 0;   // number of equivalent interface types
  }

  void Interfaces::checkAllocation()  {
  PPInterface PI1;
  int          nAlloc1,i;
    if (nInterfaces>=nAlloc)  {
      nAlloc1 = nAlloc + 100;
      PI1     = new PInterface[nAlloc1];
      for (i=0;i<nAlloc;i++)
        PI1[i] = PI[i];
      for (i=nAlloc;i<nAlloc1;i++)
        PI1[i] = NULL;
      if (PI)  delete[] PI;
      PI     = PI1;
      nAlloc = nAlloc1;
    }
  }


  mmdb::pstr printMat ( mmdb::pstr S, mmdb::pstr name,
                        mmdb::mat44 & T )  {
    sprintf ( S,"\n %s\n"
             "   %6.3f %6.3f %6.3f   %6.3f\n"
             "   %6.3f %6.3f %6.3f   %6.3f\n"
             "   %6.3f %6.3f %6.3f   %6.3f\n"
             "   %6.3f %6.3f %6.3f   %6.3f\n",
             name,T[0][0],T[0][1],T[0][2],T[0][3],
                  T[1][0],T[1][1],T[1][2],T[1][3],
                  T[2][0],T[2][1],T[2][2],T[2][3],
                  T[3][0],T[3][1],T[3][2],T[3][3] );
    return S;
  }

  bool TMatEqual ( mmdb::mat44 & T1, mmdb::mat44 & T2 )  {
  // Returns true if T1 equals to T2
  mmdb::realtype dT,dR, d;
  int      i;

    dT = 0.0;
    dR = 0.0;

    for (i=0;i<3;i++)  {
      d   = T1[i][3] - T2[i][3];
      dT += d*d;
      d   = T1[i][0] - T2[i][0];
      dR += d*d;
      d   = T1[i][1] - T2[i][1];
      dR += d*d;
      d   = T1[i][2] - T2[i][2];
      dR += d*d;
    }

    return ((dR<rot_threshold) && (dT<trans_threshold));

  }


  void Interfaces::AddContact ( PDomains D, int domain1, int domain2,
                                int nSymOpNo, mmdb::mat44 & TMatrix,
                                int cell_i, int cell_j, int cell_k )  {
  DFrame   dframe;
  PDomain D1,D2;
  int      i, pd1,pd2, ipd1,ipd2,iequiv;

    //  Check that the contact is a new one

    D1  = D->domain[domain1];
    D2  = D->domain[domain2];
    pd1 = D1->ncsParent;
    pd2 = D2->ncsParent;
    MakeDFrame ( dframe,D1,D2,TMatrix );

    iequiv = -1;
    for (i=0;(i<nInterfaces) && (iequiv<0);i++)  {
      ipd1 = D->domain[PI[i]->domain1]->ncsParent;
      ipd2 = D->domain[PI[i]->domain2]->ncsParent;
      if ((ipd1==pd1) && (ipd2==pd2))  {
        if (PI[i]->checkDFrame(dframe,true))  iequiv = i;
      }
      if ((iequiv<0) && (ipd1==pd2) && (ipd2==pd1))  {
        if (PI[i]->checkDFrame(dframe,false))  iequiv = i;
      }
    }

    if (iequiv<0)  {
      checkAllocation();
      if (!PI[nInterfaces])  PI[nInterfaces] = new Interface();
                       else  PI[nInterfaces]->Reset();
      PI[nInterfaces]->SetContact ( domain1,domain2,D,nSymOpNo,
                                    TMatrix,cell_i,cell_j,cell_k );
      nInterfaces++;
    } else
      PI[iequiv]->nOcc++;

  }


  RESULT_CODE Interfaces::CalcContacts ( mmdb::PManager MMDB,
                                         PDomains       D,
                                         int            nCellLayers )  {
  PBrick        brick[27];
  mmdb::mat44   TMatrix;
  mmdb::ivector selHnd;
  int           i,j,k, ic,jc,kc,m, n, nSymOps;
  RESULT_CODE   rc;

    rc = RESULT_Ok;

    if (MMDB->CrystReady()>=0)  {
      nSymOps = MMDB->GetNumberOfSymOps();
      if (nSymOps<=0)
        return RESULT_noSymOps; // no symmetry operations, although they
                                // should be there
    } else
      nSymOps = 0;

  //  nNCSOps = MMDB->GetNumberOfNCSMatrices();

    mmdb::GetVectorMemory ( selHnd,D->nDomains,0 );
    for (i=0;i<D->nDomains;i++)
      selHnd[i] = 0;

    D->MakeBricks();

    if (nSymOps<=0)  {
      // check only contacts in ASU

      mmdb::Mat4Init ( TMatrix );
      for (i=0;i<D->nDomains;i++)  {
        D->getBricks ( brick,D->domain[i],TMatrix );
        if (brick[0])  {
          D->domain[i]->MakeContactBricks ( MMDB,TMatrix,selHnd[i] );
          for (k=0;(k<27) && brick[k];k++)
            for (j=0;j<brick[k]->nObjects;j++)  {
              n = brick[k]->n[j];
              if (n!=i)  {
                if (D->domain[i]->isContact(TMatrix,D->domain[n],
                                            MMDB,selHnd[n]))
                  AddContact ( D,i,n,-1,TMatrix,0,0,0 );
              }
            }
        }
      }

    } else  {
      // check contacts in +/- nCellLayers neighbouring unit cells

      for (m=0;(m<nSymOps) && (!rc);m++)
        if (MMDB->GetTMatrix(TMatrix,m,0,0,0)!=mmdb::SYMOP_Ok)  {
          rc = RESULT_noSymOp; // no symmetry operation,
                               // although it should be there
        } else if (isMat4Rot(TMatrix,1.0e-2))  {
          for (ic=-nCellLayers;(ic<=nCellLayers) && (!rc);ic++)
            for (jc=-nCellLayers;(jc<=nCellLayers) && (!rc);jc++)
              for (kc=-nCellLayers;(kc<=nCellLayers) && (!rc);kc++)
                //   GetMatrix(..) calculates the transformation
                // matrix for mth symmetry operation, which
                // places atoms into unit cell, shifted by
                // ic,jc,kc in a,b,c - directions, respectively
                // (in fractional space), from the principal
                // unit cell of the coordinate file.
                if (MMDB->GetTMatrix(TMatrix,m,ic,jc,kc) !=
                    mmdb::SYMOP_Ok)  {
                    rc = RESULT_noSymOp; // no symmetry operation,
                                         // although it should be there
                } else  {
                  for (i=0;i<D->nDomains;i++)  {
                    D->getBricks ( brick,D->domain[i],TMatrix );
                    if (brick[0])  {
                      D->domain[i]->MakeContactBricks (
                                          MMDB,TMatrix,selHnd[i] );
                      for (k=0;(k<27) && brick[k];k++)
                        for (j=0;j<brick[k]->nObjects;j++)  {
                          n = brick[k]->n[j];
                          if (m || ic || jc || kc || (n!=i))  {
                            if (D->domain[i]->isContact(TMatrix,
                                      D->domain[n],MMDB,selHnd[n]))
                              AddContact ( D,n,i,m,TMatrix,ic,jc,kc );
                          }
                        }
                    }
                  }
                }
        }

    }

    D->RemoveBricks();

    for (i=0;i<D->nDomains;i++)
      if (selHnd[i])  MMDB->DeleteSelection ( selHnd[i] );
     mmdb::FreeVectorMemory ( selHnd,0 );

    return rc;

  }


  void Interfaces::AddInterface ( PInterface Interface )  {
    checkAllocation();
    if (PI[nInterfaces])  delete PI[nInterfaces];
    PI[nInterfaces] = Interface;
    nInterfaces++;
  }

  void  Interfaces::CalcInterfaceTypes ( PDomains domains )  {
  int i,j;

    nITypes = 0;

    for (i=0;i<nInterfaces;i++)
      PI[i]->type = 0;

    for (i=0;i<nInterfaces;i++)
      if (PI[i]->type<=0)  {
        nITypes++;
        PI[i]->type = nITypes;
        for (j=i+1;j<nInterfaces;j++)
          if (PI[j]->type<=0)  {
            if (PI[i]->isEquivalent(PI[j],domains))
              PI[j]->type = nITypes;
          }
      }

  }

  int Interfaces::AnalyseInterfaces ( PDomains Domains )  {
  //   Returns general information about interfaces: whether
  // they are there and whether there are overlapping macromolecular
  // and/or ligand interfaces.
  PPDomain D;
  int       i;

    if (nInterfaces<=0)
      intStatus = INTS_noInterfaces;
    else  {
      D = Domains->domain;
      intStatus = INTS_Ok;
      for (i=0;i<nInterfaces;i++)
        if (PI[i]->overlap)  {
          if ((D[PI[i]->domain1]->dclass!=DCLASS_Ligand) &&
              (D[PI[i]->domain2]->dclass!=DCLASS_Ligand))
               intStatus |= INTS_mmOverlap;
          else intStatus |= INTS_ligandOverlap;
        }
    }

    return intStatus;

  }

  PInterface Interfaces::getInterface ( int intfNo )  {
    if ((0<=intfNo) && (intfNo<nInterfaces))
         return PI[intfNo];
    else return NULL;
  }

  int Interfaces::Compare ( int i, int j )  {
  bool gt,lt;
    switch (smode)  {
      default           :
      case ISORT_Area   : // sort by decreasing the interface area
                          gt = PI[i]->intArea < PI[j]->intArea;
                          lt = PI[i]->intArea > PI[j]->intArea;
                        break;
      case ISORT_DeltaG : // sort by decreasing the interface
                          // solvation effect
                          gt = PI[i]->intDeltaG < PI[j]->intDeltaG;
                          lt = PI[i]->intDeltaG > PI[j]->intDeltaG;
                        break;
      case ISORT_IType  : // sort by increasing the interface type
                          gt = PI[i]->type > PI[j]->type;
                          lt = PI[i]->type < PI[j]->type;
    }
    if ((!lt) && (!gt))  {
      gt = PI[i]->intArea < PI[j]->intArea;
      lt = PI[i]->intArea > PI[j]->intArea;
    }
    if (gt)  return  1;
    if (lt)  return -1;
    return 0;
  }

  void Interfaces::Swap ( int i, int j )  {
  PInterface I;
    I     = PI[i];
    PI[i] = PI[j];
    PI[j] = I;
  }


  void Interfaces::Sort ( ISORT_KEY sortmode )  {
  int i,j;

    smode = sortmode;
    if (smode!=ISORT_Off)
      QuickSort::Sort ( &(PI[0]),nInterfaces );

    for (i=0;i<nInterfaces;i++)
      PI[i]->type = -PI[i]->type;

    nITypes = 0;
    for (i=0;i<nInterfaces;i++)
      if (PI[i]->type<0)  {
        nITypes++;
        for (j=i+1;j<nInterfaces;j++)
          if (PI[j]->type==PI[i]->type)
            PI[j]->type = nITypes;
        PI[i]->type = nITypes;
      }

    smode = ISORT_IType;
    QuickSort::Sort ( &(PI[0]),nInterfaces );
    smode = sortmode;

    for (i=0;i<nInterfaces;i++)
      PI[i]->id = i+1;

  }

  void Interfaces::DeleteDummyInterfaces()  {
  PInterface Interface;
  int         i,k;
    if (PI)  {
      k = 0;
      for (i=0;i<nInterfaces;i++)
        if (PI[i])  {
          if ((PI[i]->intArea1>0.0) && (PI[i]->intArea2>0.0))  {
            if (i>k)  {
              Interface = PI[k];
              PI[k]     = PI[i];
              PI[i]     = Interface;
            }
            k++;
          }
        }
      nInterfaces = k;
    }
  }


  int Interfaces::getSimilarInterface (
        PDomain A, mmdb::mat44 & Ta, int domain1, mmdb::realtype rmsdA,
        PDomain B, mmdb::mat44 & Tb, int domain2, mmdb::realtype rmsdB,
        mmdb::mat44 & Tm, PDomains Domains )  {
  //
  //   Suppose there is an interface between structure [A], which is
  // similar to [domain1], and structure [B], which is similar to
  // [domain2] :
  //
  //     [domain1] = Ta*[A] + O(rmsdA)
  //     [domain2] = Tb*[B] + O(rmsdB)
  //
  // Matrix Tm makes and interface between [A] and [B] :
  //
  //     Tm*[B] brings [B] into an interface with [A]
  //
  // This function tries to find an interface in 'this' crystal
  // (i.e. between [Domains]) that is structurally similar to
  // that between A and B. "Structurally similar" means:
  //
  //    a) the interface is made by [domain1] and [domain2] *or*
  //       structures similar to them
  //    b) the interfacing structures are found in the same relative
  //       positions as [A] and Tm*[B]
  //
  // The function returns the interface (serial) number
  // (0..nInterfaces-1) or negative if none is found.
  //
  // Domains refers to domains of 'this' crystal.
  //
  mmdb::mat44 Tm1;
  mmdb::mat44 Tas,Tbs, TTa,TTb;
  int   i,k, d1,d2;

    if (!PI)  return -12;

    mmdb::Mat4Inverse ( Tm,Tm1 );  // for checking on swapped structures

    k = -1;
    for (i=0;(i<nInterfaces) && (k<0);i++)
      if (PI[i])  {

        d1 = PI[i]->domain1;
        d2 = PI[i]->domain2;

        if (Domains->isEquivalent(domain1,d1) &&
            Domains->isEquivalent(domain2,d2))  {
          Domains->getTMatrix ( Tas,domain1,d1 ); // Tas*[domain1]->[d1]
          Domains->getTMatrix ( Tbs,domain2,d2 ); // Tbs*[domain2]->[d2]
          mmdb::Mat4Mult ( TTa,Tas,Ta );  // TTa*[A] = Tas*Ta^{-1}*[A] -> [d1]
          mmdb::Mat4Mult ( TTb,Tbs,Tb );  // TTb*[B] = Tbs*Tb^{-1}*[B] -> [d2]
          if (PI[i]->isSimilar(A,TTa,rmsdA+Domains->getRMSD(domain1,d1),
                               B,TTb,rmsdB+Domains->getRMSD(domain2,d2),
                               Tm))
            k = i;
        }

        if ((k<0) && Domains->isEquivalent(domain1,d2) &&
                     Domains->isEquivalent(domain2,d1))  {
          Domains->getTMatrix ( Tas,domain1,d2 ); // Tas*[domain1]->[d2]
          Domains->getTMatrix ( Tbs,domain2,d1 ); // Tbs*[domain2]->[d1]
          mmdb::Mat4Mult ( TTa,Tas,Ta );  // TTa*[A] = Tas*Ta^{-1}*[A] -> [d2]
          mmdb::Mat4Mult ( TTb,Tbs,Tb );  // TTb*[B] = Tbs*Tb^{-1}*[B] -> [d1]
          if (PI[i]->isSimilar(B,TTb,rmsdB+Domains->getRMSD(domain1,d2),
                               A,TTa,rmsdA+Domains->getRMSD(domain2,d1),
                               Tm1))
            k = i;
        }

      }

    return k;

  }

  void Interfaces::assignRCSBSymOps ( PRCSBData rcsbData )  {
  int i;
    if (rcsbData->nSymOps>0)  {
      for (i=0;i<nInterfaces;i++)
        PI[i]->rcsb_symop = rcsbData->rcsb_symop[PI[i]->symOpNo];
      rcsb_symops = true;
    }
  }


  mmdb::xml::PXMLObject Interfaces::getXML ( mmdb::PManager MMDB,
                                             PPDomain Domain, int as_is_param ) {
  mmdb::xml::PXMLObject xml;
  int         i;
  
  //cout<<"molecular weight"<<"\t"<< weight << "\n";
    xml = new mmdb::xml::XMLObject ( xml_pdb_entry );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_pdb_code,
                                              MMDB->GetEntryID()) );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_status,"Ok") );

    if (PI)  {
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_ints_nints,
                                                  nInterfaces ) );
      for (i=0;i<nInterfaces;i++)
        xml->AddObject ( PI[i]->getXML ( MMDB,Domain,as_is_param ) );
    } else
      xml->AddObject ( new mmdb::xml::XMLObject ( xml_ints_nints,0 ) );

    return xml;

  }


  void Interfaces::write ( mmdb::io::RFile f )  {
  int i,Version;
    Version = 2;
    f.WriteInt  ( &Version     );
    f.WriteBool ( &rcsb_symops );
    f.WriteWord ( &intStatus   );
    f.WriteInt  ( &nInterfaces );
    f.WriteInt  ( &nITypes     );
    for (i=0;i<nInterfaces;i++)  {
      PI[i]->id = i+1;
      StreamWrite ( f,PI[i] );
    }
  }

  void Interfaces::read ( mmdb::io::RFile f )  {
  int        i,Version;
    FreeMemory();
    f.ReadInt  ( &Version     );
    if (Version>1)
         f.ReadBool ( &rcsb_symops );
    else rcsb_symops = false;
    f.ReadWord ( &intStatus   );
    f.ReadInt  ( &nInterfaces );
    f.ReadInt  ( &nITypes     );
    if (nInterfaces>0)  {
      nAlloc = nInterfaces;
      PI = new PInterface[nAlloc];
      for (i=0;i<nAlloc;i++)  {
        PI[i] = NULL;
        StreamRead ( f,PI[i] );
        PI[i]->id = i+1;
      }
    }
  }

  MakeStreamFunctions(Interfaces)



  bool isFixed ( int domainId, PPInterface I, int nInterfaces )  {
  int     j;
  bool fixed = false;

    for (j=0;(j<nInterfaces) && (!fixed);j++)
      if (I[j]->fixedLigand)
        fixed = (I[j]->domain1==domainId) ||
                (I[j]->domain2==domainId);

    return fixed;

  }

}  // namespace pisa
