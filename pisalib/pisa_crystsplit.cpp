// $Id: pisa_crystsplit.cpp $
// =================================================================
//
//    06.12.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_crystsplit <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::CrystSplit
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2004-2013
//
// =================================================================
//


#include "pisa_crystsplit.h"
#include "pisa_defs.h"

namespace pisa  {

  // =========================  CrystSplit  =========================

  CrystSplit::CrystSplit()  {
    intf        = NULL;  // vector of engaged interfaces
    A           = NULL;  // list of assemblies obtained
    nAssemblies = 0;     // number of assemblies obtained
    Score       = 0;     // score 0-7
    equiv_all   = false; // true if all assemblies are equivalent
    stable_all  = false;
    orig_chains = false; // true if all ASU chains are used
  }

  CrystSplit::~CrystSplit()  {
    FreeMemory();
  }

  void CrystSplit::FreeMemory()  {
  int i;
    mmdb::FreeVectorMemory ( intf,0 );
    if (A)  {
      for (i=0;i<nAssemblies;i++)
        if (A[i])  delete A[i];
      delete[] A;
      A = NULL;
    }
    nAssemblies = 0;
  }

  void CrystSplit::Copy ( PMultimerSet multSet, PPInterface Interface,
                          int nInterfaces, int nInterfaces0,
                          PPDomain D, mmdb::mat44 & rom )  {
  int i,j,k,n;

    FreeMemory();

    nAssemblies = multSet->nMultimers;

    if (nAssemblies>0)  {

      A = new PAssembly[nAssemblies];
      for (i=0;i<nAssemblies;i++)  {
        A[i] = new Assembly();
        A[i]->Copy ( multSet->U[i],Interface,nInterfaces,
                     nInterfaces0,D,rom );
      }

    }

    if (nInterfaces0>0)  {
      mmdb::GetVectorMemory ( intf,nInterfaces0,0 );
      for (i=0;i<nInterfaces0;i++)  {
        intf[i] = 0;
        for (j=0;j<nAssemblies;j++)
          for (k=0;k<A[j]->asmSize;k++)  {
            n = 0;
            while (n<_max_n_int)
              if (A[j]->M[k]->intfl[n])  {
                if (A[j]->M[k]->intfl[n][i]>=0)  intf[i]++;
                n++;
              } else
                n = _max_n_int;
          }
        intf[i] /= 2;
      }
    }

    equiv_all  = multSet->equiv_all;
    stable_all = multSet->stable_all;

  }

  int CrystSplit::getMaxAsmSize()  {
  int i,m;
    m = 0;
    for (i=0;i<nAssemblies;i++)
      m = mmdb::IMax ( m,A[i]->asmSize );
    return m;
  }


  void CrystSplit::checkOriginalOrientations ( PDomains Domains,
                                               mmdb::ivector   icnt )  {
  PPDomain D;
  int      i;

    for (i=0;i<Domains->nNCSParents;i++)
      icnt[i] = 0;

    D = Domains->domain;

    for (i=0;i<nAssemblies;i++)
      A[i]->countOriginalOrientations ( D,icnt );

    orig_chains = true; // true if all ASU chains are used
    for (i=0;(i<Domains->nNCSParents) && orig_chains;i++)
      if (D[i]->dclass!=DCLASS_Ligand)
        orig_chains = (icnt[i]>0);

  }

  void CrystSplit::calcScore()  {
  mmdb::realtype b,bmin,bmax,bmax1;
  int      i,smin,smax;

    bmin  =  mmdb::MaxReal;
    bmax  = -mmdb::MaxReal;
    bmax1 = -mmdb::MaxReal;
    smin  =  mmdb::MaxInt;
    smax  =  0;
    for (i=0;i<nAssemblies;i++)
      if (A[i]->mmSize>0)  {
        b = -A[i]->freeEn;
        if (b<bmin)  bmin = b;
        if (b>bmax)  bmax = b;
        if ((A[i]->mmSize>1) && (b>bmax1))  bmax1 = b;
        if (A[i]->mmSize<smin)  smin = A[i]->mmSize;
        if (A[i]->nDiss>1)  {
          if (A[i]->asmSize>smax) smax = A[i]->asmSize;
        } else  {
          if (A[i]->mmSize>smax)  smax = A[i]->mmSize;
        }
      }

    if (smin>1)  {
      if (equiv_all)  {
        if (bmax<0.0)       Score = 0;
        else if (bmax<2.0)  Score = 3;
                      else  Score = 8;
      } else  {
        if (bmax<0.0)    Score = 1;
        else if (bmax<2.0)  {
          if (bmin<0.0)  Score = 4;
                   else  Score = 5;
        } else           Score = 8;
      }
    } else if (smax>1)  {
      if (bmax1<0.0)      Score = 2;
      else if (bmin<0.0)  Score = 6;
      else if (bmax<2.0)  Score = 7;
                    else  Score = 8;
    } else
      Score = 8;

  }

  mmdb::xml::PXMLObject CrystSplit::getCrystSplitXML ( int      serNo,
                                                       PDomains     D,
                                                       PInterfaces PI,
                                                       int  nCellOut ) {
  mmdb::xml::PXMLObject xml;
  int         i;

    xml = new mmdb::xml::XMLObject ( xml_cryst_split );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_cryst_split_serno,serNo));
    if (orig_chains)
          xml->AddObject ( new mmdb::xml::XMLObject(xml_orig_chains,"Yes") );
    else  xml->AddObject ( new mmdb::xml::XMLObject(xml_orig_chains,"No" ) );
    for (i=0;i<nAssemblies;i++)
      xml->AddObject ( A[i]->getAssemblyXML(D,PI,nCellOut,Score) );

    return xml;

  }


  void CrystSplit::write ( mmdb::io::RFile f, int nInterfaces )  {
  int  i;
  mmdb::byte Version;

    Version = 2;
    f.WriteByte ( &Version );

    f.WriteInt ( &nInterfaces );
    for (i=0;i<nInterfaces;i++)
      f.WriteInt ( &(intf[i]) );

    f.WriteInt ( &nAssemblies );
    for (i=0;i<nAssemblies;i++)
      A[i]->write ( f,nInterfaces );

    f.WriteInt  ( &Score      );
    f.WriteBool ( &equiv_all  );
    f.WriteBool ( &stable_all );
    f.WriteBool ( &orig_chains );

  }

  void CrystSplit::read ( mmdb::io::RFile f, int nInterfaces )  {
  int  i;
  mmdb::byte Version;

    FreeMemory();

    f.ReadByte ( &Version );

    f.ReadInt ( &nInterfaces );
    if (nInterfaces>0)  {
      mmdb::GetVectorMemory ( intf,nInterfaces,0 );
      for (i=0;i<nInterfaces;i++)
        f.ReadInt ( &(intf[i]) );
    }

    f.ReadInt ( &nAssemblies );
    if (nAssemblies>0)  {
      A = new PAssembly[nAssemblies];
      for (i=0;i<nAssemblies;i++)  {
        A[i] = new Assembly();
        A[i]->read ( f,nInterfaces );
      }
    }

    f.ReadInt  ( &Score      );
    f.ReadBool ( &equiv_all  );
    f.ReadBool ( &stable_all );
    if (Version>1)
          f.ReadBool ( &orig_chains );
    else  orig_chains = true;

  }


}  // namespace pisa
