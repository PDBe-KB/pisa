// $Id: pisa_assemblies.cpp $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_assemblies <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Assemblies
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2004-2014
//
// =================================================================
//

#include <string.h>
#include <math.h>

#include "pisa_assemblies.h"
#include "pisa_defs.h"
#include "mmdb2/mmdb_tables.h"

namespace pisa  {


  // ========================  Assemblies  ==========================

  Assemblies::Assemblies() : mmdb::io::Stream()  {
    InitAssemblies();
  }

  Assemblies::Assemblies ( mmdb::io::RPStream Object )
            : mmdb::io::Stream ( Object )  {
    InitAssemblies();
  }

  Assemblies::~Assemblies()  {
    FreeMemory();
  }

  void Assemblies::InitAssemblies()  {
    crystSplit   = NULL;  // assembly sets obtained
    nCrystSplits = 0;     // total number of assembly sets
    nCSRes       = 0;     // number of assembly sets scored less than 8
    asmStock     = NULL;  // assembly stock (bulk solution)
    nInterfaces  = 0;     // number of interfaces
    rcsb_symops  = false; // true if rcsb symops have been assigned
    orig_chains  = false; // true if all ASU chains are used
  }

  void Assemblies::FreeMemory()  {
  int i;

    if (crystSplit) {
      for (i=0;i<nCrystSplits;i++)
        if (crystSplit[i])  delete crystSplit[i];
      delete[] crystSplit;
      crystSplit = NULL;
    }
    nCrystSplits = 0;
    nCSRes  = 0;

    if (asmStock)  {
      delete asmStock;
      asmStock = NULL;
    }

  }


  DefineClass(SortAssemblies);

  class SortAssemblies : public mmdb::QuickSort  {

    public :
      SortAssemblies ();
      ~SortAssemblies() {}
      int  Compare ( int i, int j );
      void Swap    ( int i, int j );
      void Sort    ( PCrystSplit crystSplit );

    protected :
      PPAssembly    A;
      mmdb::rvector b;

  };

  SortAssemblies::SortAssemblies() : mmdb::QuickSort() {}

  int SortAssemblies::Compare ( int i, int j )  {

    if ((b[j]>0.0) && (b[i]<0.0))  return 1;

    if (A[i]->mmSize<A[j]->mmSize)  return  1;
    if (A[i]->mmSize>A[j]->mmSize)  return -1;

    if (b[i]>b[j])  return  1;
    if (b[i]<b[j])  return -1;

    return 0;

  }

  void SortAssemblies::Swap ( int i, int j )  {
  PAssembly      Ai;
  mmdb::realtype bi;
    Ai   = A[i];
    A[i] = A[j];
    A[j] = Ai;
    bi   = b[i];
    b[i] = b[j];
    b[j] = bi;
  }

  void SortAssemblies::Sort ( PCrystSplit crystSplit )  {
  int i;
    if (crystSplit->nAssemblies>1)  {
      A = crystSplit->A;
      mmdb::GetVectorMemory ( b,crystSplit->nAssemblies,0 );
      for (i=0;i<crystSplit->nAssemblies;i++)
        b[i] = -A[i]->freeEn;
      mmdb::QuickSort::Sort ( A,crystSplit->nAssemblies );
      mmdb::FreeVectorMemory ( b,0 );
    }
  }



  DefineClass(SortcrystSplits);

  class SortcrystSplits : public mmdb::QuickSort  {

    public :
      SortcrystSplits ();
      ~SortcrystSplits() {}
      int  Compare ( int i, int j );
      void Swap    ( int i, int j );
      void Sort    ( PPCrystSplit crystSplit, int nCrystSplits );

    protected :
      PPCrystSplit  A;
      mmdb::rvector b;
      mmdb::ivector s,n;

  };

  SortcrystSplits::SortcrystSplits() : mmdb::QuickSort() {}

  int SortcrystSplits::Compare ( int i, int j )  {

    if (A[i]->Score<A[j]->Score)  return -1;
    if (A[i]->Score>A[j]->Score)  return  1;

    if (A[i]->equiv_all && (!A[j]->equiv_all))  return -1;
    if ((!A[i]->equiv_all) && A[j]->equiv_all)  return  1;

    if (A[i]->stable_all && (!A[j]->stable_all))  return -1;
    if ((!A[i]->stable_all) && A[j]->stable_all)  return  1;

    if (s[i]>s[j])  return -1;
    if (s[i]<s[j])  return  1;

    if (n[i]>n[j])  return  1;
    if (n[i]<n[j])  return -1;

    if (b[i]<b[j])  return -1;
    if (b[i]>b[j])  return  1;

    return 0;

  }

  void SortcrystSplits::Swap ( int i, int j )  {
  PCrystSplit        Ai;
  mmdb::realtype bi;
  int      si;
    Ai   = A[i];
    A[i] = A[j];
    A[j] = Ai;
    bi   = b[i];
    b[i] = b[j];
    b[j] = bi;
    si   = s[i];
    s[i] = s[j];
    s[j] = si;
    si   = n[i];
    n[i] = n[j];
    n[j] = si;
  }

  void SortcrystSplits::Sort ( PPCrystSplit crystSplit, int nCrystSplits )  {

  int i,j;

    for (i=0;i<nCrystSplits;i++)
      crystSplit[i]->calcScore();

    if (nCrystSplits>1)  {
      A = crystSplit;
      mmdb::GetVectorMemory ( b,nCrystSplits,0 );
      mmdb::GetVectorMemory ( s,nCrystSplits,0 );
      mmdb::GetVectorMemory ( n,nCrystSplits,0 );
      for (i=0;i<nCrystSplits;i++)  {
        b[i] = mmdb::MaxReal;
        s[i] = 0;
        n[i] = 0;
        for (j=0;j<A[i]->nAssemblies;j++)  {
          if (A[i]->A[j]->mmSize>0)  {
            b[i] = mmdb::RMin ( b[i],-A[i]->A[j]->freeEn );
            n[i]++;
          }
          s[i] = mmdb::IMax ( s[i],A[i]->A[j]->mmSize );
        }
        if (b[i]==mmdb::MaxReal)  b[i] = 0.0;
      }
      mmdb::QuickSort::Sort ( crystSplit,nCrystSplits );
      mmdb::FreeVectorMemory ( n,0 );
      mmdb::FreeVectorMemory ( s,0 );
      mmdb::FreeVectorMemory ( b,0 );
    }

  }


  void Assemblies::Sort ( PAssembly Complex )  {
  SortAssemblies SA;
  SortcrystSplits    AS;
  int             i,j,k,m,tij,t,serNo;

    //  sort assemblies in assembly sets
    for (i=0;i<nCrystSplits;i++)
      SA.Sort ( crystSplit[i] );
    //  sort assembly sets
    AS.Sort ( crystSplit,nCrystSplits );

    //  recalculate assembly types (unique ids)

    serNo   = 1;
    nCSRes = 0;
    for (i=0;i<nCrystSplits;i++)  {
      if (crystSplit[i]->Score<8)  nCSRes++;
      for (j=0;j<crystSplit[i]->nAssemblies;j++)  {
        crystSplit[i]->A[j]->type  = -(crystSplit[i]->A[j]->type+1);
        crystSplit[i]->A[j]->serNo = serNo++;
      }
    }

    if (Complex)
      Complex->type = -(Complex->type+1);

    t = 0;
    for (i=0;i<nCrystSplits;i++)
      for (j=0;j<crystSplit[i]->nAssemblies;j++)  {
        tij = crystSplit[i]->A[j]->type;
        if (tij<0)  {
          for (k=0;k<nCrystSplits;k++)
            for (m=0;m<crystSplit[k]->nAssemblies;m++)
              if (crystSplit[k]->A[m]->type==tij)
                crystSplit[k]->A[m]->type = t;
          if (Complex)  {
            if (Complex->type==tij)
              Complex->type = t;
          }
          t++;
        }
      }

  }

  int Assemblies::getMaxAsmSize()  {
  int i,m;
    m = 0;
    for (i=0;i<nCrystSplits;i++)
      m = mmdb::IMax ( m,crystSplit[i]->getMaxAsmSize() );
    return m;
  }

  int Assemblies::getMaxAsmSize ( int crystSplitNo )  {
    if ((0<=crystSplitNo) && (crystSplitNo<nCrystSplits))
         return crystSplit[crystSplitNo]->getMaxAsmSize();
    else return 0;
  }

  int Assemblies::getNofAssembliesInSplits()  {
  int i,n;
    n = 0;
    for (i=0;i<nCSRes;i++)
      n += crystSplit[i]->nAssemblies;
    return n;
  }

  void Assemblies::calcInterfaceScores ( PPInterface interface )  {
  //   Output: interface[i]->css - complexation significance scores.
  //
  //   Interfaces are scored as
  //
  //       PBR = max(SEi/SE)       for stable assemblies
  //       PBR = max(SEi/SE)/10.0  for assemblies in "grey zone"
  //
  //   SEi - sum of stabilization energy contributed by interfaces
  //         of type i in an assembly
  //   SE  - sum of all dSEi in an assembly
  //
  mmdb::rvector  itSE;
  mmdb::ivector  intf;
  mmdb::realtype SE,deltaG,f,intf_se;
  int      i,j,k;

    mmdb::GetVectorMemory ( itSE,nInterfaces,1 );
    mmdb::GetVectorMemory ( intf,nInterfaces,0 );
    for (i=0;i<nInterfaces;i++)
      interface[i]->css = 0.0;

    for (i=0;i<nCrystSplits;i++)  {
      for (j=0;j<crystSplit[i]->nAssemblies;j++)  {
        deltaG = -crystSplit[i]->A[j]->freeEn;
        if (deltaG<0.0)      f = 1.0;
        else if (deltaG<2.0) f = 0.1;
                        else f = 0.0;
        if (f>0.0)  {
          crystSplit[i]->A[j]->getEngagedInterfaces ( intf,nInterfaces );
          SE = 0.0;
          for (k=1;k<=nInterfaces;k++)
            itSE[k] = 0.0;
          for (k=0;k<nInterfaces;k++)  {
            intf_se = intf[k]*interface[k]->stabEn;
            SE     += intf_se;
            itSE[interface[k]->type] += intf_se;
          }
          if (SE<0.0)
            for (k=0;k<nInterfaces;k++)
              interface[k]->css = mmdb::RMin ( 1.0,
                                    mmdb::RMax ( interface[k]->css,
                                      f*itSE[interface[k]->type]/SE ) );
        }
      }
    }

    mmdb::FreeVectorMemory ( itSE,1 );
    mmdb::FreeVectorMemory ( intf,0 );

  }


  void Assemblies::makeChainMapping()  {
  int i,j;
    for (i=0;i<nCrystSplits;i++)
      for (j=0;j<crystSplit[i]->nAssemblies;j++)
        crystSplit[i]->A[j]->makeChainMapping();
  }


  void Assemblies::makeOrientations ( PDomains       D,
                                      mmdb::PManager MMDB )  {
  int i,j;
    for (i=0;i<nCrystSplits;i++)
      for (j=0;j<crystSplit[i]->nAssemblies;j++)
        crystSplit[i]->A[j]->makeOrientation ( D,MMDB );
  }

  void Assemblies::checkOriginalOrientations ( PDomains D ) {
  mmdb::ivector icnt;
  int     i;
    mmdb::GetVectorMemory ( icnt,D->nDomains,0 );
    orig_chains = true; // true if all ASU chains are used
    for (i=0;i<nCrystSplits;i++)  {
      crystSplit[i]->checkOriginalOrientations ( D,icnt );
      if (!crystSplit[i]->orig_chains)
        orig_chains = false;
    }
    mmdb::FreeVectorMemory ( icnt,0 );
  }


  void Assemblies::Orth2Frac ( mmdb::PManager MMDB )  {
  int i,j;
    for (i=0;i<nCrystSplits;i++)
      for (j=0;j<crystSplit[i]->nAssemblies;j++)
        crystSplit[i]->A[j]->Orth2Frac ( MMDB );
  }


  void Assemblies::assignRCSBSymOps ( PRCSBData rcsbData )  {
  int i,j;
    if (rcsbData->nSymOps>0)  {
      rcsb_symops = true;
      for (i=0;i<nCrystSplits;i++)
        for (j=0;j<crystSplit[i]->nAssemblies;j++)
          crystSplit[i]->A[j]->assignRCSBSymOps ( rcsbData );
    } else
      rcsb_symops = false;
  }


  mmdb::xml::PXMLObject Assemblies::getAssembliesXML (
                                              mmdb::cpstr name,
                                              PDomains    D,
                                              PInterfaces PI,
                                              int         nCellOut )  {
  mmdb::xml::PXMLObject xml;
  int         i;

    xml = new mmdb::xml::XMLObject ( xml_pisa_results );
    if (name)  xml->AddObject ( new mmdb::xml::XMLObject(xml_name,name) );
         else  xml->AddObject ( new mmdb::xml::XMLObject(xml_name,"*") );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_status,"Ok") );
    xml->AddObject ( new mmdb::xml::XMLObject(xml_total_asm,nCSRes) );

    if (nCSRes<=0)
      xml->AddObject ( new mmdb::xml::XMLObject(xml_assessment,"Unstable") );
    else if (crystSplit[0]->Score>2)
      xml->AddObject ( new mmdb::xml::XMLObject(xml_assessment,"Grey") );
    else
      xml->AddObject ( new mmdb::xml::XMLObject(xml_assessment,"Stable") );

    i = -1;
    if (nCSRes>0)  {
      if ((crystSplit[0]->equiv_all) && (crystSplit[0]->stable_all))  {
        i = crystSplit[0]->A[0]->mmSize;
      }
    } else  {
      i = 1;
    }
    xml->AddObject ( new mmdb::xml::XMLObject(xml_mult_state,i) );

    if (orig_chains)
          xml->AddObject ( new mmdb::xml::XMLObject(xml_orig_chains,"Yes") );
    else  xml->AddObject ( new mmdb::xml::XMLObject(xml_orig_chains,"No" ) );

    for (i=0;i<nCSRes;i++)
      xml->AddObject ( crystSplit[i]->getCrystSplitXML(i+1,D,PI,nCellOut) );

    return xml;

  }

  static int maxNofCrystSplits = 50;

  void SetMaxNofCrystSplits ( int maxnCrystSplits )  {
    maxNofCrystSplits = maxnCrystSplits;
  }

  void Assemblies::TrimNofSets ( mmdb::pstr & warning )  {
  char S[200];
  int  i,n0;

    if (nCrystSplits>maxNofCrystSplits)  {
      for (i=maxNofCrystSplits;i<nCrystSplits;i++)
        if (crystSplit[i])  {
          delete crystSplit[i];
          crystSplit[i] = NULL;
        }
      n0       = nCrystSplits;
      nCrystSplits = maxNofCrystSplits;
      if (nCSRes>nCrystSplits)  {
        nCSRes = nCrystSplits;
        sprintf ( S,"First %i assembly sets of total %i have been "
                    "left in the solution.\n",maxNofCrystSplits,n0 );
        mmdb::CreateConcat ( warning,S );
      }
    }

  }


  void Assemblies::calcAsmStock ( PDomains D, mmdb::PManager MMDB )  {
  // This function assumes that all assembly sets are prepared and
  // that all assembly types are assigned
  mmdb::ivector  count;
  mmdb::realtype maxConc,ucVol;
  int            i,j,nTypes;

    if (asmStock)  {
      delete asmStock;
      asmStock = NULL;
    }

    nTypes = -1;
    for (i=0;i<nCrystSplits;i++)
      for (j=0;j<crystSplit[i]->nAssemblies;j++)
        nTypes = mmdb::IMax ( nTypes,crystSplit[i]->A[j]->type );

    nTypes++;
    if (nTypes>0)  {

      mmdb::GetVectorMemory ( count,nTypes,0 );
      for (i=0;i<nTypes;i++)
        count[i] = 0;

      asmStock = new AsmStock ( nTypes );
      for (i=0;i<nCrystSplits;i++)
        for (j=0;j<crystSplit[i]->nAssemblies;j++)
          asmStock->addToStock ( crystSplit[i]->A[j],count );


      maxConc = 1.0;

      ucVol = MMDB->GetCrystData()->Vol;
      if (ucVol>10.0)
        maxConc = 1.0e27/ucVol/mmdb::NAvogadro*MMDB->GetNumberOfSymOps();

      asmStock->processStock ( D,1.0e-12*maxConc,maxConc );

      mmdb::FreeVectorMemory ( count,0 );

    }

  }


  void Assemblies::write ( mmdb::io::RFile f )  {
  int        i;
  mmdb::byte Version;

    Version = 4;
    f.WriteByte ( &Version     );

    f.WriteInt  ( &nCrystSplits    );
    f.WriteInt  ( &nCSRes     );
    f.WriteInt  ( &nInterfaces );
    for (i=0;i<nCrystSplits;i++)
      crystSplit[i]->write ( f,nInterfaces );
    f.WriteBool ( &rcsb_symops );
    f.WriteBool ( &orig_chains );

    if (asmStock)  {
      i = 1;
      f.WriteInt ( &i );
      asmStock->write ( f,nInterfaces );
    } else  {
      i = 0;
      f.WriteInt ( &i );
    }

  }

  void Assemblies::read ( mmdb::io::RFile f )  {
  int        i;
  mmdb::byte Version;

    FreeMemory ();

    f.ReadByte ( &Version     );

    f.ReadInt  ( &nCrystSplits    );
    f.ReadInt  ( &nCSRes     );
    f.ReadInt  ( &nInterfaces );
    if (nCrystSplits>0)  {
      crystSplit = new PCrystSplit[nCrystSplits];
      for (i=0;i<nCrystSplits;i++)  {
        crystSplit[i] = new CrystSplit();
        crystSplit[i]->read ( f,nInterfaces );
      }
    }
    if (Version>1)
         f.ReadBool ( &rcsb_symops );
    else rcsb_symops = false;
    if (Version>2)
         f.ReadBool ( &orig_chains );
    else orig_chains = true;


    if (Version>3)  {
      f.ReadInt ( &i );
      if (i>0)  {
        asmStock = new AsmStock(0);
        asmStock->read ( f,nInterfaces );
      }
    }

  }

  MakeStreamFunctions(Assemblies)


  mmdb::cpstr getAsmStatus ( ASSMB_RC asmRC )  {

    switch (asmRC)  {
      case ASSMB_Ok             : return "Ok";
      case ASSMB_Void           : return "Not calculated";
      case ASSMB_incompleteData : return "Incomplete data";
      case ASSMB_noSymOps       : return "No symmetry operations";
      case ASSMB_noDomains      : return "No domains";
      case ASSMB_noInterfaces   : return "No interfaces";
      case ASSMB_Overlap        : return "Overlapping structures";
      case ASSMB_noSymOp        : return "No symmetry operation";
      case ASSMB_improperSymOp  : return "Improper symop";
      case ASSMB_brokenComposition     :
                           return "Broken composition in PA graph";
      case ASSMB_brokenComplementarity :
                           return "Broken complementarity in PA graph";
      case ASSMB_repeatedAssignment    :
                                  return "Repeated interface assignment";
      case ASSMB_assemblyOverflow      :
                                  return "Assembly overflow";
      case ASSMB_tooBigSystem   : return "Too big system";
      case ASSMB_timeLimit      : return "Calculations over time limit";
      default                   : return "Undocumented error";
    }

  }

}  // namespace pisa
