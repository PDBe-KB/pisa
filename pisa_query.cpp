// $Id: pisa_query.cpp $
// =================================================================
//
//    05.01.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// -----------------------------------------------------------------
//
//  **** Module  :  PISA_Query <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::QueryData
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2014
//
// =================================================================
//

#include <string.h>
#include <math.h>
#include <sys/stat.h>

#if defined(_MSC_VER) || defined (WIN32)
#include <direct.h>
#endif

#include "pisa_query.h"
#include "pisa_engine.h"
#include "pisa_types.h"

namespace pisa  {

  //  =========================  QueryData  =========================

  QueryData::QueryData() : mmdb::io::Stream()  {
    InitQueryData();
  }

  QueryData::QueryData ( mmdb::io::RPStream Object )
           : mmdb::io::Stream(Object)  {
    InitQueryData();
  }

  QueryData::~QueryData()  {

    FreeMemory();

    if (ligExcl)
      delete[] ligExcl;

  }

  void QueryData::InitQueryData()  {

    dataStatus = DS_notAnalysed;
    nSymOps    = 0;
    nNCSOps    = 0;
    coorFile   = NULL;

    resolution  = -1.0;        // not given by default
    cell_a      = 0.0;         // cell parameters
    cell_b      = 0.0;
    cell_c      = 0.0;
    cell_alpha  = 0.0;
    cell_beta   = 0.0;
    cell_gamma  = 0.0;
    OrthCode    = 0;           // orthogonalization code
    strcpy ( spaceGroup,"" );

    title        = NULL;
    warnUNKRes   = NULL;    // warinings: unknown residues
    warnExcluded = NULL;    // warinings: excluded residues
    warnOverlap  = NULL;    // warinings: interface overlap
    warnAssembly = NULL;    // warinings from assembly building

    MMDB      = NULL;
    domains   = NULL;
    PI        = NULL;
    A         = NULL;
    Complex   = NULL;
    asmStatus = ASSMB_Void;
    cmpStatus = ASSMB_Void;

    ligExcl   = NULL;
    ligKey    = LIGANDS_Auto;

    nInterfaces = 0;
    nDomains    = 0;
    nNCSParents = 0;
    nAssemblies = 0;

  }

  void  QueryData::FreeMemory()  {
    DeleteWarnings();
    if (MMDB)     delete   MMDB;
    if (domains)  delete   domains;
    if (PI)       delete   PI;
    if (A)        delete   A;
    if (Complex)  delete   Complex;
    if (coorFile) delete[] coorFile;
    if (title)    delete[] title;
    MMDB      = NULL;
    domains   = NULL;
    PI        = NULL;
    A         = NULL;
    Complex   = NULL;
    asmStatus = ASSMB_Void;
    cmpStatus = ASSMB_Void;
    coorFile  = NULL;
    title     = NULL;
  }

  void  QueryData::DeleteWarnings()  {
    if (warnUNKRes)  {
      delete[] warnUNKRes;
      warnUNKRes = NULL;
    }
    if (warnExcluded)  {
      delete[] warnExcluded;
      warnExcluded = NULL;
    }
    if (warnOverlap)  {
      delete[] warnOverlap;
      warnExcluded = NULL;
    }
    if (warnAssembly)  {
      delete[] warnAssembly;
      warnAssembly = NULL;
    }
  }


  void QueryData::setLigExclude ( mmdb::cpstr ligList )  {
    if (ligList)  {
      mmdb::CreateCopCat ( ligExcl,",",ligList,"," );
      mmdb::DelSpaces ( ligExcl );
    } else if (ligExcl)  {
      delete[] ligExcl;
      ligExcl = NULL;
    }
  }

  void QueryData::setLigKey ( LIGAND_KEY ligProc )  {
    ligKey = ligProc;
  }

  int QueryData::getStructure ( mmdb::io::RFile f,
                                mmdb::cpstr fileSpec )  {
  mmdb::PModel model;
  int          i,n;

    dataStatus = DS_Ok;
    mmdb::CreateCopy ( coorFile,fileSpec );
    if (!coorFile)          dataStatus = DS_notAnalysed;
    else if (!coorFile[0])  dataStatus = DS_notAnalysed;

    if (dataStatus!=DS_Ok)  return dataStatus;

    if (!MMDB)  MMDB = new mmdb::Manager();
    SetMMDBReadFlags ( MMDB );

    readErr = MMDB->ReadCoorFile ( f );

    if (readErr==mmdb::Error_NoError)  {

      MMDB->GetStructureTitle ( title );
      resolution = MMDB->GetResolution();

      // leave only 1st model
      n = MMDB->GetNumberOfModels();
      for (i=n;i>1;i--)  {
        model = MMDB->GetModel(i);
        if (model)  delete model;
      }
      MMDB->DeleteAltLocs   ();
      MMDB->FinishStructEdit();
      MMDB->PDBCleanup ( mmdb::PDBCLEAN_SERIAL );

      if (MMDB->GetNumberOfAtoms()<=0)
        dataStatus |= DS_noAtoms;

    } else
      dataStatus = DS_badFile;

    return dataStatus;

  }

  DATASTATUS QueryData::analyseStructure ( mmdb::cpstr  fileName,
                                           PMolRefIndex molRef,
                                           mmdb::cpstr  agentsRef )  {

  mmdb::io::File f;

    if (fileName)  {
      f.assign ( fileName,false,false,mmdb::io::GZM_CHECK );
      if (f.reset(true))  {
        analyseStructure ( f,fileName,molRef,agentsRef );
        f.shut();
      } else
        dataStatus = DS_notAnalysed;
    } else
      analyseStructure ( f,NULL,molRef,agentsRef );

    return dataStatus;

  }

  DATASTATUS QueryData::analyseStructure ( mmdb::io::RFile f,
                                           mmdb::cpstr     fileSpec,
                                           PMolRefIndex    molRef,
                                           mmdb::cpstr     agentsRef ) {
  mmdb::mat44    TMatrix;
  mmdb::realtype vol;
  int            nAADomains,nNADomains,nLigands, i;

    if (fileSpec)  {
      // special mode, when NULL assume a repeat run on the same data,
      // all stored in mmdb

      mmdb::CreateCopy ( title,"-- not given --" );
      resolution = -1.0;
      cell_a     = 0.0;
      cell_b     = 0.0;
      cell_c     = 0.0;
      cell_alpha = 0.0;
      cell_beta  = 0.0;
      cell_gamma = 0.0;
      OrthCode   = 0;
      strcpy ( spaceGroup,"" );

      nDomains    = 0;
      nNCSParents = 0;
      nInterfaces = 0;
      nAssemblies = 0;

      getStructure  ( f,fileSpec );

    }

    DeleteWarnings();

    if (dataStatus==DS_Ok)  {

      if (MMDB->isSpaceGroup())
           strcpy ( spaceGroup,MMDB->GetSpaceGroupFix() );
      else dataStatus |= DS_noSpaceGroup;

      OrthCode = -1;
      if (MMDB->isCrystInfo())
        MMDB->GetCell ( cell_a,cell_b,cell_c,
                        cell_alpha,cell_beta,cell_gamma,
                        vol,OrthCode );
      else if (dataStatus & DS_noSpaceGroup)
            dataStatus |= DS_noCrystData;
      else  dataStatus |= DS_noCellData;

      nSymOps = MMDB->GetNumberOfSymOps();
      for (i=0;i<nSymOps;i++)
        if (MMDB->GetTMatrix(TMatrix,i,0,0,0)!=mmdb::SYMOP_Ok)
          dataStatus |= DS_noSymOp;  // no symmetry operation,
                                     // although it should be there
        else if (!isMat4Rot(TMatrix,1.0e-2))
          dataStatus |= DS_badSymOp; // bad symop or orthogonalisation

      nNCSOps = MMDB->GetNumberOfNCSMatrices();

      if (MMDB->CrystReady() & mmdb::CRRDY_NoOrthCode)
        OrthCode = -2;

      // mmdb::Create default list of domains, that is, each chain
      // is a domain
      domains = new Domains();
      domains->MakeDomains ( MMDB,molRef,agentsRef,
                             warnUNKRes,warnExcluded );
      nDomains    = domains->nDomains;
      nNCSParents = domains->nNCSParents;

      if (nDomains<=0)
        dataStatus |= DS_noChains;

      domains->getNofDomains ( nAADomains,nNADomains,nLigands );
      if (nAADomains+nNADomains<=0)
        dataStatus |= DS_noMacromolecules;

    }

    /*
    printf ( " Space Group: '%s'\n"
             " N. o symops: %i\n"
             " NCS ops    : %i\n",
             spaceGroup,nSymOps,nNCSOps );
  */

    return dataStatus;

  }


  bool QueryData::areWarnings()  {
    return (warnUNKRes || warnOverlap || warnExcluded || warnAssembly);
  }

  mmdb::pstr QueryData::getDomainName ( mmdb::pstr & S, int domainNo )  {
  int k;
    if (domains)  {
      mmdb::CreateCopCat ( S,mmdb::io::GetFName(coorFile),":",
                           domains->domain[domainNo]->range );
      k = strlen(S)-1;
      if (S[k]==':')  S[k] = char(0);
    } else
      mmdb::CreateCopy ( S,"*null*" );
    return S;
  }

  mmdb::pstr QueryData::getInterfaceName ( mmdb::pstr & S, int interfaceNo )  {
  PInterface interface;
  char       N[50],sep[5];
  int        k;
    if (PI && domains)  {
      if (interfaceNo<PI->getNofInterfaces())  {
        interface = PI->getInterface ( interfaceNo );
        sprintf ( N,"[%i]",interfaceNo+1 );
        k = strlen(domains->domain[interface->domain1]->range)-1;
        if (domains->domain[interface->domain1]->range[k]==':')
             sep[0] = char(0);
        else strcpy ( sep,":" );
        mmdb::CreateCopCat ( S,mmdb::io::GetFName(coorFile),N,
                         domains->domain[interface->domain1]->range,
                         sep,
                         domains->domain[interface->domain2]->range );
        k = strlen(S)-1;
        if (S[k]==':')  S[k] = char(0);
      } else
        mmdb::CreateCopy ( S,"*error*" );
    } else
      mmdb::CreateCopy ( S,"*null*" );
    return S;
  }


  mmdb::pstr QueryData::getInterfaceName1 ( mmdb::pstr & S, int interfaceNo )  {
  PInterface interface;
    if (PI && domains)  {
      if (interfaceNo<PI->getNofInterfaces())  {
        interface = PI->getInterface ( interfaceNo );
        mmdb::CreateCopCat ( S,
                         domains->domain[interface->domain1]->range,
                         " || ",
                         domains->domain[interface->domain2]->range );
      } else
        mmdb::CreateCopy ( S,"*error*" );
    } else
      mmdb::CreateCopy ( S,"*null*" );
    return S;
  }

  mmdb::pstr QueryData::getDomainRange ( mmdb::pstr S, int domainNo,
                                         int fieldLen )  {
    if (domains)  return domains->getDomainRange ( S,domainNo,fieldLen );
            else  strcpy ( S,"*null" );
    return S;
  }

  int QueryData::getDomainClass ( int domainNo )  {
    if (domains)  return domains->getDomainClass ( domainNo );
            else  return DCLASS_None;
  }

  mmdb::pstr QueryData::getChainID ( mmdb::pstr S, int domainNo )  {
    if (domains)  return domains->getChainID ( S,domainNo );
            else  strcpy ( S,"*null" );
    return S;
  }

  mmdb::pstr QueryData::getDomainCode ( mmdb::pstr & S, int domainNo )  {
  int k;
    if (domains)  {
      mmdb::CreateCopCat ( S,mmdb::io::GetFName(coorFile),":",
                           domains->domain[domainNo]->range );
      k = strlen(S)-1;
      if (S[k]==':')  S[k] = char(0);
    } else
      mmdb::CreateCopy ( S,"*null*" );
    return S;
  }


  void QueryData::AlignSimilarDomains()  {
    if (domains && MMDB)
      domains->AlignSimilar ( MMDB );
  }

  void QueryData::SortDomains()  {
    if (domains)
      domains->SortDomains();
  }

  void QueryData::makeSSGraphs()  {
  mmdb::ivector selHnd;
  int           i;
    if (domains && MMDB)  {
      selHnd = NULL;
      domains->makeSSGraphs ( MMDB,selHnd );
      if (selHnd)  {
        for (i=0;i<nNCSParents;i++)
          if (selHnd[i]>0)  MMDB->DeleteSelection ( selHnd[i] );
        mmdb::FreeVectorMemory ( selHnd,0 );
      }
    }
  }

  void QueryData::AlignAllDomains ( int startNo,int endNo )  {
    if (domains && MMDB)
      domains->AlignAll ( MMDB,startNo,endNo );
  }

  PROSURF_RC QueryData::CalcSurfaces ( PProSurf     proSurf,
                                        PMolRefIndex molRef )  {
    return domains->CalcSurfaces ( MMDB,proSurf,molRef );
  }

  #ifdef __debug
  extern void out4 ( mmdb::mat44 & T, mmdb::pstr name );
  #endif

  void QueryData::calcDimensions()  {
  // calculate reference frames
    if (domains && MMDB)
      domains->calcDimensions ( MMDB );
  }

  RESULT_CODE QueryData::CalcContacts ( int nCellLayers )  {
  //   Alignments of similar domains and domain dimension calculations
  // must be done before calling this function
  mmdb::mat44 TMatrix;
  TFrame      frame;
  int         i,j,k;

    if ((!MMDB) || (!domains))
      return RESULT_notReadyForContacts;

    if (MMDB->CrystReady()>=0)  {
      nSymOps = MMDB->GetNumberOfSymOps();
      nNCSOps = MMDB->GetNumberOfNCSMatrices();
      if (nSymOps<=0)
        return RESULT_noSymOps; // no symmetry operations, although they
                                // should be there
    } else  {
      nSymOps = 0;
      nNCSOps = 0;
    }

    // now make sure that reference frames are properly orientated
    // in all similar domains (by enforcing where possible)
    for (i=1;i<=domains->nDTypes;i++)  {
      k = -1;
      for (j=0;(j<domains->nNCSParents) && (k<0);j++)
        if (domains->domain[j]->type==i)  k = j;
      if (k>=0)
        for (j=k+1;j<domains->nNCSParents;j++)
          if ((domains->domain[j]->type==i) &&
              domains->isAligned(k,j))  {
            domains->getTMatrix ( TMatrix,k,j );
            domains->domain[k]->getReferenceFrame ( frame,TMatrix );
            domains->domain[j]->setReferenceFrame ( frame );
          }
    }

    for (i=domains->nNCSParents;i<domains->nDomains;i++)
      domains->domain[i]->CopyParentData ( domains->domain );

    if (!PI)  PI = new Interfaces();
        else  PI->Reset();

    // find all potential contacts. We calculate all potential
    // contacts first and only then run interface calculations,
    // because the contact procedure filters out equivalent
    // contacts.

    return PI->CalcContacts ( MMDB ,domains,nCellLayers );

  }

  /*
  int QueryData::CalcInterfaces ( PCProSurf       ProSurf,
                                   PCMolRefIndex   MolRef,
                                   PCCP4SRSManager SRS,
                                   int             nCellLayers )  {
  PInterface Interface;
  int         rc,i,m;

    nInterfaces = 0;

    MMDB->DeleteAllSelections();

    // calculate domain surfaces and dimensions for
    // identification of potential interfaces
    rc = domains->CalcSurfaces ( MMDB,ProSurf,MolRef );
    if (rc)  {
      if (rc<0)  return RESULT_SurfaceCalcsFailed;
      return rc;
    }

    // find all potential contacts. We calculate all potential
    // contacts first and only then run interface calculations,
    // because the contact procedure filters out equivalent
    // contacts.

    rc = CalcContacts ( nCellLayers );

    if (!rc)  {

      m = PI->getNofInterfaces();
      for (i=0;(i<m) && (!rc);i++)  {
        Interface = PI->getInterface ( i );
        if (Interface)  // calculate the interface
          rc = Interface->calcInterface (
                                MMDB ,domains->domain[Interface->domain1],
                                domains->domain[Interface->domain2],
                                ProSurf,MolRef,SRS );
      }

      if (rc!=PROSURF_Ok)  return rc;

      SortInterfaces();

    }

    return rc;

  }
  */

  void QueryData::checkResidueBSA()  {
  //   After all surface and interface calculations are done, this
  // function checks that total BSA never exceeds ASA for each and
  // every residue. These sorts of errors may happen due to finite
  // precision of surface and interface calculations. If a discreapancy
  // is found, the function corrects residue BSAs by applying a
  // common factor to interfacing residues in all interfaces, and
  // solvation energy effect is addjusted correspondingly.
  PPInterface    interface;
  PPDomain       domain;
  mmdb::rvector  si;
  mmdb::realtype bsa,f;
  int            i,j,k, nRes;
  bool           corrected;

    if ((!PI) || (!domains))  return;

    nInterfaces = PI->getNofInterfaces();
    nDomains    = domains->nDomains;
    if (nInterfaces<=0) return;
    if (nDomains<=0)    return;

    interface = PI->getInterfaces();
    domain    = domains->domain;

    // find symmetrical interfaces and factor them with 1/2
    mmdb::GetVectorMemory ( si,nInterfaces,0 );
    for (i=0;i<nInterfaces;i++)  {
      si[i] = 1.0;
      if (interface[i]->domain1==interface[i]->domain2)  {
        k    = 0;
        nRes = domain[interface[i]->domain1]->nRes;
        for (j=0;(j<nRes) && (!k);j++)
          if (fabs(interface[i]->bsa1[j]-interface[i]->bsa2[j])>4.0)
            k = 1;
        if (!k)  si[i] = 0.5;
      }
    }

    //  identify the need for correction and make the correction
    corrected = false;
    for (i=0;i<nDomains;i++)
      for (j=0;j<domain[i]->nRes;j++)  {
        bsa = 0.0;
        for (k=0;k<nInterfaces;k++)  {
          if (interface[k]->domain1==i)
            bsa += si[k]*interface[k]->bsa1[j];
          if (interface[k]->domain2==i)
            bsa += si[k]*interface[k]->bsa2[j];
        }
        if (bsa>domain[i]->asa[j])  {
          // correction is needed
          corrected = true;
          f = domain[i]->asa[j]/bsa;
          for (k=0;k<nInterfaces;k++)  {
            if (interface[k]->domain1==i)  {
              interface[k]->bsa1[j]    *= f;
              interface[k]->solvEn1[j] *= f;
            }
            if (interface[k]->domain2==i)  {
              interface[k]->bsa2[j]    *= f;
              interface[k]->solvEn2[j] *= f;
            }
          }
        }
      }

    // recalculate interface properties
    if (corrected)
      for (i=0;i<nInterfaces;i++)  {
        interface[i]->intArea1   = 0.0;
        interface[i]->intArea2   = 0.0;
        interface[i]->intDeltaG1 = 0.0;
        interface[i]->intDeltaG2 = 0.0;
        nRes = domain[interface[i]->domain1]->nRes;
        for (j=0;j<nRes;j++)  {
          interface[i]->intArea1   += interface[i]->bsa1   [j];
          interface[i]->intDeltaG1 -= interface[i]->solvEn1[j];
        }
        nRes = domain[interface[i]->domain2]->nRes;
        for (j=0;j<nRes;j++)  {
          interface[i]->intArea2   += interface[i]->bsa2   [j];
          interface[i]->intDeltaG2 -= interface[i]->solvEn2[j];
        }
        interface[i]->intArea   = (interface[i]->intArea1 +
                                   interface[i]->intArea2)/2.0;
        interface[i]->intDeltaG = interface[i]->intDeltaG1 +
                                  interface[i]->intDeltaG2;
      }

    mmdb::FreeVectorMemory ( si,0 );

  }

  void QueryData::SortInterfaces()  {
    if (PI && domains)  {
      PI->DeleteDummyInterfaces();
      PI->CalcInterfaceTypes   ( domains );
      PI->Sort                 ( ISORT_Area );
      nInterfaces = PI->getNofInterfaces();
    }
  }


  void QueryData::makeOverlapWarnings()  {
  PInterface Interface;
  char        N[100];
  int         i;

    if (warnOverlap)  {
      delete[] warnOverlap;
      warnOverlap = NULL;
    }

    nInterfaces = PI->getNofInterfaces();

    for (i=0;i<nInterfaces;i++)  {
      Interface = PI->getInterface ( i );
      if (Interface)  {
        if (Interface->overlap)  {
          sprintf ( N," %i",i+1 );
          mmdb::CreateConcat ( warnOverlap,N );
        }
      }
    }


  }

  void QueryData::setLigandsAssemble ( bool On )  {
    if (domains)
      domains->setLigandsAssemble ( On );
  }

  void QueryData::excludeLigands()  {
    if (domains)
      domains->excludeLigands ( ligExcl );
  }


  ASSMB_RC QueryData::calcAssemblies()  {
  Assembler     assembler;
  PPInterface   interface;
  mmdb::ovector fixed,xrel;
  int           i;

    nAssemblies = 0;
    if (warnAssembly)  {
      delete[] warnAssembly;
      warnAssembly = NULL;
    }

    assembler.SetLigandKey ( ligKey );

    asmStatus = assembler.CalcAssemblies ( MMDB,domains,PI );

    if (asmStatus==ASSMB_Ok)  {
      assembler.GetAssemblies ( A );
      A->TrimNofSets ( warnAssembly );
      nAssemblies = A->nCSRes;
    } else if (A)  {
      delete A;
      A = NULL;
    }

    if (A)
      A->calcInterfaceScores ( PI->getInterfaces() );
    else  {
      nInterfaces = PI->getNofInterfaces();
      for (i=0;i<nInterfaces;i++)
        PI->getInterface(i)->css = 0.0;
    }

    mmdb::GetVectorMemory ( fixed,nInterfaces,0 );
    mmdb::GetVectorMemory ( xrel ,nInterfaces,0 );
    interface = PI->getInterfaces();

    for (i=0;i<nInterfaces;i++)  {
      fixed[i] = interface[i]->fixedLigand;
      xrel [i] = interface[i]->Xrel;
    }

    cmpStatus = assembler.AnalyseComplex ( MMDB,domains,PI );
    if (cmpStatus<ASSMB_Overlap)  {
      assembler.GetComplex ( Complex );
    } else if (Complex)  {
      delete Complex;
      Complex = NULL;
    }

    for (i=0;i<nInterfaces;i++)  {
      interface[i]->fixedLigand = fixed[i];
      interface[i]->Xrel        = xrel [i];
    }
    mmdb::FreeVectorMemory ( xrel ,0 );
    mmdb::FreeVectorMemory ( fixed,0 );

    if (!A)  A = new Assemblies();
    A->Sort ( Complex );  // a MUST call
    A->calcAsmStock ( domains,MMDB );

    return asmStatus;

  }

  RESULT_CODE QueryData::assignRCSBSymOps ( mmdb::pstr rcsbSymOpFile )  {
  PRCSBData   rcsbData;
  mmdb::pstr  spaceGroup;
  RESULT_CODE rc;

    rc = RESULT_Ok;

    if (rcsbSymOpFile && (A || PI) && MMDB)  {

      spaceGroup = MMDB->GetSpaceGroupFix();
      if (spaceGroup)  {

        rcsbData = NULL;

        if (!getRCSBData(rcsbSymOpFile,spaceGroup,rcsbData))  {
          if (A)   A ->assignRCSBSymOps ( rcsbData );
          if (PI)  PI->assignRCSBSymOps ( rcsbData );
        } else
          rc = RESULT_noRCSBData;

        if (rcsbData)  delete rcsbData;

      }

    } else
      rc = RESULT_noRCSBassignment;

    return rc;

  }

  RESULT_CODE QueryData::assignRCSBSymOps ( PPRCSBData rcsbData,
                                            int nRCSBData )  {
  mmdb::pstr  spaceGroup;
  int         i,k;
  RESULT_CODE rc;

    if (rcsbData && (A || PI) && MMDB)  {
      spaceGroup = MMDB->GetSpaceGroupFix();
      if (spaceGroup)  {
        k = -1;
        for (i=0;(i<nRCSBData) && (k<0);i++)
          if (!strcmp(spaceGroup,rcsbData[i]->spaceGroup))
            k = i;
        if (k>=0)  {
          if (A)  A ->assignRCSBSymOps ( rcsbData[k] );
          if (PI) PI->assignRCSBSymOps ( rcsbData[k] );
          rc = RESULT_Ok;
        } else
          rc = RESULT_noRCSBSpaceGroup;
      } else
        rc = RESULT_noSpaceGroup;
    } else
      rc = RESULT_noResults;

    return rc;

  }


  int QueryData::getSimilarInterface (
       PDomain DA, mmdb::mat44 & Ta, int domain1, mmdb::realtype rmsdA,
       PDomain DB, mmdb::mat44 & Tb, int domain2, mmdb::realtype rmsdB,
       mmdb::mat44 & Tm )  {
  mmdb::mat44 Ta1,Tb1;

    if (!PI)       return -10;
    if (!domains)  return -11;

    mmdb::Mat4Inverse ( Ta,Ta1 );
    mmdb::Mat4Inverse ( Tb,Tb1 );

    return PI->getSimilarInterface ( DA,Ta1,domain1,rmsdA,
                                     DB,Tb1,domain2,rmsdB,
                                     Tm ,domains );

  }



  int QueryData::getNofInterfaces()  {
    return nInterfaces;
  }

  int QueryData::getNofITypes()  {
    if (PI)  return PI->getNofITypes();
    return 0;
  }

  int QueryData::getNofDomains()  {
    return nDomains;
  }

  void QueryData::getNofDomains ( int & nAADomains,
                                  int & nNADomains,
                                  int & nLigands )  {
    if (domains)
      domains->getNofDomains ( nAADomains,nNADomains,nLigands );
    else  {
      nAADomains = -1;
      nNADomains = -1;
      nLigands   = -1;
    }
  }

  int QueryData::getNofNCSParents()  {
    return nNCSParents;
  }

  int QueryData::getNofAssemblies()  {
    return nAssemblies;
  }

  PInterface QueryData::getInterface ( int interfaceNo )  {
    if (PI)  return PI->getInterface ( interfaceNo );
    return NULL;
  }

  PPInterface QueryData::getInterfaces()  {
    if (PI)  return PI->getInterfaces();
    return NULL;
  }

  PDomain QueryData::getDomain ( int domainNo )  {
    if (domains)  return domains->getDomain ( domainNo );
    return NULL;
  }

  int QueryData::getNCSParent ( int domainNo )  {
    if (domains)  return domains->getNCSParent ( domainNo );
    return 0;
  }

  PPDomain QueryData::getDomains()  {
    if (domains)  return domains->domain;
    return NULL;
  }

  mmdb::pstr QueryData::getStructName ( mmdb::pstr & S )  {
    mmdb::CreateConcat ( S,mmdb::io::GetFName(coorFile) );
    return S;
  }

  mmdb::cpstr QueryData::getStructName()  {
    return mmdb::io::GetFName(coorFile);
  }

  void QueryData::calcDomainStats (
                          int        domainID,
                          mmdb::realtype & bsa, // buried surface in interfaces
                          mmdb::realtype & DeltaG,
                          int &      nHB,
                          int &      nSB,
                          int &      nDS,
                          int &      nAt, // no. of buried atoms
                          int &      nRes // no. of buried residues
                                   )  {
  //  Calculates domain-related properties in the whole crystal
  PPInterface interface;
  int         i;

    bsa    = 0.0;
    DeltaG = 0.0;
    nHB    = 0;
    nSB    = 0;
    nDS    = 0;
    nAt    = 0;
    nRes   = 0;

    if (PI)  {

      interface = PI->getInterfaces();

      for (i=0;i<nInterfaces;i++)  {
        if (interface[i]->domain1==domainID)  {
          bsa    += interface[i]->intArea;
          DeltaG += interface[i]->intDeltaG1;
          nHB    += interface[i]->nHBonds;
          nSB    += interface[i]->nSBridges;
          nDS    += interface[i]->nDSBonds;
          nAt    += interface[i]->nIntAtoms1;
          nRes   += interface[i]->nIntRes1;
        }
        if (interface[i]->domain2==domainID)  {
          bsa    += interface[i]->intArea;
          DeltaG += interface[i]->intDeltaG2;
          nHB    += interface[i]->nHBonds;
          nSB    += interface[i]->nSBridges;
          nDS    += interface[i]->nDSBonds;
          nAt    += interface[i]->nIntAtoms2;
          nRes   += interface[i]->nIntRes2;
        }
      }

    }

  }

  void QueryData::getExclLigandList ( mmdb::PResName & ligName,
                                      int & nLigNames )  {
    if (domains)
      domains->getExclLigandList ( ligName,nLigNames );
    else  {
      if (ligName)
        delete[] ligName;
      ligName   = NULL;
      nLigNames = 0;
    }
  }

  void QueryData::writePIData ( mmdb::io::RFile f )  {
  //   Adjust reading PI data in CPIMachine's structure search
  //  if any changes done to the code below.
  int        ds,as;
  mmdb::byte Version=2;
    f.WriteByte    ( &Version );
    if ((!title) && MMDB)
      MMDB->GetStructureTitle ( title );  // retain structure title
    ds = dataStatus;
    as = asmStatus;
    f.WriteInt     ( &ds          );
    f.WriteInt     ( &nSymOps     );
    f.WriteInt     ( &nNCSOps     );
    f.WriteInt     ( &as          );
    f.CreateWrite  ( title        );
    f.WriteTerLine ( spaceGroup   );
    f.CreateWrite  ( warnUNKRes   );
    f.CreateWrite  ( warnExcluded );
    f.CreateWrite  ( warnOverlap  );
    f.CreateWrite  ( warnAssembly );
    StreamWrite    ( f ,domains    );
    StreamWrite    ( f,PI         );
    nInterfaces = 0;
    nDomains    = 0;
    nNCSParents = 0;
    if (PI)
      nInterfaces = PI->getNofInterfaces();
    if (domains)  {
      nDomains    = domains->nDomains;
      nNCSParents = domains->nNCSParents;
    }
  }

  void QueryData::readPIData ( mmdb::io::RFile f )  {
  int        n;
  mmdb::byte Version;
    f.ReadByte   ( &Version     );
    f.ReadInt    ( &n );
    dataStatus = n;
    f.ReadInt    ( &nSymOps     );
    f.ReadInt    ( &nNCSOps     );
    f.ReadInt    ( &n  );
    asmStatus  = ASSMB_RC(n);
    f.CreateRead ( title        );
    if (Version>1)
      f.ReadTerLine ( spaceGroup );
    f.CreateRead ( warnUNKRes   );
    f.CreateRead ( warnExcluded );
    f.CreateRead ( warnOverlap  );
    f.CreateRead ( warnAssembly );
    StreamRead   ( f,domains    );
    StreamRead   ( f,PI         );
    nInterfaces = 0;
    nDomains    = 0;
    nNCSParents = 0;
    if (PI)
      nInterfaces = PI->getNofInterfaces();
    if (domains)  {
      nDomains    = domains->nDomains;
      nNCSParents = domains->nNCSParents;
    }
  }

  void shrinkMMDB ( mmdb::PManager MMDB )  {
  //  Removes unnecessary, auxiliary and temporary data from
  //  MMDB before writing the structure into a file.
  mmdb::PPAtom atom;
  int          uddHnd,selHnd,nAtoms,i,hflag;

    // Delete annotation
    MMDB->Delete ( mmdb::MMDBFCM_Title | mmdb::MMDBFCM_SecStruct |
                   mmdb::MMDBFCM_SA    | mmdb::MMDBFCM_SB        |
                   mmdb::MMDBFCM_SC    | mmdb::MMDBFCM_Footnotes );

    // 2nd model was used in interface calculations, but
    // now it's a waste
    MMDB->DeleteModel  ( 2 );

    MMDB->RemoveBonds();  // a must

    uddHnd = MMDB->GetUDDHandle ( mmdb::UDR_HIERARCHY,hydrogen_udd );
    if (uddHnd>0)  {
      MMDB->GetUDData ( uddHnd,hflag );
      if (hflag)  {
        // remove added hydrogens
        selHnd = MMDB->NewSelection();
        MMDB->Select ( selHnd,mmdb::STYPE_ATOM,0,"*",
                       mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                       "*","*","H","*",mmdb::SKEY_NEW );
        MMDB->GetSelIndex ( selHnd,atom,nAtoms );
        for (i=0;i<nAtoms;i++)
          if (atom[i]->sigOcc<0.0)  {
            delete atom[i];
            atom[i] = NULL;
          }
        MMDB->DeleteSelection ( selHnd );
        MMDB->PutUDData ( uddHnd,0 );  // flag "no added hydrogens"
      }
    }

    MMDB->FinishStructEdit();

  }

  void QueryData::writeStructure ( mmdb::io::RFile f )  {
  mmdb::byte Version=1;

    shrinkMMDB ( MMDB );

    f.WriteByte ( &Version );
    StreamWrite ( f,MMDB   );  // write out structure
    
  }

  void QueryData::readStructure ( mmdb::io::RFile f )  {
  mmdb::byte Version;
    f.ReadByte ( &Version );
    StreamRead ( f,MMDB   );
  }


  void QueryData::writeAssemblies ( mmdb::io::RFile f )  {
  int        k;
  mmdb::byte Version=1;
    f.WriteByte ( &Version );
    k = asmStatus;
    f.WriteInt  ( &k  );
    StreamWrite ( f,A );
    if (A)  nAssemblies = A->nCSRes;
      else  nAssemblies = 0;
    if (Complex)  {
      f.WriteInt  ( &nInterfaces );
      Complex->write ( f,nInterfaces );
    } else  {
      k = -1;
      f.WriteInt  ( &k );
    }
  }

  void QueryData::readAssemblies ( mmdb::io::RFile f )  {
  int        k;
  mmdb::byte Version;
    f.ReadByte ( &Version );
    f.ReadInt  ( &k );
    asmStatus = (ASSMB_RC)k;
    StreamRead ( f,A );
    if (A)  nAssemblies = A->nCSRes;
      else  nAssemblies = 0;
    f.ReadInt  ( &k );
    if (k>=0)  {
      nInterfaces = k;
      if (!Complex)  Complex = new Assembly();
      Complex->read ( f,k );
    } else if (Complex)  {
      delete Complex;
      Complex = NULL;
    }
  }


  void QueryData::outputSIData ( mmdb::io::RFile f, int nCellLayers )  {
  PInterface Interface;
  char       S[500];
  int        i,nInterfaces;

    if (domains)  {

      f.WriteLine (
    "\n Interfacing structure(s)\n"
    " -----------------------------------------------------------------------\n"
    "     |       |   Structure  |     Surface  |              |\n"
    "  ## | Range |  Nat  | Nres |  Nat  | Nres | Surface area |    Delta G\n"
    " ----+-------+-------+------+-------+------+--------------+-------------"
                );
      for (i=0;i<domains->nDomains;i++)  {
        sprintf ( S," %3i |   %2s  | %5i | %4i | %5i | %4i |"
                    "  %10.1f  | %10.1f",
                  i+1 ,domains->domain[i]->range,
                      domains->domain[i]->nAtoms,
                      domains->domain[i]->nRes,
                      domains->domain[i]->nSurfAtoms,
                      domains->domain[i]->nSurfRes,
                      domains->domain[i]->surfArea,
          domains->domain[i]->DeltaG );
        f.WriteLine ( S );
      }
      f.WriteLine (
    " -----------------------------------------------------------------------\n"
                  );

    }

    if (PI)  {

      f.WriteLine (
    " Interfaces (Part 1)\n"
    " -------------------------------------------------------------------\n"
    "  ## | Range |  Nat  | Nres | Surface area |   Delta G  | nHB | nSB \n"
    " ----+-------+-------+------+--------------+------------+-----+-----"
                );
      nInterfaces = PI->getNofInterfaces();
      for (i=0;i<nInterfaces;i++)  {
        Interface = PI->getInterface ( i );
        sprintf ( S," %3i | %2s %2s | %5i | %4i |  %10.1f  |"
                    " %10.1f | %3i | %3i ",
                  i+1,
           domains->domain[Interface->domain1]->range,
           domains->domain[Interface->domain2]->range,
           Interface->nIntAtoms1+Interface->nIntAtoms1,
           Interface->nIntRes1+Interface->nIntRes1,
           Interface->intArea,
           Interface->intDeltaG,
           Interface->nHBonds,
           Interface->nSBridges );
        f.WriteLine ( S );
      }

      f.WriteLine (
    " -------------------------------------------------------------------\n"
    "\n"
    " Interfaces (Part 2)\n"
    " -------------------------------------------------------------------\n"
    "  ## | Range | Cell |   Symmetry operation\n"
    " ----+-------+------+-----------------------------------------------"
                );
      for (i=0;i<nInterfaces;i++)  {
        Interface = PI->getInterface ( i );
        sprintf ( S," %3i | %2s %2s |  %1i%1i%1i | %s ",
                    i+1,
           domains->domain[Interface->domain1]->range,
           domains->domain[Interface->domain2]->range,
           Interface->cell_i+nCellLayers,
           Interface->cell_j+nCellLayers,
           Interface->cell_k+nCellLayers,
           Interface->symOp );
        f.WriteLine ( S );
      }
      f.WriteLine (
    " -------------------------------------------------------------------\n"
                  );

    }

  }

  void QueryData::outputAsmData ( mmdb::io::RFile f )  {
  PAssembly     Asm;
  char          S[500],N[50];
  mmdb::ivector intf;
  mmdb::pstr    F;
  int           i,j,k,nA;

    F    = NULL;
    intf = NULL;

    if (A)  {

      if (A->nCSRes<=0)
         f.WriteLine ( "\n"
                       "   STRUCTURES DO NOT COMPLEXATE IN COLUTION"
                       "\n\n\n" );
      else if (!A->crystSplit[0]->A[0]->stable)
         f.WriteLine ( "\n"
                       "   NO STABLE ASSEMBLIES ARE FOUND\n\n\n" );

      mmdb::GetVectorMemory ( intf,A->nInterfaces,0 );

      for (i=0;i<A->nCSRes;i++)  {

        sprintf ( S,"\n\n  PQS set #%i\n"
                    " ==============\n",i+1 );
        f.WriteLine ( S );

        nA = A->crystSplit[i]->nAssemblies;
        for (j=0;j<nA;j++)  {
          Asm = A->crystSplit[i]->A[j];
          if (Asm->stable)  strcpy ( N,"stable" );
                      else  strcpy ( N,"unstable" );
          sprintf ( S,"\n   Assembly #%i (%s)\n"
                      "   ~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
            "   Type Id: %i   Size: %i    No. in unit cell: %i\n"
            "   No. of dissociating parts   : %i\n"
            "   Solvent-accessible area     : %10.1f\n"
            "   Buried surface area         : %10.1f\n"
            "   Solvation energy gain       : %10.1f\n"
            "   Free energy of dissociation : %10.1f\n"
            "   Ground free energy          : %10.1f\n"
            "   Dissociation interface area : %10.1f\n",
            j+1,N,
             Asm->type+1,Asm->asmSize,Asm->nUC,Asm->nDiss,Asm->asa,
             Asm->bsa   ,Asm->seGain ,Asm->freeEn,Asm->freeEn0,
             Asm->dissIntArea );
          f.WriteLine ( S );
          f.Write ( "   Formula                 : " );
          f.WriteLine ( Asm->getFormula ( F ,domains,0,NULL ) );
          f.Write ( "   Composition             : " );
          f.WriteLine ( Asm->getComposition ( F ,domains,0,NULL ) );
          f.Write ( "   Dissociation            : " );
          f.WriteLine ( Asm->getDissPattern ( F ,domains,0,NULL ) );

          Asm->getEngagedInterfaces ( intf,A->nInterfaces );
          f.Write ( "   Engaged interfaces      : " );
          for (k=0;k<A->nInterfaces;k++)
            if (intf[k]>0)  {
              if (intf[k]>1)  sprintf ( S," %ix[%i]",intf[k],k+1 );
                        else  sprintf ( S," [%i]",k+1 );
              f.Write ( S );
            }
          f.LF();
          f.Write ( "   Dissociating interfaces : " );
          Asm->getDissInterfaces ( intf,A->nInterfaces );
          for (k=0;k<A->nInterfaces;k++)
            if (intf[k]>0)  {
              if (intf[k]>1)  sprintf ( S," %ix[%i]",intf[k],k+1 );
                        else  sprintf ( S," [%i]",k+1 );
              f.Write ( S );
            }
          f.LF();

        }

      }

      mmdb::FreeVectorMemory ( intf,0 );

    } else
      f.WriteLine ( "\n"
                    "   NO ASSEMBLIES HAS BEEN CALCULATED\n" );

    if (F)  delete[] F;

  }

  void QueryData::outputAssemblies ( mmdb::pstr outDir,
                                     int maxNofAsmSets )  {
  mmdb::PManager AsmStructure;
  mmdb::pstr     pqsDir,asmFName;
  char           N[50];
  int            i,j,na;

    pqsDir   = NULL;
    asmFName = NULL;
    if (A)  {
      na = mmdb::IMin ( A->nCSRes,maxNofAsmSets );
      for (i=0;i<na;i++)  {
        sprintf ( N,"%i",i+1 );
  #if defined(_MSC_VER) || defined (WIN32)
        mmdb::CreateCopCat ( pqsDir,outDir,"pqsset",N,"\\" );
        if (!_mkdir(pqsDir))  {
  #else
        mmdb::CreateCopCat ( pqsDir,outDir,"pqsset",N,"/" );
        if (!mkdir(pqsDir,06777))  {
  #endif
          for (j=0;j<A->crystSplit[i]->nAssemblies;j++)  {
            AsmStructure = A->crystSplit[i]->A[j]->getAssemblyStructure (
                                                       MMDB ,domains );
            if (AsmStructure)  {
              sprintf ( N,"%i",j+1 );
              mmdb::CreateCopCat ( asmFName,pqsDir,"assembly",N,".pdb" );
              AsmStructure->WritePDBASCII ( asmFName );
            }
          }
        }
      }
    }

  }

  mmdb::xml::PXMLObject QueryData::getInterfacesXML()  {
    if (PI && MMDB && domains)
          return PI->getXML ( MMDB ,domains->domain );
    else  return NULL;
  }

  void QueryData::read ( mmdb::io::RFile f )  {
  int  n,Version;

    FreeMemory();

    f.ReadInt     ( &Version     );
    f.ReadInt     ( &n           );
    dataStatus = n;
    f.ReadInt     ( &nSymOps     );
    f.ReadInt     ( &nNCSOps     );
    f.CreateRead  ( coorFile     );

    f.ReadReal    ( &resolution  );
    f.ReadReal    ( &cell_a      );
    f.ReadReal    ( &cell_b      );
    f.ReadReal    ( &cell_c      );
    f.ReadReal    ( &cell_alpha  );
    f.ReadReal    ( &cell_beta   );
    f.ReadReal    ( &cell_gamma  );
    f.ReadInt     ( &OrthCode    );
    f.ReadTerLine ( spaceGroup   );
    f.CreateRead  ( title        );
    f.CreateRead  ( warnUNKRes   );
    f.CreateRead  ( warnExcluded );
    f.CreateRead  ( warnOverlap  );
    f.CreateRead  ( warnAssembly );

    f.ReadInt     ( &nInterfaces );
    f.ReadInt     ( &nDomains    );
    f.ReadInt     ( &nNCSParents );
    f.ReadInt     ( &nAssemblies );

  }

  void  QueryData::write ( mmdb::io::RFile f )  {
  int n;
  int Version = 1;

    f.WriteInt     ( &Version     );
    n = dataStatus;
    f.WriteInt     ( &n           );
    f.WriteInt     ( &nSymOps     );
    f.WriteInt     ( &nNCSOps     );
    f.CreateWrite  ( coorFile     );

    f.WriteReal    ( &resolution  );
    f.WriteReal    ( &cell_a      );
    f.WriteReal    ( &cell_b      );
    f.WriteReal    ( &cell_c      );
    f.WriteReal    ( &cell_alpha  );
    f.WriteReal    ( &cell_beta   );
    f.WriteReal    ( &cell_gamma  );
    f.WriteInt     ( &OrthCode    );
    f.WriteTerLine ( spaceGroup   );
    f.CreateWrite  ( title        );
    f.CreateWrite  ( warnUNKRes   );
    f.CreateWrite  ( warnExcluded );
    f.CreateWrite  ( warnOverlap  );
    f.CreateWrite  ( warnAssembly );

    f.WriteInt     ( &nInterfaces );
    f.WriteInt     ( &nDomains    );
    f.WriteInt     ( &nNCSParents );
    f.WriteInt     ( &nAssemblies );

  }

  MakeStreamFunctions(QueryData)


  static mmdb::pstr sym_lib = NULL;

  void SetSyminfoLib ( mmdb::pstr syminfo_lib )  {
    sym_lib = syminfo_lib;
  }

  void SetMMDBReadFlags ( mmdb::PManager MMDB )  {
    MMDB->SetSyminfoLib ( sym_lib );
    MMDB->SetFlag ( mmdb::MMDBF_IgnoreDuplSeqNum       |
                    mmdb::MMDBF_IgnoreNonCoorPDBErrors |
                    mmdb::MMDBF_IgnoreBlankLines       |
                    mmdb::MMDBF_FixSpaceGroup );
  }


  bool CheckDataStatus ( DATASTATUS ds, RESULT_CODE & rc,
                         mmdb::pstr & msg )  {

    if (ds & pisa::DS_notAnalysed)   {
      msg = mmdb::pstr("structure was not analysed (failure)");
      rc  = pisa::RESULT_notAnalysed;
    } else if (ds & pisa::DS_noAtoms)  {
      msg = mmdb::pstr("no atoms found in the file\n");
      rc  = pisa::RESULT_noAtoms;
    } else if (ds & pisa::DS_noChains)  {
      msg = mmdb::pstr("no chains found in the file\n");
      rc  = pisa::RESULT_noCoordinateData;
    } else if (ds & pisa::DS_noMacromolecules)  {
      msg = mmdb::pstr("no macromolecules found in the file\n");
      rc  = pisa::RESULT_noMacromolecules;
    } else if (ds & pisa::DS_badFile)  {
      msg = mmdb::pstr("bad file (cannot parse)\n");
      rc  = pisa::RESULT_badCoordinateFile;
    } else  {
      msg = mmdb::pstr("");
      rc = pisa::RESULT_Ok;
    }

    return (rc!=pisa::RESULT_Ok);

  }

}  // namespace pisa
