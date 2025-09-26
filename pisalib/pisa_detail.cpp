// $Id: pisa_detail.cpp $
// =================================================================
//
//    04.03.16   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_lists <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Lists
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2016
//
// =================================================================
//

#include <string.h>
#include <math.h>

#include "pisa_detail.h"
#include "pisa_types.h"
#include "pisa_defs.h"


#include<iostream>
using namespace std;

namespace pisa  {

  // =======================  Detail  ==========================

  Detail::Detail ( mmdb::cpstr confPath ) : Lists ( confPath )  {
    InitDetail();
  }

  Detail::~Detail()  {
  }

  void Detail::InitDetail()  {
  }

  void Detail::set_As_Is_Key_xml ( AS_IS_KEY as_is_Proc )  {
    asisKey = as_is_Proc;
  }

  void Detail::printInterfaceSummary ( mmdb::io::RFile f,
                                       PInterface interface )  {
  PDomain        D1,D2;
  char           S[2000],name1[20],name2[20],dclass1[20],dclass2[20];
  char           symOp[200], SE1[20],SE2[20];
  mmdb::realtype nA1,nA2,nR1,nR2;
  int            symopNo,i,k;

    D1 = query->getDomain ( interface->domain1 );
    D2 = query->getDomain ( interface->domain2 );

    if (query->PI->rcsb_symops)  symopNo = interface->rcsb_symop;
                           else  symopNo = interface->symOpNo+1;

    nA1 = D1->nAtoms;  nA1 /= 100.0;
    nA2 = D2->nAtoms;  nA2 /= 100.0;
    nR1 = D1->nRes;    nR1 /= 100.0;
    nR2 = D2->nRes;    nR2 /= 100.0;

    strcpy ( symOp,interface->symOp );
    i = strlen(symOp);
    k = (20-i)/2;
    while (k>0)  {
      symOp[i++] = ' ';
      k--;
    }
    symOp[i] = char(0);

    if (D1->dclass==DCLASS_Protein)
          sprintf ( SE1,"%10.1f",D1->DeltaG );
    else  strcpy  ( SE1,"N/A" );

    if (D2->dclass==DCLASS_Protein)
          sprintf ( SE2,"%10.1f",D2->DeltaG );
    else  strcpy  ( SE2,"N/A" );


    sprintf ( S,
  " 1. Interface Summary\n"
  " --------------------------.-------------------.-------------------\n"
  "                           |    Structure 1    |    Structure 2 \n"
  " --------------------------+-------------------+-------------------\n"
  "           Selection range | %15s   | %15s \n"
  "                     Class |   %10s      |   %10s \n"
  "        Symmetry operation |       X,Y,Z       | %20s \n"
  "               Symmetry ID |       1_555       |      %i_%i%i%i \n"
  " --------------------------+-------------------+-------------------\n"
  "    Atoms in the interface | %8i (%5.1f%%) | %8i (%5.1f%%)\n"
  "            on the surface | %8i (%5.1f%%) | %8i (%5.1f%%)\n"
  "                     total | %8i (100.0%%) | %8i (100.0%%)\n"
  " --------------------------+-------------------+-------------------\n"
  " Residues in the interface | %8i (%5.1f%%) | %8i (%5.1f%%)\n"
  "            on the surface | %8i (%5.1f%%) | %8i (%5.1f%%)\n"
  "                     total | %8i (100.0%%) | %8i (100.0%%)\n"
  " --------------------------+-------------------+-------------------\n"
  "         Buried ASA, sq. A | %8.1f (%5.1f%%) | %8.1f (%5.1f%%)\n"
  "          Total ASA, sq. A | %8.1f (100.0%%) | %8.1f (100.0%%)\n"
  " --------------------------+-------------------+-------------------\n"
  " Solvation energy kcal/mol | %11s       | %11s \n"
  "         SE gain, kcal/mol | %11.1f       | %11.1f \n"
  " --------------------------'-------------------'-------------------\n"
  "\n",
      D1->getDomainRange(name1,11), D2->getDomainRange(name2,11),
      D1->getDomainClass(dclass1),  D2->getDomainClass(dclass2),
      symOp,symopNo,
      interface->cell_i+nCellOut,
      interface->cell_j+nCellOut,
      interface->cell_k+nCellOut,
      interface->nIntAtoms1 ,interface->nIntAtoms1/nA1,
      interface->nIntAtoms2 ,interface->nIntAtoms2/nA2,
      D1->nSurfAtoms,D1->nSurfAtoms/nA1,
      D2->nSurfAtoms,D2->nSurfAtoms/nA2,
      D1->nAtoms, D2->nAtoms,
      interface->nIntRes1 ,interface->nIntRes1/nR1,
      interface->nIntRes2 ,interface->nIntRes2/nR2,
      D1->nSurfRes,D1->nSurfRes/nR1,
      D2->nSurfRes,D2->nSurfRes/nR2,
      D1->nRes, D2->nRes,
      interface->intArea1,100.0*interface->intArea1/D1->surfArea,
      interface->intArea2,100.0*interface->intArea2/D2->surfArea,
      D1->surfArea,D2->surfArea,
      SE1,SE2,
      interface->intDeltaG1,interface->intDeltaG2
   );

    f.Write ( S );


  }

  mmdb::PAtom Detail::GetAtom ( int serNum )  {
  mmdb::PPAtom A;
  mmdb::PAtom  atom;
  int          nAtoms,i;

    if (!query->MMDB)  return NULL;

    query->MMDB->GetAtomTable ( A,nAtoms );

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

  void  Detail::printBondTable ( mmdb::io::RFile f,
                                 mmdb::cpstr     Title,
                                 mmdb::cpstr     noBondMsg,
                                 PInterface      interface,
                                 PIBond          Bond,
                                 int             nBonds )  {
  UNUSED_ARGUMENT(interface);
  mmdb::PAtom   atom;
  char          S[500],S1[100],S2[100];
  mmdb::ChainID chID;
  int           i;

    if (nBonds<=0)  {
      f.WriteLine ( noBondMsg );
      return;
    }

    sprintf ( S,
    " %s\n"
    " ----.------------------.-------.------------------\n"
    "  ## |    Structure 1   | Dist. |    Structure 2\n"
    " ----+------------------+-------+------------------\n",
      Title );
    f.Write ( S );

    for (i=0;i<nBonds;i++)  {
      atom = GetAtom ( Bond[i].serNum1 );
      if (atom)  {
        strcpy ( chID,atom->GetChainID() );
        if (!chID[0])  strcpy ( chID,"-" );
        sprintf ( S1,"%1s:%3s%4i%1s[%4s]",
                     chID,atom->GetResName(),
                     atom->GetSeqNum(),atom->GetInsCode(),
                     atom->name );
      } else
        strcpy ( S1,"                   " );
      atom = GetAtom ( Bond[i].serNum2 );
      if (atom)  {
        strcpy ( chID,atom->GetChainID() );
        if (!chID[0])  strcpy ( chID,"-" );
        sprintf ( S2,"%1s:%3s%4i%1s[%4s]",
                     chID,atom->GetResName(),
                     atom->GetSeqNum(),atom->GetInsCode(),
                     atom->name );
      } else
        strcpy ( S2,"                   " );

      sprintf ( S," %3i | %s |  %3.1f  | %s",i+1,S1,Bond[i].dist,S2 );

      sprintf ( S2," res:[%3i : %3i]",Bond[i].res1,Bond[i].res2 );
      strcat ( S,S2 );

      f.WriteLine ( S );
    }

    f.WriteLine (
       " ----'------------------'-------'------------------\n" );

  }


  void  Detail::printIntResTable ( mmdb::io::RFile f,
                                   mmdb::cpstr     Title,
                                   mmdb::cpstr     noResMsg,
                                   PInterface      interface,
                                   int             structNo )  {
  PDomain         D;
  mmdb::PManager  MMDB;
  mmdb::PPResidue Res;
  mmdb::PAtom     atom;
  mmdb::rvector   bsa,solvEn;
  char            S[500],rid[50],sb[50];
  mmdb::psvector  stab_s,destab_s;
  mmdb::ChainID   chID;
  PSREffect       sre_stab,sre_destab;
  int             nRes,nStab,nDestab;
  int             selHndRes,selHndResS,selHndResI;
  int             selHndResH,selHndResB;
  int             i,j;

    if (structNo<=1)  {
      D          = query->getDomain ( interface->domain1 );
      bsa        = interface->bsa1;
      solvEn     = interface->solvEn1;
      sre_stab   = interface->sre_stab1;
      sre_destab = interface->sre_destab1;
      nStab      = interface->nStab1;
      nDestab    = interface->nDestab1;
    } else  {
      D          = query->getDomain ( interface->domain2 );
      bsa        = interface->bsa2;
      solvEn     = interface->solvEn2;
      sre_stab   = interface->sre_stab2;
      sre_destab = interface->sre_destab2;
      nStab      = interface->nStab2;
      nDestab    = interface->nDestab2;
    }

    stab_s   = new mmdb::pstr[nStab];
    destab_s = new mmdb::pstr[nDestab];
    for (i=0;i<nStab;i++)    stab_s  [i] = NULL;
    for (i=0;i<nDestab;i++)  destab_s[i] = NULL;

    MMDB = query->MMDB;
    selHndRes  = MMDB->NewSelection();
    selHndResS = MMDB->NewSelection();
    selHndResI = MMDB->NewSelection();
    selHndResH = MMDB->NewSelection();
    selHndResB = MMDB->NewSelection();

    D->SelectDomain   ( selHndRes,MMDB,mmdb::STYPE_RESIDUE,
                                       mmdb::SKEY_NEW,1,false );
    MMDB->GetSelIndex ( selHndRes,Res,nRes );

    if (nRes<=0)  {
      f.WriteLine ( noResMsg );
      return;
    }

    MMDB->Select ( selHndResS,mmdb::STYPE_RESIDUE,
                   D->selHndSurf,mmdb::SKEY_NEW );

    if (structNo<=1)  i = interface->selHndInt1;
                else  i = interface->selHndInt2;
    MMDB->Select ( selHndResI,mmdb::STYPE_RESIDUE,i,mmdb::SKEY_NEW );

    for (i=0;i<interface->nHBonds;i++)  {
      if (structNo<=1)  atom = GetAtom(interface->HBond[i].serNum1);
                  else  atom = GetAtom(interface->HBond[i].serNum2);
      if (atom)
        MMDB->SelectResidue ( selHndResH,atom->GetResidue(),
                              mmdb::STYPE_RESIDUE,mmdb::SKEY_OR,false );
    }
    MMDB->MakeSelIndex ( selHndResH );

    for (i=0;i<interface->nSBridges;i++)  {
      if (structNo<=1)
            atom = GetAtom ( interface->SBridge[i].serNum1 );
      else  atom = GetAtom ( interface->SBridge[i].serNum2 );
      if (atom)
        MMDB->SelectResidue ( selHndResB,atom->GetResidue(),
                              mmdb::STYPE_RESIDUE,mmdb::SKEY_OR,false );
    }
    MMDB->MakeSelIndex ( selHndResB );

    sprintf ( S,
    " %s\n"
    " -----.-.------------.--.----------------------\n"
    "   ## |I|  Struct. %i |HS|   ASA    BSA  DeltaG\n"
    " -----+-+------------+--+----------------------\n",
      Title,structNo );
    f.Write ( S );

    for (i=0;i<nRes;i++)  {

      if (Res[i]->isInSelection(selHndResI))
           strcpy ( rid,"I" );
      else if (Res[i]->isInSelection(selHndResS))
           strcpy ( rid,"s" );
      else strcpy ( rid," " );

      strcpy ( sb,"  " );
      if (Res[i]->isInSelection(selHndResH))  sb[0] = 'H';
      if (Res[i]->isInSelection(selHndResB))  sb[1] = 'S';

      strcpy ( chID,Res[i]->GetChainID() );
      if (!chID[0])  strcpy ( chID,"-" );

      sprintf ( S,"%5i |%1s| %1s:%3s%4i%1s |%2s| %6.2f %6.2f %6.2f" ,
                  i+1,rid,chID,Res[i]->GetResName(),
                  Res[i]->GetSeqNum(),Res[i]->GetInsCode(),sb,
                  D->asa[i],bsa[i],solvEn[i] );

      f.WriteLine ( S );

      for (j=0;j<nStab;j++)
        if (sre_stab[j].resNo==i)  {
          sprintf ( rid," %6.2f\n",sre_stab[j].effect );
          mmdb::CreateCopCat ( stab_s[j]," ",S,rid );
        }

      for (j=0;j<nDestab;j++)
        if (sre_destab[j].resNo==i)  {
          sprintf ( rid," %6.2f\n",sre_destab[j].effect );
          mmdb::CreateCopCat ( destab_s[j]," ",S,rid );
        }

    }

    f.WriteLine (
      " -----'-'------------'--'----------------------\n" );

    if (stab_s[0])  {
      sprintf ( S,
        " Residues with most stabilizing effect\n"
        " ------.-.------------.--.-----------------------------\n"
        "    ## |I|  Struct. %i |HS|   ASA    BSA  DeltaG  Total\n"
        " ------+-+------------+--+-----------------------------\n",
        structNo );
      f.Write ( S );
      for (j=0;j<nStab;j++)
        if (stab_s[j])  {
          f.Write ( stab_s[j] );
          delete[] stab_s[j];
        }
      f.WriteLine (
        " ------'-'------------'--'-----------------------------\n" );
    }
    delete[] stab_s;

    if (destab_s[0])  {
      sprintf ( S,
        " Residues with most destabilizing effect\n"
        " ------.-.------------.--.-----------------------------\n"
        "    ## |I|  Struct. %i |HS|   ASA    BSA  DeltaG  Total\n"
        " ------+-+------------+--+-----------------------------\n",
        structNo );
      f.Write ( S );
      for (j=0;j<nDestab;j++)
        if (destab_s[j])  {
          f.Write ( destab_s[j] );
          delete[] destab_s[j];
        }
      f.WriteLine (
        " ------'-'------------'--'-----------------------------\n" );
    }
    delete[] destab_s;

    MMDB->DeleteSelection ( selHndRes  );
    MMDB->DeleteSelection ( selHndResS );
    MMDB->DeleteSelection ( selHndResI );
    MMDB->DeleteSelection ( selHndResH );
    MMDB->DeleteSelection ( selHndResB );

  }


  RESULT_CODE Detail::DetailInterface ( mmdb::cpstr sessionName,
                                        mmdb::cpstr fileName,
                                        int serialNo )  {
  mmdb::io::File f;
  PInterface     interface;
  char           S[100];
  RESULT_CODE    rc;

    //   1. Check configuration

    if (ConfStatus()!=CFG_Configured)
      return RESULT_ConfigurationError;

    //   2. Check session directory and results

    switch (checkCrDir(sessionName))  {

      case SDIR_doesntExist : return RESULT_SessionDoesntExist;
      case SDIR_noResults   : return RESULT_noSessionResults;

      default : ;

    }

    //   3. Open output file

    rc = StartTextOutput ( f,fileName );
    if (rc!=RESULT_Ok)  return rc;

    //   4. Read and check data

    rc = readPIData();
    if (rc!=RESULT_Ok)  return rc;

    if (query->getNofInterfaces()<=0)  {
      f.WriteLine ( "\n NO INTERFACES FOUND\n" );
      return RESULT_interfaceNoOutOfRange;
    }

    interface = query->getInterface ( serialNo-1 );
    if (!interface)
      return RESULT_interfaceNoOutOfRange;

    rc = readStructure();
    if (rc!=RESULT_Ok)  return rc;

    //   5. Print output tables

    sprintf ( S," INTERFACE No. %i\n",interface->id );
    f.WriteLine ( S );

    printInterfaceSummary ( f,interface );
    printBondTable        ( f," 2. Hydrogen Bonds",
                              "  2. No Hydrogen Bonds found\n",
                              interface,interface->HBond,
                              interface->nHBonds );
    printBondTable        ( f," 3. Salt Bridges",
                              "  3. No Salt Bridges found\n",
                              interface,interface->SBridge,
                              interface->nSBridges );
    printIntResTable      ( f," 4. Interfacing Residues: Structure 1",
                              "  4. No Interfacing residues in "
                              "structure 1 found",
                              interface,1 );
    printIntResTable      ( f," 5. Interfacing Residues: Structure 2",
                              "  5. No Interfacing Residues in "
                              "Structure 2 found",
                              interface,2 );

    f.shut();

    return RESULT_Ok;

  }


  void Detail::printMonomerSummary ( mmdb::io::RFile f, int domainNo )  {
  PDomain D;
  char     S[1000],DeltaG[20];
  mmdb::realtype cryst_bsa,cryst_DeltaG;
  mmdb::realtype totDeltaG;
  int      cryst_nHB,cryst_nSB,cryst_nDS,cryst_nAt,cryst_nRes;
  int      i;

    D = query->getDomain ( domainNo );

    query->calcDomainStats ( domainNo,
                             cryst_bsa,cryst_DeltaG,cryst_nHB,
                             cryst_nSB,cryst_nDS,cryst_nAt,cryst_nRes );

    totDeltaG = 0.0;
    for (i=0;i<D->nRes;i++)
      totDeltaG -= D->solvEn[i];

    if (D->dclass==DCLASS_Protein)
          sprintf ( DeltaG,"%8.1f(*)",D->DeltaG );
    else  strcpy  ( DeltaG,"N/A   " );

    sprintf ( S,
  " 1. Monomer Summary\n"
  " ------------------.---------------------------------\n"
  "                   | Interfaces  Surface    Total\n"
  " ------------------+---------------------------------\n"
  " Atoms             |  %8i %8i  %8i\n"
  " Residues          |  %8i %8i  %8i\n"
  " Area, sq.A        |  %8.1f %8.1f\n"
  " Delta G, kcal/mol |  %8.1f %8.1f %12s\n"
  " Hydrogen bonds    |  %8i\n"
  " Salt bridges      |  %8i\n"
  " Disulphide bonds  |  %8i\n"
  " ------------------'---------------------------------\n",
    cryst_nAt   ,D->nSurfAtoms,D->nAtoms,
    cryst_nRes  ,D->nSurfRes  ,D->nRes,
    cryst_bsa   ,D->surfArea  ,
    cryst_DeltaG,totDeltaG    ,DeltaG,
    cryst_nHB,
    cryst_nSB,
    cryst_nDS
           );

    f.Write ( S );
    if (D->dclass==DCLASS_Protein)
         f.WriteLine ( " (*) solvation energy gain of folding\n" );
    else f.LF();

  }

  void Detail::printResidueData ( mmdb::io::RFile f, int domainNo )  {
  PDomain         D;
  mmdb::PManager  MMDB;
  mmdb::PPResidue Res;
  PPInterface     interface;
  char            S[500];
  mmdb::ChainID   chID;
  mmdb::mat44     T;
  mmdb::rvector   symf;
  mmdb::realtype  ssum0,bsa,bsa0;
  mmdb::realtype  DeltaG0,DeltaGi,DeltaGi0;
  int             selHndRes;
  int             nRes,nInterfaces;
  int             i,j;

    D    = query->getDomain ( domainNo );
    MMDB = query->MMDB;
    selHndRes  = MMDB->NewSelection();

    D->SelectDomain   ( selHndRes,MMDB,mmdb::STYPE_RESIDUE,mmdb::SKEY_NEW,1,false );
    MMDB->GetSelIndex ( selHndRes,Res,nRes );

    nInterfaces = query->getNofInterfaces();
    interface   = query->getInterfaces();
    mmdb::GetVectorMemory ( symf,nInterfaces,0 );
    for (i=0;i<nInterfaces;i++)  {
      symf[i] = -1.0;
      if ((interface[i]->domain1==domainNo) ||
          (interface[i]->domain2==domainNo))  {
        symf[i] = 1.0;
        if (interface[i]->domain1==interface[i]->domain2) {
          // this may be a symmetry Pi-interface, in which case
          // we should half its interface area on output
          mmdb::Mat4Mult ( T,interface[i]->TMatrix,
                             interface[i]->TMatrix );
          if (mmdb::isMat4Unit(T,0.00001,true))
            symf[i] = 0.5;
        }
      }
    }


    f.Write (
  " 2. Residue Accessibility and Solvation Energy Effect\n"
  " -----.------------.---------------------------------------------\n"
  "   ## |   Residue  |       ASA       BSA      Delta G   Delta Gi\n"
  " -----+------------+---------------------------------------------\n"
            );


    ssum0    = 0.0;
    bsa0     = 0.0;
    DeltaG0  = 0.0;
    DeltaGi0 = 0.0;

    for (i=0;i<nRes;i++)  {

      bsa     = 0.0;
      DeltaGi = 0.0;
      for (j=0;j<nInterfaces;j++)
        if (symf[j]>0.0)  {
          if (interface[j]->domain1==domainNo)  {
            bsa     += interface[j]->bsa1[i];
            DeltaGi -= interface[j]->solvEn1[i];
          }
          if (interface[j]->domain2==domainNo)  {
            bsa     += interface[j]->bsa2[i];
            DeltaGi -= interface[j]->solvEn2[i];
          }
          bsa      *= symf[j];
          DeltaGi  *= symf[j];
          bsa0     += bsa;
          DeltaGi0 += DeltaGi;
        }
      ssum0   += D->asa[i];
      DeltaG0 += D->solvEn[i];

      strcpy ( chID,Res[i]->GetChainID() );
      if (!chID[0])  strcpy ( chID,"-" );

      sprintf ( S,"%5i | %1s:%3s%4i%1s | %9.2f %9.2f   %9.2f %9.2f",
                  i+1,chID,Res[i]->GetResName(),
                  Res[i]->GetSeqNum(),Res[i]->GetInsCode(),
                  D->asa[i],bsa,D->solvEn[i],DeltaGi );
      f.WriteLine ( S );

    }

    sprintf ( S,
  " -----'------------+---------------------------------------------\n"
  "             Total | %9.2f %9.2f   %9.2f %9.2f\n"
  " ------------------'---------------------------------------------\n"
  "  ASA:      accessible surface area, sq.A\n"
  "  BSA:      buried (in all interfaces) surface area, sq.A\n"
  "  Delta G:  solvation energy contribution, kcal/mol\n"
  "  Delta Gi: solvation energy effect of all interfaces, kcal/mol\n",
      ssum0,bsa0,DeltaG0,DeltaGi0 );
    f.WriteLine ( S );

    mmdb::FreeVectorMemory ( symf,0 );

  }


  RESULT_CODE Detail::DetailMonomer ( mmdb::cpstr sessionName,
                                      mmdb::cpstr fileName,
                                      int serialNo )  {
  mmdb::io::File f;
  char           S[200],dName[100];
  RESULT_CODE    rc;

    //   1. Check configuration

    if (ConfStatus()!=CFG_Configured)
      return RESULT_ConfigurationError;

    //   2. Check session directory and results

    switch (checkCrDir(sessionName))  {

      case SDIR_doesntExist : return RESULT_SessionDoesntExist;
      case SDIR_noResults   : return RESULT_noSessionResults;

      default : ;

    }

    //   3. Open output file

    rc = StartTextOutput ( f,fileName );
    if (rc!=RESULT_Ok)  return rc;

    //   4. Read and check data

    rc = readPIData();
    if (rc!=RESULT_Ok)  return rc;

    if (query->getNofDomains()<=0)  {
      f.WriteLine ( "\n NO MONOMERS FOUND\n" );
      return RESULT_monomerNoOutOfRange;
    }

    if (!query->getDomain(serialNo-1))
      return RESULT_monomerNoOutOfRange;

    rc = readStructure();
    if (rc!=RESULT_Ok)  return rc;

    sprintf ( S," MONOMER No. %i: %s\n",serialNo,
                query->getDomainRange(dName,serialNo-1,0) );
    f.WriteLine ( S );

    printMonomerSummary ( f,serialNo-1 );
    printResidueData    ( f,serialNo-1 );

    f.shut();

    return RESULT_Ok;

  }


  #define _leadStr  mmdb::pstr("                       ")

  void Detail::printAssemblySummary ( mmdb::io::RFile f, PAssembly A )  {
  char       S[2000];
  mmdb::pstr F,C,P;

    F = NULL;
    C = NULL;
    P = NULL;

    sprintf ( S,
    " 1. Assembly Summary\n"
    " ---------------------------------------------------------------\n"
    " Multimeric state         %4i\n"
    " Copies in unit cell      %4i\n"
    " Symmetry number          %4i\n"
    " Surface area             %9.1f sq.A\n"
    " Buried area              %9.1f sq.A\n"
    " Delta G diss             %9.1f kcal/mol\n"
    " Entropy of diss-n        %9.1f kcal/mol\n"
    " Ground Delta G diss      %9.1f kcal/mol\n"
    " Ground Entropy of diss-n %9.1f kcal/mol\n"
    " Formula:                 %s\n"
    " Composition:             %s\n"
    " Dissociation pattern:    %s\n"
    " ---------------------------------------------------------------\n",
      A->mmSize, A->nUC, A->symNumber, A->asa, A->bsa,
      A->freeEn,A->entropy,A->freeEn0,A->entropy0,
      A->getFormula(F,query->domains,60,_leadStr),
      A->getComposition(C,query->domains,60,_leadStr),
      A->getDissPattern(P,query->domains,60,_leadStr)
       );
    f.WriteLine ( S );

    if (F)  delete[] F;
    if (C)  delete[] C;
    if (P)  delete[] P;

  }

  void Detail::printEngagedInterfaces ( mmdb::io::RFile f, PAssembly A )  {
  PInterface     interface;
  mmdb::ivector  intf,dintf;
  char           S[1000],name1[100],name2[100],diss[5];
  int            nInterfaces,i,symopNo;

    nInterfaces = query->getNofInterfaces();

    mmdb::GetVectorMemory   ( intf,nInterfaces,0 );
    A->getEngagedInterfaces ( intf,nInterfaces );

    mmdb::GetVectorMemory ( dintf,nInterfaces,0 );
    A->getDissInterfaces  ( dintf,nInterfaces   );

    f.Write (
   " 2. List of Engaged Interfaces\n"
   " ------.-------------.-------------.---.-"
          "-------.------------------------\n"
   " ## Id |   Monomer1  |   Monomer2  | D | "
          "Sym.Id |   Area  DeltaG Nhb Nsb\n"
   " ------+-------------+-------------+---+-"
          "-------+------------------------\n"
      );


    for (i=0;i<nInterfaces;i++)
      if (intf[i]>0)  {
        interface = query->getInterface ( i );
        if (interface)  {

          if (query->PI->rcsb_symops)  symopNo = interface->rcsb_symop;
                                 else  symopNo = interface->symOpNo+1;

          if (dintf[i]>0)  strcpy ( diss,"*" );
                     else  strcpy ( diss," " );

          sprintf ( S,
            " %2i %2i | %11s | %11s | %1s | %2i_%1i%1i%1i | "
            "%7.1f %6.1f %3i %3i ",
            i+1,interface->type,
            query->getDomainRange(name1,interface->domain1,11),
            query->getDomainRange(name2,interface->domain2,11),
            diss,symopNo,
            interface->cell_i+nCellOut,
            interface->cell_j+nCellOut,
            interface->cell_k+nCellOut,
            interface->intArea,
            interface->intDeltaG,
            interface->nHBonds,
            interface->nSBridges
                  );

          f.WriteLine ( S );

        }
      }

    f.WriteLine (
   " ------'-------------'-------------'---'-"
          "-------'------------------------\n"
   " ##:   serial number             Monomer: selection range\n"
   " Id:   interface type id         D:       dissociating interface\n"
   " Area: sq. angstrom              DeltaG:  kcal/mol\n"
   " Nhb:  no. of hydrogen bonds     Nsb:     no. of salt bridges\n"
                );

    mmdb::FreeVectorMemory ( intf ,0 );
    mmdb::FreeVectorMemory ( dintf,0 );

  }

  void Detail::printAssembledMonomers ( mmdb::io::RFile f, PAssembly A )  {
  PDomain        domain;
  PPDomain       domains;
  PAsmUnit       M;
  PPInterface    interface;
  char           S[1000],name[100];
  mmdb::realtype bsa,DeltaG;
  int            i,nOcc,nInterfaces,nHB,nSB,nDS,nAt,nRes,symopNo;

    f.Write (
   " 3. List of Assembled Monomers\n"
   " -------.-------------.-------------------"
          "----------------------------------------\n"
   "  ## Id |   Monomer   |     ASA      BSA  DeltaG Nhb Nsb Sym.Id "
  "Symmetry op-n\n"
   " -------+-------------+-------------------"
          "----------------------------------------\n"
      );

    nInterfaces = query->getNofInterfaces();
    interface   = query->getInterfaces();
    domains     = query->getDomains();

    for (i=0;i<A->asmSize;i++)  {
      M      = A->M[i];
      domain = query->getDomain ( M->id );
      if (domain) {

        A->calcStats ( M->id,nOcc,bsa,DeltaG,nHB,nSB,nDS,nAt,nRes,
                       interface,domains,nInterfaces );

        if (query->A->rcsb_symops)  symopNo = M->rcsb_symop;
                              else  symopNo = M->symOpNo+1;

        sprintf ( S,
   "%4i %2i | %11s | %8.1f %8.1f %6.1f %3i %3i %2i_%1i%1i%1i %s",
          i+1,M->id+1,query->getDomainRange(name,M->id,11),
          domain->surfArea,bsa,DeltaG,nHB,nSB,
          symopNo,M->cell_i+nCellOut,M->cell_j+nCellOut,
          M->cell_k+nCellOut,M->symOp );

        f.WriteLine ( S );
      }
    }

    f.WriteLine (
   " -------'-------------'-------------------"
      "----------------------------------------\n"
   "  Id:  monomer serial no.    "
   "ASA/BSA: total/buried surface area, sq.A\n"
   "  DeltaG: solvation energy effect upon complexation, kcal/mol\n"  );

  }

  void Detail::printTransformations ( mmdb::io::RFile f, PAssembly A ) {
  PDomain        domain;
  PAsmUnit      M;
  mmdb::ivector  symopNo,cell_i,cell_j,cell_k;
  char           S[1000];
  int            i,j,k;
  bool           B;

    f.Write (
   " 4. Transformation matrices\n"
   " --------.-------------------------------------------\n"
   "  Sym.Id |        Rotation            Translation\n"
      );

    mmdb::GetVectorMemory ( symopNo,A->asmSize,0 );
    mmdb::GetVectorMemory ( cell_i ,A->asmSize,0 );
    mmdb::GetVectorMemory ( cell_j ,A->asmSize,0 );
    mmdb::GetVectorMemory ( cell_k ,A->asmSize,0 );

    k = 0;
    for (i=0;i<A->asmSize;i++)  {
      M      = A->M[i];
      domain = query->getDomain ( M->id );
      if (domain) {

        if (query->A->rcsb_symops)  symopNo[k] = M->rcsb_symop;
                              else  symopNo[k] = M->symOpNo+1;
        cell_i[k] = M->cell_i+nCellOut;
        cell_j[k] = M->cell_j+nCellOut;
        cell_k[k] = M->cell_k+nCellOut;
        B = false;
        for (j=0;(j<k) && (!B);j++)
          B = (symopNo[j]==symopNo[k]) &&
              (cell_i[j]==cell_i[k]) &&
              (cell_j[j]==cell_j[k]) &&
              (cell_k[j]==cell_k[k]);

        if (!B)  {
          f.WriteLine (
   " --------+-------------------------------------------" );

          sprintf ( S,
   "         |  %6.3f %6.3f %6.3f  %12.3f\n"
   " %3i_%1i%1i%1i |  %6.3f %6.3f %6.3f  %12.3f\n"
   "         |  %6.3f %6.3f %6.3f  %12.3f",
           M->T[0][0],M->T[0][1],M->T[0][2], M->T[0][3],
           symopNo[k], cell_i[k], cell_j[k], cell_k[k],
           M->T[1][0],M->T[1][1],M->T[1][2], M->T[1][3],
           M->T[2][0],M->T[2][1],M->T[2][2], M->T[2][3] );
          f.WriteLine ( S );
          k++;
        }
      }

    }


    f.WriteLine (
   " --------'-------------------------------------------\n" );

  }


  RESULT_CODE Detail::DetailAssembly ( mmdb::cpstr sessionName,
                                       mmdb::cpstr fileName,
                                       int  serialNo )  {
  mmdb::io::File f;
  PAssembly      A;
  char           S[200];
  int            i,j;
  RESULT_CODE    rc;


    //   1. Check configuration

    if (ConfStatus()!=CFG_Configured)
      return RESULT_ConfigurationError;


    //   2. Check session directory and results

    switch (checkCrDir(sessionName))  {

      case SDIR_doesntExist : return RESULT_SessionDoesntExist;
      case SDIR_noResults   : return RESULT_noSessionResults;

      default : ;

    }


    //   3. Read and check data

    rc = readPIData();
    if (rc!=RESULT_Ok)  return rc;

    rc = readAssemblies();
    if (rc!=RESULT_Ok)  return rc;

    if ((!query->A) || (query->asmStatus!=ASSMB_Ok))  {
      f.assign ( "stdout",true );
      f.rewrite();
      makeNoAssembliesPage ( f );
      f.shut();
      return RESULT_assemblyNoOutOfRange;
    }

    A = NULL;
    for (i=0;(i<query->A->nCSRes) && (!A);i++)
      for (j=0;(j<query->A->crystSplit[i]->nAssemblies) && (!A);j++)
        if (query->A->crystSplit[i]->A[j]->serNo==serialNo)
           A = query->A->crystSplit[i]->A[j];

    if (!A)  return RESULT_assemblyNoOutOfRange;


    //   4. Open output file

    rc = StartTextOutput ( f,fileName );
    if (rc!=RESULT_Ok)  return rc;

    sprintf ( S," Assembly No. %i, type %i\n",serialNo,A->type+1 );
    f.WriteLine ( S );

    printAssemblySummary   ( f,A );
    printEngagedInterfaces ( f,A );
    printAssembledMonomers ( f,A );
    printTransformations   ( f,A );

    f.shut();

    return RESULT_Ok;

  }


  DefineStructure(SRem350);

  struct SRem350  {
    mmdb::pstr  c;
    mmdb::mat44 T;
  };

  void  Detail::printRemark350 ( mmdb::io::RFile  f, PAssembly  A,
                                 bool fractional )  {
  PSRem350       r;
  mmdb::pstr     p;
  mmdb::ChainID  chID;
  char           L[1000],LL[100];
  mmdb::mat44    T;
  int            i,j,k,m, ii,jj;
  bool           firstLine;
  char           c;

    r = new SRem350[A->asmSize];
    for (i=0;i<A->asmSize;i++)
      r[i].c = NULL;

    k = 0;
    for (i=0;i<A->asmSize;i++)
      if (query->getDomainClass(A->M[i]->id)!=DCLASS_Ligand)  {
        if (fractional)
          query->MMDB->Orth2Frac ( A->M[i]->T,T );
        else
          mmdb::Mat4Copy ( A->M[i]->T,T );
        m = -1;
        for (j=0;(j<k) && (m<0);j++)  {
          m = j;
          for (ii=0;(ii<3) && (m==j);ii++)
            for (jj=0;(jj<4) && (m==j);jj++)
              if (fabs(T[ii][jj]-r[j].T[ii][jj])>0.00001)
                m = -1;
        }
        if (m<0)  {
          mmdb::Mat4Copy   ( T,r[k].T );
          mmdb::CreateCopy ( r[k].c,query->getChainID(chID,A->M[i]->id) );
          k++;
        } else
          mmdb::CreateConcat ( r[m].c,
                               query->getChainID(chID,A->M[i]->id) );
//          mmdb::CreateConcat ( r[m].c,", ",
//                               query->getChainID(chID,A->M[i]->id) );
      }

    sprintf ( L,
"REMARK 350\n" 
"REMARK 350 BIOMOLECULE: %i\n" 
"REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: %s\n"
"REMARK 350 SOFTWARE USED: PISA\n"
"REMARK 350 TOTAL BURIED SURFACE AREA: %.0f ANGSTROM**2\n"
"REMARK 350 SURFACE AREA OF THE COMPLEX: %.0f ANGSTROM**2\n"
"REMARK 350 CHANGE IN SOLVENT FREE ENERGY: %.1f KCAL/MOL\n",
      A->serNo,A->sizeName(LL),A->bsa,A->asa,A->seGain );
    f.Write ( L );

    for (i=0;i<k;i++)  {
      p = r[i].c;
      ii = strlen(p);
      for (j=0;j<ii;j++)
        for (m=j+1;m<ii;m++)
          if (p[j]>p[m])  {
            c    = p[j];
            p[j] = p[m];
            p[m] = c;
          }
      r[i].c = new char[4*ii];
      r[i].c[0] = p[0];
      m = 1;
      for (j=1;j<ii;j++)  {
        r[i].c[m++] = ',';
        r[i].c[m++] = ' ';
        r[i].c[m++] = p[j];
      }
      r[i].c[m] = char(0);
      delete[] p;
      p = r[i].c;
      firstLine = true;
      while (*p)  {
        j = strlen(p);
        if (j>26)
          p[26] = char(0);
        if (firstLine)  {
          f.Write ( "REMARK 350 APPLY THE FOLLOWING TO CHAINS: " );
          firstLine = false;        
        } else
          f.Write ( "REMARK 350                    AND CHAINS: " );
        f.WriteLine ( p );
        if (j>26)  p = &(p[27]);
             else  p = &(p[j]);
      }
      sprintf ( L,
        "REMARK 350   BIOMT1 %3i %9.6f %9.6f %9.6f     %10.5f\n"
        "REMARK 350   BIOMT2 %3i %9.6f %9.6f %9.6f     %10.5f\n"
        "REMARK 350   BIOMT3 %3i %9.6f %9.6f %9.6f     %10.5f\n",
        i+1,r[i].T[0][0],r[i].T[0][1],r[i].T[0][2],r[i].T[0][3],
        i+1,r[i].T[1][0],r[i].T[1][1],r[i].T[1][2],r[i].T[1][3],
        i+1,r[i].T[2][0],r[i].T[2][1],r[i].T[2][2],r[i].T[2][3] );
      f.Write ( L );
    }

    for (i=0;i<A->asmSize;i++)
      if (r[i].c)  delete[] r[i].c;
    delete[] r;

  }


  void  Detail::printRemark350 ( mmdb::io::RFile  f, PDomain  D,
                                 int domainNo )  {
  char  L[1000],LL[100];

    sprintf ( L,
"REMARK 350\n" 
"REMARK 350 BIOMOLECULE: %i\n" 
"REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: MONOMERIC\n"
"REMARK 350 SOFTWARE USED: PISA\n"
"REMARK 350 SURFACE AREA OF THE MONOMER: %.0f ANGSTROM**2\n"
"REMARK 350 FREE ENERGY OF FOLDING: %.1f KCAL/MOL\n"
"REMARK 350 APPLY THE FOLLOWING TO CHAINS: %s\n"
"REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000\n"
"REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000\n"
"REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000\n",
      domainNo,D->surfArea,D->DeltaG,D->getChainID(LL) );
    f.Write ( L );

  }

  mmdb::cpstr rem300_text =
"\n"
"REMARK 300 SEE REMARK 350 FOR THE AUTHOR PROVIDED AND/OR PROGRAM\n"
"REMARK 300 GENERATED ASSEMBLY INFORMATION FOR THE STRUCTURE IN\n"
"REMARK 300 THIS ENTRY. THE REMARK MAY ALSO  PROVIDE INFORMATION ON\n"
"REMARK 300 BURIED SURFACE AREA.\n";

  mmdb::cpstr rem350_header =
"REMARK 350\n"
"REMARK 350 COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN\n"
"REMARK 350 BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE\n"
"REMARK 350 MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS\n"
"REMARK 350 GIVEN BELOW. BOTH NON-CRYSTALLOGRAPHIC AND\n"
"REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.\n";


  void Detail::printRemark350M ( mmdb::io::RFile f )  {
  // all assemblies are monomeric
  PDomain        D;
  char           S[200];
  int            i,j;

    f.Write ( "REMARK 300\n"
              "REMARK 300 BIOMOLECULE:" );

    i = 0;
    for (j=0;j<query->domains->nDomains;j++)  {
      D = query->domains->domain[j];
      if (D)  {
        if (D->dclass!=DCLASS_Ligand)  {
          i++;
          if (i==1)  sprintf ( S," %i" ,i );
               else  sprintf ( S,", %i",i );
          f.Write ( S );
        }
      }
    }

    f.Write ( rem300_text   );
    f.Write ( rem350_header );

    i = 0;
    for (j=0;j<query->domains->nDomains;j++)  {
      D = query->domains->domain[j];
      if (D)  {
        if (D->dclass!=DCLASS_Ligand)  {
          i++;
          printRemark350 ( f,D,i );
        }
      }
    }

  }

  RESULT_CODE Detail::Remark350 ( mmdb::cpstr sessionName,
                                  mmdb::cpstr fileName,
                                  int  serialNo )  {
  mmdb::io::File f;
  PAssembly      A;
  char           S[200];
  int            i,j,sNo;
  RESULT_CODE    rc;

    if (sessionName)  {

      //   1. Check configuration

      if (ConfStatus()!=CFG_Configured)
        return RESULT_ConfigurationError;


      //   2. Check session directory and results

      switch (checkCrDir(sessionName))  {

        case SDIR_doesntExist : return RESULT_SessionDoesntExist;
        case SDIR_noResults   : return RESULT_noSessionResults;

        default : ;

      }


      //   3. Read and check data

      rc = readPIData();
      if (rc!=RESULT_Ok)  return rc;

      rc = readStructure();
      if (rc!=RESULT_Ok)  return rc;

      rc = readAssemblies();
      if (rc!=RESULT_Ok)  return rc;

    }
    
    if (((!query->A) && (serialNo<mmdb::MaxInt4)) ||
         (query->asmStatus!=ASSMB_Ok))  {
      if (sessionName)
        f.assign ( "stdout",true );
      else  {
        rc = StartTextOutput ( f,fileName );
        if (rc!=RESULT_Ok)  return rc;
      }
      f.LF();
      makeNoAssembliesPage ( f );
      f.shut();
      return RESULT_assemblyNoOutOfRange;
    }
    
    if (serialNo==mmdb::MaxInt4)  {
      sNo = -query->A->nCSRes-1;
      if (query->A->nCSRes>0)  {
        if (query->A->crystSplit[0]->Score<3)
          sNo = -1;
      }
    } else
      sNo = serialNo;

    if (sNo<0)  {
    
      // make complete remark 350 for the whole 1st split
      
      i = -sNo-1;

      if (i>=query->A->nCSRes)  {
        if (serialNo==mmdb::MaxInt4)  {
           rc = StartTextOutput ( f,fileName );
           if (rc!=RESULT_Ok)  return rc;
           printRemark350M ( f );
           f.shut();
        } else
          return RESULT_assemblyNoOutOfRange;
      } else  {

        rc = StartTextOutput ( f,fileName );
        if (rc!=RESULT_Ok)  return rc;
  
        f.Write ( "REMARK 300\n"
                  "REMARK 300 BIOMOLECULE:" );

        for (j=0;j<query->A->crystSplit[i]->nAssemblies;j++)  {
          A = query->A->crystSplit[i]->A[j];
          if (A)  {
            if (j==0)  sprintf ( S," %i" ,A->serNo );
                 else  sprintf ( S,", %i",A->serNo );
            f.Write ( S );
          }
        }

        f.Write ( rem300_text   );
        f.Write ( rem350_header );

        for (j=0;j<query->A->crystSplit[i]->nAssemblies;j++)  {
          A = query->A->crystSplit[i]->A[j];
          if (A)
            printRemark350 ( f,A,false );
        }

        f.shut();

      }

    } else  {
    
      if (serialNo>=1)  {
  
        A = NULL;
        for (i=0;(i<query->A->nCSRes) && (!A);i++)
          for (j=0;(j<query->A->crystSplit[i]->nAssemblies) && (!A);j++)
            if (query->A->crystSplit[i]->A[j]->serNo==serialNo)
               A = query->A->crystSplit[i]->A[j];
  
        if (!A)  return RESULT_assemblyNoOutOfRange;
  
        sprintf ( S," Assembly No. %i, type %i\n",serialNo,A->type+1 );
  
      } else if (query->Complex)  {
  
        A = query->Complex;
        sprintf ( S," Complex represented by the content of asu\n" );
  
      } else  {
  
        if (sessionName)  f.assign ( "stdout",true );
        else  {
          rc = StartTextOutput ( f,fileName );
          if (rc!=RESULT_Ok)  return rc;
        }
  
        f.WriteLine ( "\n"
        " Analysis of complex, represented by ASU, was not performed\n"
        " (possibly because it is monomeric).\n"
        "\n" );
        f.shut();
  
        return RESULT_assemblyNoOutOfRange;
  
      }
  
      //   4. Open output file
  
      rc = StartTextOutput ( f,fileName );
      if (rc!=RESULT_Ok)  return rc;
  
      f.WriteLine ( S );
  
      f.WriteLine ( " 1. Remark 350 in orthogonal coordinates\n" );
      f.Write     ( rem350_header );
      printRemark350 ( f,A,false );
   
      f.WriteLine ( " 2. Remark 350 in fractional coordinates\n" );
      f.Write     ( rem350_header );
      printRemark350 ( f,A,true  );
  
      f.shut();
      
    }

    return RESULT_Ok;

  }


  RESULT_CODE Detail::MakeInterfacesXML ( mmdb::cpstr sessionName,
                                          mmdb::cpstr fileName )  {
  mmdb::xml::PXMLObject xml;
  RESULT_CODE           rc;
  PInterface Interface;

    if (sessionName)  {
      
      //   1. Check configuration

      if (ConfStatus()!=CFG_Configured)
        return RESULT_ConfigurationError;

      //   2. Check session directory and results

      switch (checkCrDir(sessionName))  {

        case SDIR_doesntExist : return RESULT_SessionDoesntExist;
        case SDIR_noResults   : return RESULT_noSessionResults;

        default : ;

      }

      //   3. Read and check data

      rc = readPIData();
      if (rc!=RESULT_Ok)  return rc;

      rc = readStructure();
      if (rc!=RESULT_Ok)  return rc;

    }
    
    xml = query->getInterfacesXML(asisKey);
    if (!xml)  {
      makeNoInterfacesXML ( sessionName,fileName );
      return RESULT_Ok;
    } else  {
      if (xml->WriteObject(fileName,0,2)!=mmdb::xml::XMLRC_Ok)
            rc = RESULT_cantWriteFile;
      else  rc = RESULT_Ok;
      delete xml;
    }

    return rc;

  }

  RESULT_CODE Detail::MakeAssembliesXML ( mmdb::cpstr sessionName,
                                          mmdb::cpstr fileName )  {
  mmdb::xml::PXMLObject xml,complex_xml;
  RESULT_CODE           rc;

    if (sessionName)  {

      
      //   1. Check configuration

      if (ConfStatus()!=CFG_Configured)
        return RESULT_ConfigurationError;

      //   2. Check session directory and results

      switch (checkCrDir(sessionName))  {

        case SDIR_doesntExist : return RESULT_SessionDoesntExist;
        case SDIR_noResults   : return RESULT_noSessionResults;

        default : ;

      }

      //   3. Read and check data

      rc = readPIData();
      if (rc!=RESULT_Ok)  return rc;

      rc = readAssemblies();
      if (rc!=RESULT_Ok)  return rc;

    }


    if ((!query->A) || (query->asmStatus!=ASSMB_Ok))
          xml = makeNoAssembliesXML ( sessionName );
    else  xml = query->A->getAssembliesXML ( sessionName,
                                             query->domains,
                                             query->PI,
                                             nCellOut );


    if (query->Complex)  {

      complex_xml = new mmdb::xml::XMLObject ( xml_asu_complex );
      
      if(asisKey==AS_IS_off){
      complex_xml->AddObject (
      query->Complex->getAssemblyXML ( query->domains,query->PI,
				       nCellOut,-1) );
      }
      if (asisKey==AS_IS_on){
	complex_xml->AddObject (
	     query->Complex->getasisAssemblyXML ( query->domains,query->PI,
		      nCellOut) );	
		      }
      xml->AddObject ( complex_xml );
    }

    if (xml->WriteObject(fileName,0,2)!=mmdb::xml::XMLRC_Ok)
          rc = RESULT_cantWriteFile;
    else  rc = RESULT_Ok;

    delete xml;

    return rc;

  }

}  // namespace pisa
