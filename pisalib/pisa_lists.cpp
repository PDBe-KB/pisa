// $Id: pisa_lists.cpp $
// =================================================================
//
//    15.03.19   <--  Date of Last Modification.
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
//  (C) E. Krissinel 2007-2019
//
// =================================================================
//

#include <string.h>

#include "pisa_lists.h"
#include "pisa_types.h"

namespace pisa  {

  // ========================  Lists  ===========================

  Lists::Lists ( mmdb::cpstr confPath ) : Data ( confPath )  {
    InitLists();
  }

  Lists::~Lists()  {}

  void Lists::InitLists()  {}


  RESULT_CODE  Lists::ListInterfaces ( mmdb::cpstr sessionName,
                                       mmdb::cpstr fileName )  {
  mmdb::io::File f;
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

    rc = readPIData();
    if (rc!=RESULT_Ok)  return rc;

    rc = StartTextOutput ( f,fileName );
    if (rc!=RESULT_Ok)  return rc;

    rc = ListInterfaces ( f );
    f.shut();

    return rc;

  }


  RESULT_CODE Lists::ListInterfaces (  mmdb::io::RFile f ) {
  PInterface interface;
  char       S[1000],name1[20],name2[20],fix[3];
  int        nInterfaces,i, symopNo;

    nInterfaces = query->getNofInterfaces();

    if (nInterfaces<=0)  {
      f.WriteLine ( "\n NO INTERFACES FOUND\n" );
    } else  {

        f.Write (
       " LIST OF INTERFACES\n"
       " -------.-------------.----------------------------------"
          "-------.---------------------------\n"
       "  ## Id |   Monomer1  |   Monomer2    Symmetry operation "
          "Sym.Id |   Area  DeltaG Nhb Nsb Nds\n"
       " -------+-------------+----------------------------------"
          "-------+---------------------------\n"
      );

      for (i=0;i<nInterfaces;i++)  {
        interface = query->getInterface ( i );
        if (interface)  {

          if (query->PI->rcsb_symops)  symopNo = interface->rcsb_symop;
                                 else  symopNo = interface->symOpNo+1;

          if (interface->fixedLigand)  strcpy ( fix,"f" );
                                 else  strcpy ( fix," " );
          if (interface->Xrel)  {
            if (fix[0]==' ')  fix[0] = 'x';
                        else  fix[0] = '#';
          }

          sprintf ( S,
            " %3i %2i%1s| %11s | %11s %20s %2i_%1i%1i%1i | "
            "%7.1f %6.1f %3i %3i %3i",
            i+1,interface->type,fix,
            query->getDomainRange(name1,interface->domain1,11),
            query->getDomainRange(name2,interface->domain2,11),
            interface->symOp,symopNo,
            interface->cell_i+nCellOut,
            interface->cell_j+nCellOut,
            interface->cell_k+nCellOut,
            interface->intArea,
            interface->intDeltaG,
            interface->nHBonds,
            interface->nSBridges,
            interface->nDSBonds
                  );

          f.WriteLine ( S );

        }
      }

      f.WriteLine (
       " -------'-------------'----------------------------------"
          "-------'---------------------------\n"
       " ##:  serial number          Monomer: selection range   "
       " Sym.Op: applies to 2nd monomer\n"
       " Id:  interface type id      Area:    sq. angstrom      "
       " DeltaG:  kcal/mol\n"
       " Nhb: no. of hydrogen bonds  Nsb:     no. of salt bridges\n"
       " Nds: no. of disulfide bonds\n"
                  );

    }

    return RESULT_Ok;

  }


  RESULT_CODE Lists::ListInterfaces_csv (  mmdb::io::RFile f ) {
  PInterface interface;
  char       S[1000],name1[20],name2[20],fix[3];
  int        nInterfaces,i, symopNo;

    nInterfaces = query->getNofInterfaces();

    if (nInterfaces<=0)  {
      f.WriteLine ( "NO INTERFACES FOUND\n" );
    } else  {

      f.WriteLine (
        " LIST OF INTERFACES\n\n"
        " \"## Id\",\"Monomer1\",\"Monomer2\",\"Symmetry operation\","
        "\"Sym.Id\",\"Area\",\"DeltaG\",\"Nhb\",\"Nsb\",\"Nds\""
      );

      for (i=0;i<nInterfaces;i++)  {
        interface = query->getInterface ( i );
        if (interface)  {

          if (query->PI->rcsb_symops)  symopNo = interface->rcsb_symop;
                                 else  symopNo = interface->symOpNo+1;

          if (interface->fixedLigand)  strcpy ( fix,"f" );
                                 else  strcpy ( fix," " );
          if (interface->Xrel)  {
            if (fix[0]==' ')  fix[0] = 'x';
                        else  fix[0] = '#';
          }

          sprintf ( S,
            " \"%3i\",\"%2i%1s\",\"%11s\",\"%11s\",\"%20s\","
            "\"%2i_%1i%1i%1i\",\"%7.1f\",\"%6.1f\",\"%3i\",\"%3i\",\"%3i",
            i+1,interface->type,fix,
            query->getDomainRange(name1,interface->domain1,11),
            query->getDomainRange(name2,interface->domain2,11),
            interface->symOp,symopNo,
            interface->cell_i+nCellOut,
            interface->cell_j+nCellOut,
            interface->cell_k+nCellOut,
            interface->intArea,
            interface->intDeltaG,
            interface->nHBonds,
            interface->nSBridges,
            interface->nDSBonds
                  );

          f.WriteLine ( S );

        }

      }

    }

    return RESULT_Ok;

  }


  RESULT_CODE Lists::ListMonomers ( mmdb::cpstr sessionName,
                                    mmdb::cpstr fileName )  {
  mmdb::io::File f;
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

    rc = readPIData();
    if (rc!=RESULT_Ok)  return rc;

    rc = StartTextOutput ( f,fileName );
    if (rc!=RESULT_Ok)  return rc;

    rc = ListMonomers ( f );
    f.shut();

    return rc;

  }

  RESULT_CODE Lists::ListMonomers ( mmdb::io::RFile f )  {
  PDomain domain;
  char    S[1000],name[20],dclass[20],DeltaG[20];
  int     nMonomers,i;

    f.Write (
     " LIST OF MONOMERS\n"
     " -------.-------------.-----------------------------------"
        "-------------\n"
     "  ## Id |   Monomer   |  Class    Nat Nres   Sat Sres     "
        "Area   DeltaG\n"
     " -------+-------------+-----------------------------------"
        "-------------\n"
    );

    nMonomers = query->getNofDomains();
    for (i=0;i<nMonomers;i++)  {
      domain = query->getDomain ( i );
      if (domain)  {

        if (domain->dclass==DCLASS_Protein)  {
          strcpy  ( dclass,"Protein" );
          sprintf ( DeltaG,"%7.1f",domain->DeltaG );
        } else  {
          DeltaG[0] = char(0);
          if (domain->dclass==DCLASS_DNA)      strcpy ( dclass,"DNA" );
          else if (domain->dclass==DCLASS_RNA) strcpy ( dclass,"RNA" );
                                          else strcpy ( dclass,"Ligand " );
        }

        sprintf ( S,
          " %3i %2i | %11s | %7s %5i %4i %5i %4i %9.1f %7s",
          i+1,domain->type,
          query->getDomainRange(name,i,11),dclass,
          domain->nAtoms,
          domain->nRes,
          domain->nSurfAtoms,
          domain->nSurfRes,
          domain->surfArea,
          DeltaG
                );

        f.WriteLine ( S );

      }
    }

    f.WriteLine (
     " -------'-------------'----------------------------------"
        "--------------\n"
     " ##:     serial number     Nat:  total atoms    "
     "Nres: total residues\n"
     " Id:     monomer type id   Sat:  surface atoms  "
     "Sres: surface residues\n"
     " Monomer: selection range  Area: sq. angstrom   "
     "DeltaG: kcal/mol\n"
    );

    return RESULT_Ok;

  }

  RESULT_CODE Lists::ListMonomers_csv ( mmdb::io::RFile f )  {
  PDomain domain;
  char    S[1000],name[20],dclass[20],DeltaG[20];
  int     nMonomers,i;

    f.WriteLine (
      " LIST OF MONOMERS\n\n"
      " \"## Id\",\"Monomer\",\"Class\",\"Nat Nres\",\"Sat Sres\","
      "\"Area\",\"DeltaG\"" );

    nMonomers = query->getNofDomains();
    for (i=0;i<nMonomers;i++)  {
      domain = query->getDomain ( i );
      if (domain)  {

        if (domain->dclass==DCLASS_Protein)  {
          strcpy  ( dclass,"Protein" );
          sprintf ( DeltaG,"%7.1f",domain->DeltaG );
        } else  {
          DeltaG[0] = char(0);
          if (domain->dclass==DCLASS_DNA)      strcpy ( dclass,"DNA" );
          else if (domain->dclass==DCLASS_RNA) strcpy ( dclass,"RNA" );
                                          else strcpy ( dclass,"Ligand " );
        }

        sprintf ( S,
          " \"%3i\",\"%2i\",\"%11s\",\"%7s\",\"%5i\",\"%4i\",\"%5i\","
          "\"%4i\",\"%9.1f\",\"%7s",
          i+1,domain->type,
          query->getDomainRange(name,i,11),dclass,
          domain->nAtoms,
          domain->nRes,
          domain->nSurfAtoms,
          domain->nSurfRes,
          domain->surfArea,
          DeltaG
                );

        f.WriteLine ( S );

      }
    }

    return RESULT_Ok;

  }


  void writePQSTableHeader (  mmdb::io::RFile f, int key )  {

    switch (key)  {

      default:
      case 0: f.Write (  // crystal solutions
  " ----.-----.------------------------------------------------.---------------\n"
  " Set |  No | Size  Id      ASA       BSA   DGdiss0    mG0   | Formula\n"
  " ----+-----+------------------------------------------------+---------------\n"
              );
              break;

      case 1: f.Write (  // as-is complex
  " -----------------------------------------------.---------------\n"
  " Size  Id      ASA       BSA   DGdiss0    mG0   | Formula\n"
  " -----------------------------------------------+---------------\n"
              );
              break;

      case 2: f.Write (  // stock
  " ----.------------------------------------------------.---------------\n"
  "  No | Size  Id      ASA       BSA   DGdiss0    mG0   | Formula\n"
  " ----+------------------------------------------------+---------------\n"
              );

    }

  }


  void writePQSTableHeader_csv (  mmdb::io::RFile f, int key )  {

    switch (key)  {

      default:
      case 0: f.WriteLine (  // crystal solutions
                     " \"Set\",\"No\",\"Size\",\"Id\",\"ASA\",\"BSA\","
                     "\"DGdiss0\",\"mG0\",\"Formula\"" );
              break;

      case 1: f.WriteLine (  // as-is complex
    " \"Size\",\"Id\",\"ASA\",\"BSA\",\"DGdiss0\",\"mG0\",\"Formula\"" );
              break;

      case 2: f.WriteLine (  // stock
                " \"No\",\"Size\",\"Id\",\"ASA\",\"BSA\",\"DGdiss0\","
                "\"mG0\",\"Formula\"" );

    }

  }

  void writePQSTableSeparator (  mmdb::io::RFile f, int key )  {

    switch (key)  {

      default:
      case  0: f.Write (  // crystal solutions
  " ----+-----+------------------------------------------------+---------------\n"
              );
              break;

      case 10: f.Write (  // crystal solutions
  " ----'-----'------------------------------------------------'---------------\n"
              );
              break;

      case  1: f.Write (  // as-is complex
  " -----------------------------------------------+---------------\n"
              );
              break;

      case 11: f.Write (  // as-is complex
  " -----------------------------------------------'---------------\n"
              );
              break;

      case  2: f.Write (  // stock
  " ----+------------------------------------------------+---------------\n"
              );
              break;

      case 12: f.Write (  // stock
  " ----'------------------------------------------------'---------------\n"
              );

    }

  }


  void Lists::writePQSTableAssembly (  mmdb::io::RFile f, PAssembly A,
                                       int setNo, int serNo,
                                       mmdb::pstr & F )  {
  char S[1000];

    if (!A)  return;

    if (setNo>=0)  {
      A->getFormula ( F,query->domains,87,
      "     |     |                                                | " );
      sprintf ( S," %3i |%4i |%4i %4i %9.1f %9.1f %8.1f %8.1f | %s",
                  setNo,serNo,A->mmSize,A->type+1,A->asa,A->bsa,
                  A->freeEn,A->freeEn0/mmdb::IMax(1,A->mmSize),F );
    } else if (serNo>=0)  {

      A->getFormula ( F,query->domains,87,
      "     |                                                | " );
      sprintf ( S,"%4i |%4i %4i %9.1f %9.1f %8.1f %8.1f | %s",
                  serNo,A->mmSize,A->type+1,A->asa,A->bsa,
                  A->freeEn,A->freeEn0/mmdb::IMax(1,A->mmSize),F );
    } else  {
      A->getFormula ( F,query->domains,87,
      "                                                | " );
      sprintf ( S,"%4i %4i %9.1f %9.1f %8.1f %8.1f | %s",
                  A->mmSize,A->type+1,A->asa,A->bsa,
                  A->freeEn,A->freeEn0/mmdb::IMax(1,A->mmSize),F );
    }

    f.WriteLine ( S );

  }


  void Lists::writePQSTableAssembly_csv (  mmdb::io::RFile f,
                                           PAssembly A,
                                           int setNo, int serNo,
                                           mmdb::pstr & F )  {
  char S[1000];

    if (!A)  return;

    if (setNo>=0)  {
      A->getFormula ( F,query->domains,10000,"" );
      sprintf ( S," \"%3i\",\"%4i\",\"%4i\",\"%4i\",\"%9.1f\","
                  "\"%9.1f\",\"%8.1f\",\"%8.1f\",\"%s\"",
                  setNo,serNo,A->mmSize,A->type+1,A->asa,A->bsa,
                  A->freeEn,A->freeEn0/mmdb::IMax(1,A->mmSize),F );
    } else if (serNo>=0)  {

      A->getFormula ( F,query->domains,10000,"" );
      sprintf ( S," \"%4i\",\"%4i\",\"%4i\",\"%9.1f\",\"%9.1f\","
                  "\"%8.1f\",\"%8.1f\",\"%s\"",
                  serNo,A->mmSize,A->type+1,A->asa,A->bsa,
                  A->freeEn,A->freeEn0/mmdb::IMax(1,A->mmSize),F );
    } else  {
      A->getFormula ( F,query->domains,10000,"" );
      sprintf ( S," \"%4i\",\"%4i\",\"%9.1f\",\"%9.1f\",\"%8.1f\","
                  "\"%8.1f\",\"%s\"",
                  A->mmSize,A->type+1,A->asa,A->bsa,
                  A->freeEn,A->freeEn0/mmdb::IMax(1,A->mmSize),F );
    }

    f.WriteLine ( S );

  }


  void Lists::makePQSTable ( mmdb::io::RFile f,
                             int Score1, int Score2 )  {
  PCrystSplit crystSplit;
  PAssembly   A;
  mmdb::pstr  F;
  int         i,j;

    F          = NULL;
    crystSplit = NULL;

    if (Score1==-1)  {
      // printout for "as is: complex is requested

      writePQSTableHeader ( f,1 );
      writePQSTableAssembly ( f,query->Complex,-1,-1,F );
      writePQSTableSeparator ( f,11 );

    } else if (Score1==-2)  {
      // print assembly stock

      writePQSTableHeader ( f,2 );
      if (query->A)  {
        if (query->A->asmStock)
          for (i=0;i<query->A->asmStock->nAssemblies;i++)
            writePQSTableAssembly ( f,query->A->asmStock->A[i],-1,i+1,F );
      }
      writePQSTableSeparator ( f,12 );

    } else  {
      // print crystal solutions

      writePQSTableHeader ( f,0 );

      for (i=0;i<query->A->nCSRes;i++)
        if ((Score1<=query->A->crystSplit[i]->Score) &&
            (query->A->crystSplit[i]->Score<=Score2))  {
          if (crystSplit)
            writePQSTableSeparator ( f,0 );
          crystSplit = query->A->crystSplit[i];
          for (j=0;j<crystSplit->nAssemblies;j++)  {
            A = crystSplit->A[j];
            if (A)
              writePQSTableAssembly ( f,A,i+1,A->serNo,F );
          }
        }

      writePQSTableSeparator ( f,10 );

    }

    if (F) delete[] F;

  }


  void Lists::makePQSTable_csv ( mmdb::io::RFile f,
                                 int Score1, int Score2 )  {
  PCrystSplit crystSplit;
  PAssembly   A;
  mmdb::pstr  F;
  int         i,j;

    F          = NULL;
    crystSplit = NULL;

    if (Score1==-1)  {
      // printout for "as is: complex is requested

      writePQSTableHeader_csv ( f,1 );
      writePQSTableAssembly_csv ( f,query->Complex,-1,-1,F );

    } else if (Score1==-2)  {
      // print assembly stock

      writePQSTableHeader_csv ( f,2 );
      if (query->A)  {
        if (query->A->asmStock)
          for (i=0;i<query->A->asmStock->nAssemblies;i++)
            writePQSTableAssembly_csv ( f,query->A->asmStock->A[i],
                                        -1,i+1,F );
      }

    } else  {
      // print crystal solutions

      writePQSTableHeader_csv ( f,0 );

      for (i=0;i<query->A->nCSRes;i++)
        if ((Score1<=query->A->crystSplit[i]->Score) &&
            (query->A->crystSplit[i]->Score<=Score2))  {
//          if (crystSplit)
//            writePQSTableSeparator ( f,0 );
          crystSplit = query->A->crystSplit[i];
          for (j=0;j<crystSplit->nAssemblies;j++)  {
            A = crystSplit->A[j];
            if (A)
              writePQSTableAssembly_csv ( f,A,i+1,A->serNo,F );
          }
        }

    }

    if (F) delete[] F;

  }


  void Lists::makeAssembliesPage ( mmdb::io::RFile f )  {
  PAssemblies    A;
  mmdb::PResName ligName;
  mmdb::pstr     S,S1;
  int            nLigNames;
  int            i,k,Score1,Score2;

    S  = NULL;
    S1 = NULL;
    A  = query->A;

    ligName   = NULL;
    nLigNames = 0;
    query->getExclLigandList ( ligName,nLigNames );

    if (nLigNames>0)  {
      f.Write ( " In the oligomeric state analysis, ligands\n" );
      k = 0;
      for (i=0;i<nLigNames;i++)  {
        if (!k)  f.Write ( "       " );
        f.Write ( ligName[i] );
        if (i<nLigNames-1)
          f.Write ( ", " );
        k++;
        if (k>18)  {
          f.Write ( "\n" );
          k = 0;
        }
      }
      if (k) f.Write ( "\n" );
      f.Write ( " were excluded.\n\n" );
      delete[] ligName;
    }

    if (A->nCSRes<=0)  {

      // There are no assemblies in the list.

      f.Write (
        " Analysis of protein interfaces has not revealed any specific\n"
        " specific interactions that could result in the formation of\n"
        " stable quaternary structures. Most probably, " );
      if (query->domains->nNCSParents<2)  {
        f.Write ( " structure\n\n                " );
        f.Write ( query->getDomainCode(S,0) );
        f.Write ( "\n\n does " );
      } else {
        f.Write ( " structures\n\n" );
        mmdb::CreateCopy ( S,"    " );
        for (i=0;i<query->domains->nNCSParents;i++)  {
          mmdb::CreateConcat ( S,query->getDomainCode(S1,i) );
          if (i<query->domains->nNCSParents-1)
            mmdb::CreateConcat ( S,", " );
          if (strlen(S)>60)  {
            f.WriteLine ( S );
            mmdb::CreateCopy ( S,"    " );
          }
        }
        if (strlen(S)>4)
          f.WriteLine ( S );
        f.Write ( "\n\n do " );
      }
      f.WriteLine ( "not complexate in solution.\n" );

    } else  {

      if (A->crystSplit[0]->Score>2)  {
        f.Write (
        " Analysis of protein interfaces has not revealed any strong\n"
        " indications that the analysed " );
        if (query->domains->nNCSParents<2)
             f.Write ( "structure" );
        else f.Write ( "structures" );
        f.Write ( " may complexate\n"
        " in solution. The following results show some quaternary\n"
        " structures that can be formed from symmetry and protein\n"
        " affinity considerations.\n" );
      }

      k = 0;
      while (k<A->nCSRes)  {

        Score1 = A->crystSplit[k]->Score;
        if (Score1<=2)  {
          Score1 = 0;
          Score2 = 2;
        } else if (Score1<=5)  {
          Score1 = 3;
          Score2 = 5;
        } else  {
          Score1 = 6;
          Score2 = 7;
        }

        i = 0;
        while (k<A->nCSRes)
          if (A->crystSplit[k]->Score<=Score2)  {
            i += A->crystSplit[k]->nAssemblies;
            k++;
          } else
            break;

        if (Score2<=2)  {
          f.Write (
        " Analysis of protein interfaces suggests that the following\n"
        " quaternary " );
          if (i<2)  f.Write ( "structure is"   );
              else  f.Write ( "structures are" );
          f.Write ( " stable in solution\n" );
        } else if (Score2<=5)  {
          f.Write ( "\n"
        " The following quaternary " );
          if (i<2)  f.Write ( "structure falls" );
              else  f.Write ( "structures fall" );
          f.Write ( " into a grey region\n"
        " of complexation criteria. " );
          if (i<2)  f.Write ( "This structure"   );
              else  f.Write ( "These structures" );
          f.Write ( " may or may not be\n"
        " stable in solution.\n" );
        } else  {
          f.Write ( "\n"
        " The following sets of quaternary " );
          if (i<2)  f.Write ( "structure" );
              else  f.Write ( "structures" );
          f.Write ( " may be formed\n"
        " from crystallographic considerations, however " );
          if (i<2)  f.Write ( "it does" );
              else  f.Write ( "they do" );
          f.Write ( " not\n"
        " appear to be stable in solution.\n" );
        }

        makePQSTable ( f,Score1,Score2 );

      }

    }

    f.LF();

  }


  void Lists::makeAssembliesPage_csv ( mmdb::io::RFile f )  {
  PAssemblies    A;
  mmdb::PResName ligName;
  int            nLigNames;
  int            i,k,Score1,Score2;

    A = query->A;

    ligName   = NULL;
    nLigNames = 0;
    query->getExclLigandList ( ligName,nLigNames );

    if (A->nCSRes<=0)  { // There are no assemblies in the list.
      f.WriteLine ( "No assemblies found" );
    } else  {

      if (A->crystSplit[0]->Score>2)
        f.WriteLine ( "No stable assemblies found" );

      k = 0;
      while (k<A->nCSRes)  {

        Score1 = A->crystSplit[k]->Score;
        if (Score1<=2)  {
          Score1 = 0;
          Score2 = 2;
        } else if (Score1<=5)  {
          Score1 = 3;
          Score2 = 5;
        } else  {
          Score1 = 6;
          Score2 = 7;
        }

        i = 0;
        while (k<A->nCSRes)
          if (A->crystSplit[k]->Score<=Score2)  {
            i += A->crystSplit[k]->nAssemblies;
            k++;
          } else
            break;

        if (Score2<=2)  {
          f.WriteLine ( "\n Stable assemblies:\n" );
        } else if (Score2<=5)  {
          f.WriteLine ( "\n\n Marginally stable assemblies:\n" );
        } else  {
          f.WriteLine ( "\n\n Marginally unstable assemblies:\n" );
        }

        makePQSTable_csv ( f,Score1,Score2 );

      }

    }

    f.LF();

  }


  void Lists::makeComplexAsIsPage (  mmdb::io::RFile f )  {

    if (query->Complex)  {

      f.WriteLine (
      " Analysis of complex, represented by ASU or bare content\n"
      " of the coordinate file" );
      makePQSTable ( f,-1,-1 );

      f.Write (
      " Set: PQS set number   Size: macromolecular oligomeric state\n"
      " Id:  Assembly id      ASA:  accessible surface area, sq.A\n"
      " BSA: buried SA, sq.A\n"
      " DGdiss0: standard free energy of dissociation, kcal/mol\n"
      " mG0:     standard chemical potential, kcal/mol\n"
      "\n"    );

    } else  {

      f.WriteLine ( "\n"
      " Analysis of complex, represented by ASU, was not performed\n"
      " (possibly because it is monomeric).\n"
      "\n" );

    }

  }


  void Lists::makeComplexAsIsPage_csv (  mmdb::io::RFile f )  {

    if (query->Complex)  {

      f.WriteLine (
      " Analysis of complex, represented by ASU or bare content\n"
      " of the coordinate file" );
      makePQSTable_csv ( f,-1,-1 );

    } else  {

      f.WriteLine ( "\n"
      " Analysis of complex, represented by ASU, was not performed\n"
      " (possibly because it is monomeric).\n"
      "\n" );

    }

  }


  RESULT_CODE Lists::ListAssemblies ( mmdb::cpstr sessionName,
                                      mmdb::cpstr fileName ) {
  mmdb::io::File f;
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

    rc = readPIData();
    if (rc!=RESULT_Ok)  return rc;

    rc = readAssemblies();
    if (rc!=RESULT_Ok)  return rc;

    rc = StartTextOutput ( f,fileName );
    if (rc!=RESULT_Ok)  return rc;

    rc = ListAssemblies ( f );
    f.shut();

    return rc;

  }

  RESULT_CODE Lists::ListAssemblies ( mmdb::io::RFile f ) {

    if ((!query->A) || (query->asmStatus!=ASSMB_Ok))
         makeNoAssembliesPage ( f );
    else makeAssembliesPage   ( f );
    makeComplexAsIsPage ( f );

    return RESULT_Ok;

  }

  RESULT_CODE Lists::ListAssemblies_csv ( mmdb::io::RFile f ) {

    if ((!query->A) || (query->asmStatus!=ASSMB_Ok))
         makeNoAssembliesPage ( f );
    else makeAssembliesPage_csv ( f );
    makeComplexAsIsPage_csv ( f );

    return RESULT_Ok;

  }


  void Lists::makeStockPage (  mmdb::io::RFile f )  {
  char S[100];
  int  i,j;

    f.WriteLine ( " Assembly stock" );
    makePQSTable ( f,-2,-2 );

    f.Write (
      " No:  Stock position  Size: macromolecular oligomeric state\n"
      " Id:  Assembly id      ASA:  accessible surface area, sq.A\n"
      " BSA: buried SA, sq.A\n"
      " DGdiss0: standard free energy of dissociation, kcal/mol\n"
      " mG0:     standard chemical potential, kcal/mol\n"
      "\n"  );

    f.WriteLine ( "\n Aggregation index" );
    f.Write ( " ASU Conc " );
    for (i=0;i<query->A->asmStock->nAssemblies;i++)  {
      sprintf ( S, "   Asm %-3i",i+1 );
      f.Write ( S );
    }
    f.LF();
    f.Write ( "----------" );
    for (i=0;i<query->A->asmStock->nAssemblies;i++)
      f.Write ( "----------" );
    f.LF();
    for (j=0;j<query->A->asmStock->nConc;j++)  {
      sprintf ( S,"%10.3g",query->A->asmStock->getRelativeConc(j) );
      f.Write ( S );
      for (i=0;i<query->A->asmStock->nAssemblies;i++)  {
        sprintf ( S,"%10.3g",query->A->asmStock->getAggregationIndex(i,j) );
//        sprintf ( S,"%10.3g",query->A->asmStock->asmConc[i][j] );
        f.Write ( S );
      }
      f.LF();
    }


  }




  RESULT_CODE Lists::ListStock ( mmdb::cpstr sessionName,
                                 mmdb::cpstr fileName ) {
  mmdb::io::File f;
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

    rc = readPIData();
    if (rc!=RESULT_Ok)  return rc;

    rc = readAssemblies();
    if (rc!=RESULT_Ok)  return rc;

    rc = StartTextOutput ( f,fileName );
    if (rc!=RESULT_Ok)  return rc;

    rc = ListStock ( f );
    f.shut();

    return rc;

  }


  RESULT_CODE Lists::ListStock (  mmdb::io::RFile f ) {

    if ((!query->A) || (query->asmStatus!=ASSMB_Ok))
         makeNoAssembliesPage ( f );
    else makeStockPage ( f );

    return RESULT_Ok;

  }

}  // namespace pisa
