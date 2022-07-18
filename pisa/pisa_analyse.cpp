// $Id: pisa_analyse.cpp $
// =================================================================
//
//    20.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_analyse <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Analyse
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2013
//
// =================================================================
//

#include <string.h>

#include "pisa_analyse.h"
//#include "pisalib/pisa_defs.h"

#include "ccp4srs/ccp4srs_defs.h"

namespace pisa  {

  // =======================  Analyse  ==========================

  Analyse::Analyse ( mmdb::cpstr confPath ) : Data ( confPath )  {
    InitAnalyse();
  }

  Analyse::~Analyse()  {
    if (ligExcl)
      delete[] ligExcl;
  }

  void Analyse::InitAnalyse()  {
    ligExcl = NULL;
    ligKey  = LIGANDS_Auto;
  }

  void Analyse::setLigExclude ( mmdb::pstr ligList )  {
    mmdb::CreateCopy ( ligExcl,ligList );
  }

  void Analyse::setLigKey ( LIGAND_KEY ligProc )  {
    ligKey = ligProc;
  }

  RESULT_CODE Analyse::analyse ( mmdb::pstr sessionName,
                                 mmdb::pstr coorFile )  {
  PInterface          interface;
  mmdb::pstr          FN;
  int                 nContacts,i,k;
  RESULT_CODE         rc;
  ccp4srs::CCP4SRS_RC srsrc;
  PROSURF_RC          psrc;
  DATASTATUS          ds;

    //   1. Check configuration

   if (ConfStatus()!=CFG_Configured)
      return RESULT_ConfigurationError;

    rc = readAssemblyParameters();
    if (rc!=RESULT_Ok)
      return rc;


    //   2. mmdb::Create session directory

    switch (makeCrDir(sessionName))  {

      case SDIR_Created : printf ( "\n ... session '%s' Created\n",
                                   sessionName );
                      break;

      case SDIR_Exists  : printf ( "\n ... session '%s' exists, "
                                  "will be reused\n",sessionName );
                      break;

      case SDIR_Fail    : return  RESULT_SessionDirFailure;

      default : ;

    }


    //   3. Load dictionaries

    switch (openDictionaries(false))  {
      case DICT_noSRS         : return  RESULT_noSRS;
      case DICT_noMolRefIndex : return  RESULT_noMolRefIndex;
      case DICT_noMolRefData  : return  RESULT_noMolRefData;
      case DICT_noAgentsData  : return  RESULT_noAgentsData;
      case DICT_Ok :
      default      : ;
    }

    //   4. Prepare surface calculations module

    makeProSurf();


    //   5. Analyse coordinate file

    FN = NULL;

    if (!query)  query = new QueryData();

    query->setLigKey     ( ligKey  );
    query->setLigExclude ( ligExcl );

    ds = query->analyseStructure ( coorFile,molRef,agentsRef );

    if (ds & DS_notAnalysed)           rc = RESULT_notAnalysed;
    else if (ds & DS_noAtoms)          rc = RESULT_noAtoms;
    else if (ds & DS_noChains)         rc = RESULT_noCoordinateData;
    else if (ds & DS_noMacromolecules) rc = RESULT_noMacromolecules;
    else if (ds & DS_badFile)  {
      rc = RESULT_badCoordinateFile;
//      resultCode = ds-DS_badFile;
    } else {

      printf ( " ... structure read and understood\n" );

      srsrc = getEnergyTypes ( query->MMDB );

      if (srsrc==ccp4srs::CCP4SRS_Ok)  {
        // serial numbers will be used to uniquely identify atom
        // in the structure.
        query->MMDB->PDBCleanup    ( mmdb::PDBCLEAN_SERIAL );
        query->calcDimensions      ();
        query->AlignSimilarDomains ();
        query->SortDomains         ();
        query->AlignAllDomains     ( 0,-1 );
        printf ( " ... structural equivalence analysed\n" );

        query->MMDB->DeleteAllSelections();

        // calculate domain surfaces and dimensions for
        // identification of potential interfaces
        psrc = query->CalcSurfaces ( proSurf,molRef );
        if (psrc==PROSURF_Ok)  {

          printf ( " ... molecular surfaces analysed\n" );

          // find all potential contacts. We calculate all potential
          // contacts first and only then run interface calculations,
          // because the contact procedure filters out equivalent
          // contacts.
          rc = query->CalcContacts ( nCellLayers );
          if (rc==RESULT_Ok)  {
            nContacts = query->PI->getNofInterfaces();
            if (nContacts>0)
                  printf ( " ... %i potential contacts identified, "
                           "now processing:\n",nContacts );
            else  printf ( " ... no crystal contacts found\n" );
            k = 0;
            for (i=0;i<nContacts;i++)  {
              interface = query->getInterface ( i );
              if (interface)  {  // calculate the interface
                psrc = interface->calcInterface ( query->MMDB,
                              query->domains->domain[interface->domain1],
                              query->domains->domain[interface->domain2],
                              proSurf,molRef,SRS );
                if (psrc!=PROSURF_Ok)  {
                  printf ( "\n +++ calculations failed for "
                           "contact #%i\n",i+1 );
                  rc = RESULT_InterfaceCalcFailure;
                  k = 0;
                } else  {
                  if (k<1)  printf ( "     " );
                  printf ( "|" );
                  fflush ( stdout );
                  k++;
                  if (k>50) {
                    printf ( "\n" );
                    k = 0;
                  }
                }
              }
            }
            if (k>0)  printf ( "\n" );

            if (rc==RESULT_Ok)  {

              query->checkResidueBSA();
              query->SortInterfaces ();

              if (nContacts>0)  {
                nContacts = query->getNofInterfaces();
                printf ( " ... total %i interfaces calculated\n",
                         nContacts );
              }
              // we proceed with any number of contacts in order to get
              // results for monomeric "as is" assemblies included
              printf ( " ... assembly analysis: " );
              fflush ( stdout );
              query->excludeLigands();
              query->calcAssemblies();
              mmdb::CreateCopCat ( FN,PIStoreDir,rcsb_symops_file );
              rc = query->assignRCSBSymOps ( FN );
              if (rc!=RESULT_Ok)  {
                printResultMessage ( rc );
                rc = RESULT_Ok;
              }

              switch (query->asmStatus)  {
                case ASSMB_Ok           : printf ( "done\n" ); break;
                case ASSMB_noInterfaces : printf (
                     "not performed as no macromolecular interfaces\n"
  "                         were found.\n" );
                                     break;
                case ASSMB_noSymOps     : printf (
                      "not performed as no crystal data were found.\n" );
                                     break;

		//GDL: replacing with 'as is', as there is not assembly analysis performed
		default :  printf ( "assembly 'as is'/ %s\n",
                                    getAsmStatus(query->asmStatus) );
              }

              writeParams();

              if ((writePIData()!=RESULT_Ok) ||
                  (writeStructure()!=RESULT_Ok) ||
                  (writeAssemblies()!=RESULT_Ok))
                rc = RESULT_cantWriteResults;

              if (rc==RESULT_Ok)
                printf ( " ... results written\n" );

            }

          }

        } else if (psrc<0)
          rc = RESULT_SurfaceCalcsFailed;

      } else
        rc = RESULT_SRSFailed;
    }

    if (FN)  delete[] FN;

    return rc;

  }


  #ifdef  _dimer_special

  int Analyse::DimerSpecial ( mmdb::pstr sessionName, mmdb::pstr pdbCode,
                                   mmdb::pstr pdbDir, int & intfNo )  {
  CFile        f;
  PCAssembly   A;
  PPCInterface Interface;
  PPCDomain    Domain;
  ivector      dintf;
  mmdb::pstr         coorFile;
  ChainID      ch1,ch2;
  mat44        tm1;
  realtype     DG0;
  int          rc, nInterfaces, i,k, Nhb,Nsb,Nds, dk1,dk2, di1,di2;

    dintf    = NULL;
    coorFile = NULL;

    //  1. Analyse PDB entry with ligands

    if (pdbDir[strlen(pdbDir)-1]=='/')
         mmdb::CreateCopCat ( coorFile,pdbDir,"pdb",pdbCode,".ent" );
    else mmdb::CreateCopCat ( coorFile,pdbDir,"/pdb",pdbCode,".ent" );

    rc = Analyse ( sessionName,coorFile );
    if (coorFile)  delete[] coorFile;
    if (rc!=RESULT_Ok)  return rc;


    //  2. Check whether the topmost assembly is a dimer

    rc = RESULT_noDimerFound;
    if (query->A->nAsmRes<=0)  {
      printf ( "\n ... no stable assemblies found\n" );
    } else if (query->A->asmSet[0]->Score>2)  {
      printf ( "\n ... only grey-zone assemblies found\n" );
    } else if (!query->A->asmSet[0]->equiv_all)  {
      printf ( "\n ... assembly type mixture\n" );
    } else  {
      A = query->A->asmSet[0]->A[0];
      if (A->mmSize!=2)  {
        printf ( "\n ... no dimers found (mmSize=%i)\n",A->mmSize );
      } else  {
        DG0 = -A->dissEn - A->entropy;
        query->setLigandsAssemble ( False );
        //query->CalcAssemblies();
        writeParams();
        if ((writePIData()!=RESULT_Ok) ||
            (writeStructure()!=RESULT_Ok) ||
            (writeAssemblies()!=RESULT_Ok))  {
           rc = RESULT_cantWriteResults;
        } else if (query->A->nAsmRes<=0)  {
          printf ( "\n ... no stable assemblies found w/o ligands\n" );
        } else if (query->A->asmSet[0]->Score>2)  {
          printf ( "\n ... only grey-zone assemblies found w/o "
                   "ligands\n" );
        } else if (!query->A->asmSet[0]->equiv_all)  {
          printf ( "\n ... assembly type mixture w/o ligands\n" );
        } else  {

          A = query->A->asmSet[0]->A[0];
          if (A->mmSize!=2)  {
            printf ( "\n ... no dimers found w/o ligands\n" );
          } else  {

            nInterfaces = query->getNofInterfaces();
            Interface   = query->getInterfaces   ();
            GetVectorMemory      ( dintf,nInterfaces,0 );
            A->getDissInterfaces ( dintf,nInterfaces   );
            k = -1;
            for (i=0;(i<nInterfaces) && (k<0);i++)
              if (dintf[i]>0)  k = i;

            if (k>=0)  {
              Nhb = Interface[k]->nHBonds;
              Nsb = Interface[k]->nSBridges;
              Nds = Interface[k]->nDSBonds;

              Domain = query->getDomains();
              dk1 = Domain[Interface[k]->domain1]->ncsParent;
              dk2 = Domain[Interface[k]->domain2]->ncsParent;
              query->getChainID ( ch1,dk1 );
              query->getChainID ( ch2,dk2 );
              if (ch2[0]<ch1[0])  {
                query->getChainID ( ch1,dk2 );
                query->getChainID ( ch2,dk1 );
              }

              printf ( "\n"
    " ---------------------------------------------------------------\n"
    " %7.0f %6.0f  %6.1f %5.1f %5.1f %4i %3i %3i   %4s:%1s,%1s\n"
    " ---------------------------------------------------------------\n",
                round_figure(A->asa,-1),round_figure(A->bsa,-1),
                A->seGain,-A->dissEn-A->entropy,DG0,
                Nhb,Nsb,Nds,
                pdbCode,ch1,ch2 );
              printf ( " %4s %1s %1s\n\n",pdbCode,ch1,ch2 );

              if (intfNo<0)  {
                coorFile = new char[1000];
                f.assign ( "_dimer.data",True );
                f.append();
                sprintf ( coorFile,
    " %4i %6.0f %6.0f  %6.1f %5.1f %5.1f %4i %3i %3i   %4s:%1s,%1s",
                  -intfNo,round_figure(A->asa,-1),round_figure(A->bsa,-1),
                  A->seGain,-A->dissEn-A->entropy,DG0,
                  Nhb,Nsb,Nds,pdbCode,ch1,ch2 );
                f.WriteLine ( coorFile );
                f.shut();
                f.assign ( "_dimer.list",True );
                f.append();
                sprintf ( coorFile,"%4s %1s %1s",pdbCode,ch1,ch2 );
                f.WriteLine ( coorFile );
                f.shut();
                delete[] coorFile;
              }

              Mat4Inverse ( Interface[k]->TMatrix,tm1 );
              intfNo = -1;
              for (i=nInterfaces-1;(i>=0) && (intfNo<0);i--)
                if ((i!=k) && (Interface[i]->dclass1!=DCLASS_Ligand) &&
                              (Interface[i]->dclass2!=DCLASS_Ligand))  {
                  di1 = Domain[Interface[i]->domain1]->ncsParent;
                  di2 = Domain[Interface[i]->domain2]->ncsParent;
                  if ((dk1==di1) && (dk2==di2))  {
                    if (!isMat4Eq(Interface[k]->TMatrix,
                                  Interface[i]->TMatrix,
                                  1.0e-4,True))
                      intfNo = i;
                  } else if ((dk1==di2) && (dk2==di1))  {
                    if (!isMat4Eq(tm1,Interface[i]->TMatrix,
                                  1.0e-4,True))
                      intfNo = i;
                  }
                }

              if (intfNo<0)  {
                intfNo = k+1;
                printf ( " ... interface #%i chosen - "
                         "same as the dimer\n",intfNo );
              } else  {
                intfNo++;
                printf ( " ... interface #%i chosen - "
                         "different of the dimer (#%i)\n",intfNo,k+1 );
              }

              rc = RESULT_Ok;

            } else
              printf ( " ++++ no dimer interface found "
                       "(program error)\n" );

            FreeVectorMemory ( dintf,0 );

          }
        }
      }
    }

    return rc;

  }


#endif

}  // namespace pisa
