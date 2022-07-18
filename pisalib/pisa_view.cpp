// $Id: pisa_view.cpp $
// =================================================================
//
//    14.03.19   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_view <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::View
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2019
//
// =================================================================
//

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "pisa_view.h"
#include "pisa_defs.h"
#include "mmdb2/mmdb_manager.h"

namespace pisa  {

  // =========================  View  ===========================

  View::View ( mmdb::cpstr confPath ) : Data ( confPath )  {
    InitPISAView();
  }

  View::~View()  {}

  void View::InitPISAView()  {}

  int  View::launchRasmol ( mmdb::pstr FN )  {
  mmdb::pstr S;
  int  rc;

    S = NULL;
    mmdb::CreateCopCat ( S,rasmol_com," -script ",FN );
    rc = system ( S );
    delete[] S;
    return rc;

  }

  int  View::launchMG ( mmdb::pstr FN )  {
  mmdb::pstr S;
  int  rc;

    S = NULL;
    mmdb::CreateCopCat ( S,ccp4mg_com," -pict ",FN );
    rc = system ( S );
    delete[] S;
    return rc;

  }


  RESULT_CODE View::ViewInput ( mmdb::cpstr   sessionName,
                                VIEWER_TARGET target,
                                mmdb::cpstr   fileName )  {
  mmdb::io::File f;
  mmdb::PManager M;
  mmdb::PPAtom   atom;
  mmdb::pstr     FN;
  int            selHnd,nAtoms;
  int            i;
  RESULT_CODE    rc;

    M = NULL;

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

    rc = readStructure();
    if (rc!=RESULT_Ok)  return rc;

    FN = NULL;

    selHnd = query->MMDB->NewSelection();
    query->MMDB->SelectAtoms ( selHnd,query->MMDB->GetFirstModelNum(),
               "*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
               "*","*","*","*",mmdb::SKEY_NEW );
    query->MMDB->GetSelIndex  ( selHnd,atom,nAtoms );

    if (target==TARGET_Download)  {

      f.assign ( fileName,true,false,mmdb::io::GZM_CHECK );
      f.rewrite();

    } else if (target==TARGET_CIF)  {

      M = new mmdb::Manager();

    } else if (target==TARGET_Rasmol)  {

      f.assign ( makeFName(FN,rasmol_file),true );
      if (!f.rewrite())
        return RESULT_cantWriteRasmolData;

      f.Write ( "load pdb inline\n"
                "select all\n"
                "wireframe off\n"
                "spacefill\n"
                "select atomno=1\n"
                "colour " COLOUR_domain1 "\n"
                "select atomno=2\n"
                "colour " COLOUR_surface1 "\n"
                "select all\n"
                "exit\n" );
    } else  {

      f.assign ( makeFName(FN,ccp4mg_file),true );
      if (!f.rewrite())
        return RESULT_cantWriteMGData;

      f.Write ( "MolData (\n"
                "     filename = ['INLINE', '" );
      f.Write ( sessionID );
      f.Write ( "', '']\n"
                "      )\n"
                "Wizard ( drawing_style = 'PISA_data:interface' )\n"
                "\n"
                "Inline (name =  '" );
      f.Write ( sessionID );
      f.Write ( "' , data = '''\n" );

    }

    for (i=0;i<nAtoms;i++)  {
      if (target==TARGET_CIF)  M->PutAtom ( 0,atom[i] );
                         else  atom[i]->PDBASCIIDump ( f );
    }

    if (target==TARGET_CIF)  {
      M->WriteCIFASCII ( fileName,mmdb::io::GZM_CHECK );
    } else if (target==TARGET_CCP4MG)  {
      f.Write ( "''' )\n" );
      f.shut();
      launchMG ( FN );
    } else  {
      f.shut();
      if (target==TARGET_Rasmol)
        launchRasmol ( FN );
    }

    if (FN)  delete[] FN;
    if (M)   delete   M;

    return rc;

  }


  RESULT_CODE View::ViewMonomer ( mmdb::cpstr   sessionName,
                                  int           serialNo,
                                  VIEWER_TARGET target,
                                  mmdb::cpstr   fileName )  {
  mmdb::io::File f;
  mmdb::PManager M;
  mmdb::PPAtom   atom;
  PDomain        Domain;
  PQueryData     query0;
  mmdb::pstr     FN;
  mmdb::realtype x,y,z;
  int            selHnd,nAtoms;
  int            i,sN;
  RESULT_CODE    rc;


    //   1. Check configuration

    if (ConfStatus()!=CFG_Configured)
      return RESULT_ConfigurationError;

    query0 = query;
    query  = NULL;
    M      = NULL;

    //   2. Check session directory and results

    switch (checkCrDir(sessionName))  {
      case SDIR_doesntExist : query = query0;
                            return RESULT_SessionDoesntExist;
      case SDIR_noResults   : query = query0;
                            return RESULT_noSessionResults;
      default : ;
    }

    rc = readPIData();
    if (rc!=RESULT_Ok)  {
      if (query)  delete query;
      query = query0;
      return rc;
    }

    if ((serialNo<1) || (serialNo>query->getNofDomains()))  {
      if (query)  delete query;
      query = query0;
      return RESULT_monomerNoOutOfRange;
    }

    rc = readStructure();
    if (rc!=RESULT_Ok)  {
      if (query)  delete query;
      query = query0;
      return rc;
    }

    FN = NULL;

    Domain = query->getDomain ( serialNo-1 );
    if (!Domain)  {
      if (query)  delete query;
      query = query0;
      return RESULT_monomerNoOutOfRange;
    }

    // cut out domain from resMMDB
    selHnd = Domain->SelectDomain ( query->MMDB,mmdb::STYPE_ATOM,1,false );
    query->MMDB->GetSelIndex  ( selHnd,atom,nAtoms );

    if (target==TARGET_Download)  {

      f.assign ( fileName,true,false,mmdb::io::GZM_CHECK );
      f.rewrite();

    } else if (target==TARGET_CIF)  {

      M = new mmdb::Manager();

    } else if (target==TARGET_Rasmol)  {

      f.assign ( makeFName(FN,rasmol_file),true );
      if (!f.rewrite())  {
        if (query)  delete query;
        query = query0;
        return RESULT_cantWriteRasmolData;
      }

      f.Write ( "load pdb inline\n"
                "select all\n"
                "wireframe off\n"
                "spacefill\n"
                "select atomno=1\n"
                "colour " COLOUR_domain1 "\n"
                "select atomno=2\n"
                "colour " COLOUR_surface1 "\n"
                "select all\n"
                "exit\n" );
    } else  {

      f.assign ( makeFName(FN,ccp4mg_file),true );
      if (!f.rewrite())  {
        if (query)  delete query;
        query = query0;
        return RESULT_cantWriteMGData;
      }

      f.Write ( "MolData (\n"
                "     filename = ['INLINE', '" );
      f.Write ( sessionID );
      f.Write ( "', '']\n"
                "      )\n"
                "Wizard ( drawing_style = 'PISA_data:interface' )\n"
                "\n"
                "Inline (name =  '" );
      f.Write ( sessionID );
      f.Write ( "' , data = '''\n" );

    }

    for (i=0;i<nAtoms;i++)  {
      sN = atom[i]->serNum;  // make consequitive serial numbers as
      if (target!=TARGET_Download)  {
        if (atom[i]->isInSelection(Domain->selHndSurf))
              atom[i]->serNum = 2;
        else  atom[i]->serNum = 1;
      }
      x = atom[i]->x;
      y = atom[i]->y;
      z = atom[i]->z;
      atom[i]->Transform ( Domain->ncs_m );
      if (target==TARGET_CIF)  M->PutAtom ( 0,atom[i] );
                         else  atom[i]->PDBASCIIDump ( f );
      atom[i]->x = x;
      atom[i]->y = y;
      atom[i]->z = z;
      atom[i]->serNum = sN;  // restore serial number
    }

    if (target==TARGET_CIF)  {
      M->WriteCIFASCII ( fileName,mmdb::io::GZM_CHECK );
    } else if (target==TARGET_CCP4MG)  {
      f.Write ( "''' )\n" );
      f.shut();
      launchMG ( FN );
    } else  {
      f.shut();
      if (target==TARGET_Rasmol)
        launchRasmol ( FN );
    }

    if (FN)     delete[] FN;
    if (M)      delete   M;
    if (query)  delete   query;
    query = query0;

    return rc;

  }


  RESULT_CODE View::ViewInterface ( mmdb::cpstr   sessionName,
                                    int           serialNo,
                                    VIEWER_TARGET target,
                                    mmdb::cpstr   fileName )  {
  mmdb::io::File f;
  mmdb::PManager M;
  mmdb::PPAtom   atom;
  mmdb::PPChain  chain1,chain2;
  PQueryData     query0;
  PDomain        PID;
  PInterface     interface;
  mmdb::PChainID chID1,chID2;
  mmdb::pstr     FN;
  mmdb::ChainID  chID;
  mmdb::mat44    tm;
  mmdb::realtype x,y,z;
  int            selHnd,selHndChain,nAtoms,nChains1,nChains2;
  int            i,j,sN;
  RESULT_CODE    rc;
  bool           B;

    query0 = query;
    query  = NULL;
    M      = NULL;

    //   1. Check configuration

    if (ConfStatus()!=CFG_Configured)  {
      if (query)  delete query;
      query = query0;
      return RESULT_ConfigurationError;
    }

    //   2. Check session directory and results

    switch (checkCrDir(sessionName))  {

      case SDIR_doesntExist : return RESULT_SessionDoesntExist;
      case SDIR_noResults   : return RESULT_noSessionResults;

      default : ;

    }

    rc = readPIData();
    if (rc!=RESULT_Ok)  {
      if (query)  delete query;
      query = query0;
      return rc;
    }

    if ((serialNo<1) || (serialNo>query->getNofInterfaces()))  {
      if (query)  delete query;
      query = query0;
      return RESULT_interfaceNoOutOfRange;
    }

    rc = readStructure();
    if (rc!=RESULT_Ok)  {
      if (query)  delete query;
      query = query0;
      return rc;
    }

    FN = NULL;

    interface = query->getInterface ( serialNo-1 );
    if (!interface)  {
      if (query)  delete query;
      query = query0;
      return RESULT_interfaceNoOutOfRange;
    }

    // output 1st domain
    PID    = query->getDomain  ( interface->domain1 );
    selHnd = PID->SelectDomain ( query->MMDB,mmdb::STYPE_ATOM,1,false );
    query->MMDB->GetSelIndex   ( selHnd,atom,nAtoms );

    if (target==TARGET_Download)  {

      f.assign ( fileName,true,false,mmdb::io::GZM_CHECK );
      f.rewrite();

    } else if (target==TARGET_CIF)  {

      M = new mmdb::Manager();

    } else if (target==TARGET_Rasmol)  {

      f.assign ( makeFName(FN,rasmol_file),true );
      if (!f.rewrite())
        return RESULT_cantWriteRasmolData;

      f.Write ( "load pdb inline\n"
                "select all\n"
                "wireframe off\n"
                "spacefill\n"
                "select atomno=1\n"
                "colour " COLOUR_domain1 "\n"
                "select atomno=2\n"
                "colour " COLOUR_surface1 "\n"
                "select atomno=3\n"
                "colour " COLOUR_interface1 "\n"
                "select atomno=4\n"
                "colour " COLOUR_domain2 "\n"
                "select atomno=5\n"
                "colour " COLOUR_surface2 "\n"
                "select atomno=6\n"
                "colour " COLOUR_interface2 "\n"
                "select all\n"
                "exit\n" );
    } else  {

      f.assign ( makeFName(FN,ccp4mg_file),true );
      if (!f.rewrite())
        return RESULT_cantWriteMGData;

      f.Write ( "MolData (\n"
                "     filename = ['INLINE', '" );
      f.Write ( sessionID );
      f.Write ( "', '']\n"
                "      )\n"
                "Wizard ( drawing_style = 'PISA_data:interface' )\n"
                "\n"
                "Inline (name =  '" );
      f.Write ( sessionID );
      f.Write ( "' , data = '''\n" );

    }

    for (i=0;i<nAtoms;i++)  {
      sN = atom[i]->serNum;
      if (target!=TARGET_Download)  {
        if (atom[i]->isInSelection(interface->selHndInt1))
             atom[i]->serNum = 3;
        else if (atom[i]->isInSelection(PID->selHndSurf))
             atom[i]->serNum = 2;
        else atom[i]->serNum = 1;
      }
      x = atom[i]->x;
      y = atom[i]->y;
      z = atom[i]->z;
      atom[i]->Transform ( PID->ncs_m );
      if (target==TARGET_CIF)  M->PutAtom ( 0,atom[i] );
                         else  atom[i]->PDBASCIIDump ( f );
      atom[i]->x = x;
      atom[i]->y = y;
      atom[i]->z = z;
      atom[i]->serNum = sN;
    }

    selHndChain = query->MMDB->NewSelection();
    query->MMDB->Select ( selHndChain,mmdb::STYPE_CHAIN,
                          selHnd,mmdb::SKEY_NEW );
    query->MMDB->GetSelIndex ( selHndChain,chain1,nChains1 );
    chID1 = new mmdb::ChainID[nChains1];
    for (i=0;i<nChains1;i++)
      strcpy ( chID1[i],chain1[i]->GetChainID() );

    // output 2nd domain
    PID    = query->getDomain  ( interface->domain2 );
    selHnd = PID->SelectDomain ( query->MMDB,mmdb::STYPE_ATOM,1,false );
    query->MMDB->GetSelIndex   ( selHnd,atom,nAtoms );

    query->MMDB->Select ( selHndChain,mmdb::STYPE_CHAIN,
                          selHnd,mmdb::SKEY_NEW );
    query->MMDB->GetSelIndex ( selHndChain,chain2,nChains2 );
    chID2 = new mmdb::ChainID[nChains2];
    for (i=0;i<nChains2;i++)
      strcpy ( chID2[i],chain2[i]->GetChainID() );

    strcpy ( chID,"A" );
    for (i=0;i<nChains2;i++)
      do  {
        B = false;
        for (j=0;(j<i) && (!B);j++)
          B = !strcmp(chain2[i]->GetChainID(),chain2[j]->GetChainID());
        for (j=0;(j<nChains1) && (!B);j++)
          B = !strcmp(chain2[i]->GetChainID(),chID1[j]);
        if (B) {
          chID[0] = char(int(chID[0])+1);
          chain2[i]->SetChainID ( chID );
        }
      } while (B);

    mmdb::Mat4Mult ( tm,interface->TMatrix,PID->ncs_m );
    for (i=0;i<nAtoms;i++)  {
      sN = atom[i]->serNum;
      if (target!=TARGET_Download)  {
        if (atom[i]->isInSelection(interface->selHndInt2))
             atom[i]->serNum = 6;
        else if (atom[i]->isInSelection(PID->selHndSurf))
             atom[i]->serNum = 5;
        else atom[i]->serNum = 4;
      }
      x = atom[i]->x;
      y = atom[i]->y;
      z = atom[i]->z;
      atom[i]->Transform ( tm );
      if (target==TARGET_CIF)  M->PutAtom ( 0,atom[i] );
                         else  atom[i]->PDBASCIIDump ( f );
      atom[i]->x = x;
      atom[i]->y = y;
      atom[i]->z = z;
      atom[i]->serNum = sN;
    }

    for (i=0;i<nChains2;i++)
      chain2[i]->SetChainID ( chID2[i] );

    delete[] chID1;
    delete[] chID2;

    query->MMDB->DeleteSelection ( selHnd      );
    query->MMDB->DeleteSelection ( selHndChain );

    if (target==TARGET_CIF)  {
      M->WriteCIFASCII ( fileName,mmdb::io::GZM_CHECK );
    } else if (target==TARGET_CCP4MG)  {
      f.Write ( "''' )\n" );
      f.shut();
      launchMG ( FN );
    } else  {
      f.shut();
      if (target==TARGET_Rasmol)
        launchRasmol ( FN );
    }

    if (FN)     delete[] FN;
    if (M)      delete   M;
    if (query)  delete   query;
    query = query0;

    return rc;

  }


  RESULT_CODE View::ViewAssembly ( mmdb::cpstr   sessionName,
                                   int           serialNo,
                                   bool          dissociated,
                                   VIEWER_TARGET target,
                                   mmdb::cpstr   fileName )  {
  mmdb::io::File  f;
  mmdb::PManager  M;
  mmdb::PPAtom    atom;
  PQueryData      query0;
  PAssembly       A;
  mmdb::rvector   dx,dy,dz, wt;
  mmdb::ivector   selHnd;
  mmdb::pstr      FN;
  mmdb::ChainID   chnId,chId0;
  mmdb::realtype  x,y,z, w, x1,x2, y1,y2, z1,z2, R;
  int             nAtoms,i,j,k,dNo;
  RESULT_CODE     rc;

    query0 = query;
    query  = NULL;
    M      = NULL;

    //   1. Check configuration

    if (ConfStatus()!=CFG_Configured)  {
      if (query)  delete query;
      query = query0;
      return RESULT_ConfigurationError;
    }

    //   2. Check session directory and results

    switch (checkCrDir(sessionName))  {

      case SDIR_doesntExist : return RESULT_SessionDoesntExist;
      case SDIR_noResults   : return RESULT_noSessionResults;

      default : ;

    }

    rc = readPIData();
    if (rc!=RESULT_Ok)  {
      if (query)  delete query;
      query = query0;
      return rc;
    }

    rc = readAssemblies();
    if (rc!=RESULT_Ok)  {
      if (query)  delete query;
      query = query0;
      return rc;
    }

    if (serialNo>0)  {
      if ((!query->A) || (query->asmStatus!=ASSMB_Ok))  {
        f.assign ( fileName,true,false,mmdb::io::GZM_CHECK );
        f.rewrite();
        makeNoAssembliesPage ( f );
        f.shut();
        if (query)  delete query;
        query = query0;
        return RESULT_assemblyNoOutOfRange;
      }

      A = NULL;
      k = 0;
      for (i=0;(i<query->A->nCSRes) && (!A);i++)
        for (j=0;(j<query->A->crystSplit[i]->nAssemblies) && (!A);j++)
          if (query->A->crystSplit[i]->A[j]->serNo==serialNo)
             A = query->A->crystSplit[i]->A[j];

      if (!A)  {
        if (query)  delete query;
        query = query0;
        return RESULT_assemblyNoOutOfRange;
      }

    } else  {
      A = query->Complex;
      if ((!A) || (query->asmStatus!=ASSMB_Ok))  {
        f.assign ( fileName,true,false,mmdb::io::GZM_CHECK );
        f.rewrite();
        makeNoAssembliesPage ( f );
        f.shut();
        if (query)  delete query;
        query = query0;
        return RESULT_assemblyNoOutOfRange;
      }
    }

    rc = readStructure();
    if (rc!=RESULT_Ok)  {
      if (query)  delete query;
      query = query0;
      return rc;
    }

    FN = NULL;

    if (target==TARGET_Download)  {

      f.assign ( fileName,true,false,mmdb::io::GZM_CHECK );
      f.rewrite();

    } else if (target==TARGET_CIF)  {

      M = new mmdb::Manager();

    } else if (target==TARGET_Rasmol)  {

      f.assign ( makeFName(FN,rasmol_file),true );
      if (!f.rewrite())
        return RESULT_cantWriteRasmolData;

      f.Write ( "load pdb inline\n"
                "select all\n"
                "wireframe off\n"
                "spacefill\n"
                "colour chains\n"
                "exit\n" );

    } else  {

      f.assign ( makeFName(FN,ccp4mg_file),true );
      if (!f.rewrite())
        return RESULT_cantWriteMGData;

      f.Write ( "MolData (\n"
                "     filename = ['INLINE', '" );
      f.Write ( sessionID );
      f.Write ( "', '']\n"
                "      )\n"
                "Wizard ( drawing_style = 'PISA_data:assembly' )\n"
                "\n"
                "Inline (name =  '" );
      f.Write ( sessionID );
      f.Write ( "' , data = '''\n" );

    }


    mmdb::GetVectorMemory ( selHnd,query->domains->nDomains,0 );
    for (i=0;i<query->domains->nDomains;i++)
      selHnd[i] = 0;

    if (dissociated)  {
      mmdb::GetVectorMemory ( dx,A->nDiss,1 );
      mmdb::GetVectorMemory ( dy,A->nDiss,1 );
      mmdb::GetVectorMemory ( dz,A->nDiss,1 );
      mmdb::GetVectorMemory ( wt,A->nDiss,1 );
      for (i=1;i<=A->nDiss;i++)  {
        dx[i] = 0.0;   // mass centers of dissociation units
        dy[i] = 0.0;
        dz[i] = 0.0;
        wt[i] = 0.0;   // weights of dissociation units
      }
      x1 =  mmdb::MaxReal;
      x2 = -mmdb::MaxReal;
      y1 =  mmdb::MaxReal;
      y2 = -mmdb::MaxReal;
      z1 =  mmdb::MaxReal;
      z2 = -mmdb::MaxReal;
      for (i=0;i<A->asmSize;i++)  {
        dNo = A->M[i]->id;
        j   = A->M[i]->dissAsm;
        k   = query->domains->domain[dNo]->ncsParent;
        x   = query->domains->domain[k]->mx;
        y   = query->domains->domain[k]->my;
        z   = query->domains->domain[k]->mz;
        w   = query->domains->domain[k]->weight;
        mmdb::TransformXYZ ( A->M[i]->T, x,y,z );
        dx[j] += w*x;
        dy[j] += w*y;
        dz[j] += w*z;
        wt[j] += w;
        x1 = mmdb::RMin ( x1,x );
        x2 = mmdb::RMax ( x2,x );
        y1 = mmdb::RMin ( y1,y );
        y2 = mmdb::RMax ( y2,y );
        z1 = mmdb::RMin ( z1,z );
        z2 = mmdb::RMax ( z2,z );
      }
      x = 0.0;  // assembly mass center
      y = 0.0;
      z = 0.0;
      w = 0.0;  // assembly weight
      for (i=1;i<=A->nDiss;i++)
        if (wt[i]>0.0)  {
          x += dx[i];
          y += dy[i];
          z += dz[i];
          w += wt[i];
          dx[i] /= wt[i];
          dy[i] /= wt[i];
          dz[i] /= wt[i];
        }
      if (w>0.0)  {
        x /= w;
        y /= w;
        z /= w;
        x2 -= x1;
        y2 -= y1;
        z2 -= z1;
        R  = mmdb::RMax ( query->domains->Dmax,sqrt(x2*x2+y2*y2+z2*z2) );
        R /= 5.0;
        for (i=1;i<=A->nDiss;i++)
          if (wt[i]>0.0)  {
            dx[i] -= x;
            dy[i] -= y;
            dz[i] -= z;
            w = sqrt ( dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i] );
            if (w>0.2)  {
              w  = R/w;
              dx[i] *= w;
              dy[i] *= w;
              dz[i] *= w;
            } else  {
              dx[i] = 0.0;
              dy[i] = 0.0;
              dz[i] = 0.0;
            }
          }
      }
    }

    strcpy ( chnId,"A" );
    for (i=0;i<A->asmSize;i++)  {
      dNo = A->M[i]->id;
      if (!selHnd[dNo])
        selHnd[dNo] = query->domains->domain[dNo]->SelectDomain (
                                   query->MMDB,mmdb::STYPE_ATOM,1,false );
      query->MMDB->GetSelIndex ( selHnd[dNo],atom,nAtoms );
      if (nAtoms>0)  {
        strcpy ( chId0,atom[0]->GetChain()->GetChainID() );
        atom[0]->GetChain()->SetChainID ( chnId );
        chnId[0]++;
        if (chnId[0]>'Z')
          chnId[0] = 'A';
      }
      k = A->M[i]->dissAsm;
      for (j=0;j<nAtoms;j++)  {
        x = atom[j]->x;
        y = atom[j]->y;
        z = atom[j]->z;
        atom[j]->Transform ( A->M[i]->T );
        if (dissociated)  {
          atom[j]->x += dx[k];
          atom[j]->y += dy[k];
          atom[j]->z += dz[k];
        }
        if (target==TARGET_CIF)  M->PutAtom ( 0,atom[j] );
                           else  atom[j]->PDBASCIIDump ( f );
        atom[j]->x = x;
        atom[j]->y = y;
        atom[j]->z = z;
      }
      if (nAtoms>0)
        atom[0]->GetChain()->SetChainID ( chId0 );
    }

    if (dissociated)  {
      mmdb::FreeVectorMemory ( dx,1 );
      mmdb::FreeVectorMemory ( dy,1 );
      mmdb::FreeVectorMemory ( dz,1 );
      mmdb::FreeVectorMemory ( wt,1 );
    }

    for (i=0;i<query->domains->nDomains;i++)
      if (selHnd[i])
        query->MMDB->DeleteSelection ( selHnd[i] );

    mmdb::FreeVectorMemory ( selHnd,0 );

    if (target==TARGET_CIF)  {
      M->WriteCIFASCII ( fileName,mmdb::io::GZM_CHECK );
    } else if (target==TARGET_CCP4MG)  {
      f.Write ( "''' )\n" );
      f.shut();
      launchMG ( FN );
    } else  {
      f.shut();
      if (target==TARGET_Rasmol)
        launchRasmol ( FN );
    }

    if (FN)     delete[] FN;
    if (M)      delete   M;
    if (query)  delete   query;
    query = query0;

    return rc;

  }

}  // namespace pisa;
