// $Id: pisa_molref.cpp $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_molref <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::AtmRefData
//       ~~~~~~~~~  pisa::MolRefData
//                  pisa::MolRefIndex
//
//  (C) E. Krissinel, 2006-2014
//
// =================================================================
//

#include <string.h>

#include "mmdb2/mmdb_atom.h"
#include "mmdb2/mmdb_tables.h"

#include "pisa_prosurf.h"
#include "pisa_molref.h"

#include "ccp4srs/ccp4srs_manager.h"
#include "ccp4srs/ccp4srs_defs.h"

namespace pisa  {

  mmdb::cpstr molref_index_file = "molref.idx";
  mmdb::cpstr molref_data_file  = "molref.rdt";


  // =========================  AtmRefData  =========================

  int AtmRefData::readRefData ( mmdb::mmcif::PLoop mmCIFLoop, int aNo )  {
  mmdb::pstr p;
  int        rc,rc0;

    rc0 = 0;
    p = mmCIFLoop->GetString ( "pdb_atom_name",aNo,rc );
    if (!p)  {
      name[0] = char(0);
      rc0 = 1;
    } else
      strncpy ( name,p,sizeof(name) );

    p = mmCIFLoop->GetString ( "pdb_atom_name_old",aNo,rc );
    if (p)  strncpy ( old_name,p,sizeof(old_name) );
      else  strcpy  ( old_name,name );

    if (mmCIFLoop->GetReal(radius,"radius"  ,aNo,false))  rc0 = 2;
    if (mmCIFLoop->GetReal(refASA,"ref_asa" ,aNo,false))  rc0 = 3;
    if (mmCIFLoop->GetReal(asp   ,"asp"     ,aNo,false))  rc0 = 4;
    if (mmCIFLoop->GetInteger(aspId,"asp_id",aNo,false))  rc0 = 5;

    return rc0;

  }

  void AtmRefData::allocRefData ( mmdb::pstr atname,
                                  mmdb::pstr old_atname )  {
    strncpy ( name,atname,sizeof(name) );
    name[4] = char(0);
    strncpy ( old_name,old_atname,sizeof(old_name) );
    old_name[4] = char(0);
    radius  = 0.0;
    refASA  = 0.0;
    asp     = 0.0;
    aspId   = 0;
  }

  void AtmRefData::swap ( RAtmRefData a )  {
  char           S[sizeof(name)];
  mmdb::realtype V;
  int            i;

    strcpy ( S,name      );
    strcpy ( name,a.name );
    strcpy ( a.name,S    );

    strcpy ( S,old_name          );
    strcpy ( old_name,a.old_name );
    strcpy ( a.old_name,S        );

    V = radius;   radius = a.radius;    a.radius = V;
    V = refASA;   refASA = a.refASA;    a.refASA = V;
    V = asp;      asp    = a.asp;       a.asp    = V;
    i = aspId;    aspId  = a.aspId;     a.aspId  = i;

  }

  void AtmRefData::write ( mmdb::io::RFile f )  {
    f.WriteFile  ( name    ,sizeof(name)     );
    f.WriteFile  ( old_name,sizeof(old_name) );
    f.WriteFloat ( &radius );
    f.WriteFloat ( &refASA );
    f.WriteFloat ( &asp    );
    f.WriteInt   ( &aspId  );
  }

  void AtmRefData::read ( mmdb::io::RFile f )  {
    f.ReadFile  ( name    ,sizeof(name)     );
    f.ReadFile  ( old_name,sizeof(old_name) );
    f.ReadFloat ( &radius );
    f.ReadFloat ( &refASA );
    f.ReadFloat ( &asp    );
    f.ReadInt   ( &aspId  );
  }


  // =========================  MolRefData  =========================

  MolRefData::MolRefData()  {
    InitMolRefData();
  }

  MolRefData::~MolRefData()  {
    FreeMemory();
  }

  void MolRefData::InitMolRefData()  {
    name[0]    = char(0);
    refName[0] = char(0);
    classId    = 0;
    symNumber  = 1;
    nAtoms     = 0;
    atom       = NULL;
  }

  void MolRefData::FreeMemory()  {
    if (atom)  delete[] atom;
    atom   = NULL;
    nAtoms = 0;
  }

  int MolRefData::readRefData ( mmdb::mmcif::PData mmCIF )  {
  mmdb::mmcif::PLoop   mmCIFLoop;
  mmdb::mmcif::PStruct mmCIFStruct;
  mmdb::pstr           p;
  int                  i,rc;

    FreeMemory();

    mmCIFLoop = mmCIF->GetLoop ( "_ref_data" );
    if (!mmCIFLoop)  return -1;

    strncpy ( name,mmCIF->GetDataName(),sizeof(name) );
    name[sizeof(name)-1] = char(0);
    i = 0;
    while (name[i])  {
      if (name[i]=='_')  name[i] = ' ';
      i++;
    }

    mmCIFStruct = mmCIF->GetStructure ( "_classification" );
    if (mmCIFStruct)  {
      p = mmCIFStruct->GetString ( "ref_name",rc );
      if (p)  strncpy ( refName,p,sizeof(refName) );
        else  refName[0] = char(0);
      if (mmCIFStruct->GetInteger(classId,"id",false))
        classId = MRCLASS_Unknown;
      if (mmCIFStruct->GetInteger(symNumber,"sym_number",false))
        symNumber = 1;
    } else  {
      refName[0] = char(0);
      classId    = MRCLASS_Unknown;
      symNumber  = 1;
    }


    nAtoms = mmCIFLoop ->GetLoopLength();

    if (nAtoms>0)  {

      rc = 0;

      atom = new AtmRefData[nAtoms];
      for (i=0;i<nAtoms;i++)
        if (atom[i].readRefData(mmCIFLoop,i))  {
          printf ( " *** problems in atom #%i structure '%s'\n",
                   i,mmCIF->GetDataName() );
          rc = 1;
        }

      SortAtoms();

      return rc;

    } else
      rc = -1;

    return rc;

  }

  void alignName ( mmdb::cpstr molName, mmdb::ResName resName )  {
  int l1,l2;
    l1 = strlen(molName);
    l2 = 2;
    while (l2>=0)  {
      l1--;
      if (l1>=0) resName[l2--] = molName[l1];
            else resName[l2--] = ' ';
    }
    resName[3] = char(0);
  }

  int MolRefData::allocRefData ( ccp4srs::PMonomer monomer )  {
  ccp4srs::PAtom atomi;
  mmdb::AtomName anamei;
  int            i,j,rc, mon_natoms;

    FreeMemory();

    alignName ( monomer->ID(),name );
    strcpy    ( refName,name );
    classId = MRCLASS_Ligand;

    mon_natoms = monomer->n_atoms();
    nAtoms     = 0;
    for (i=0;i<mon_natoms;i++)  {
      atomi = monomer->atom(i);
      if (atomi)  {
        atomi->name_pdb ( anamei );
        if (anamei[1]!='H')  nAtoms++;
        if (!strcmp(anamei," OXT"))
          classId = MRCLASS_ModAminoacid;
      }
    }

    if (mon_natoms<3)  symNumber = 4;
                 else  symNumber = 1;

    if (nAtoms>0)  {

      rc = 0;

      atom = new AtmRefData[nAtoms];
      j    = 0;
      for (i=0;i<mon_natoms;i++)  {
        atomi = monomer->atom(i);
        if (atomi)  {
          atomi->name_pdb ( anamei );
          if (anamei[1]!='H')
            atom[j++].allocRefData ( anamei,anamei ); // 2nd should be
                                     // the old PDB atom name
        }
      }

    } else
      rc = -1;

    return rc;

  }


  void MolRefData::SortAtoms()  {
  int i,j;

    for (i=0;i<nAtoms;i++)
      for (j=i+1;j<nAtoms;j++)
        if (strcmp(atom[i].name,atom[j].name)>0)
          atom[i].swap ( atom[j] );

  }


  void MolRefData::CheckRepeats()  {
  int i,j;

    for (i=0;i<nAtoms;i++)
      for (j=i+1;j<nAtoms;j++)
        if (!strcmp(atom[i].name,atom[j].name))
          printf ( " *** atom %s in compound %s is entered more "
                   "than once.\n",atom[i].name,name );

  }

  void MolRefData::CheckASPRange()  {
  int i;
    for (i=0;i<nAtoms;i++)
      if ((atom[i].aspId>=nASPs) ||
          (atom[i].aspId<-2) || (atom[i].aspId==-1))
        printf ( " *** atom %s in compound %s has wrong ASP Id %i\n",
                 atom[i].name,name,atom[i].aspId );
  }


  PAtmRefData  MolRefData::getAtmRefData ( mmdb::pstr atmName )  {
  int     l1,l2,l,k;

    l1 = 0;
    l2 = nAtoms-1;
    do  {
      l = (l1+l2)/2;
      k = strcmp ( atmName,atom[l].name );
      if (k<0)      l2 = l;
      else if (k>0) l1 = l;
               else  break;
    } while (l1<l2-1);

    if (k!=0)  {
      if (l==l1)  l = l2;
            else  l = l1;
      k = strcmp ( atmName,atom[l].name );
      if (k)  {
        // Search old names. They are not alphabetically ordered,
        // so the search is a primitive look-up.
        k = 0;
        for (l=0;(l<nAtoms) && (!k);l++)
          if (!strcmp(atmName,atom[l].old_name))
            k = l + 1;
        if (!k)  return NULL;
        l = k - 1;
      }
    }

    return &(atom[l]);

  }

  void MolRefData::write ( mmdb::io::RFile f )  {
  int i;
    f.WriteFile ( name   ,sizeof(name)    );
    f.WriteFile ( refName,sizeof(refName) );
    f.WriteInt  ( &classId   );
    f.WriteInt  ( &symNumber );
    f.WriteInt  ( &nAtoms    );
    for (i=0;i<nAtoms;i++)
      atom[i].write ( f );
  }

  void MolRefData::read ( mmdb::io::RFile f )  {
  int i;
    FreeMemory();
    f.ReadFile ( name   ,sizeof(name)    );
    f.ReadFile ( refName,sizeof(refName) );
    f.ReadInt  ( &classId   );
    f.ReadInt  ( &symNumber );
    f.ReadInt  ( &nAtoms    );
    if (nAtoms>0)  {
      atom = new AtmRefData[nAtoms];
      for (i=0;i<nAtoms;i++)
        atom[i].read ( f );
    }
  }


  // =========================  MolRefIndex  =======================

  void MolRefEntry::Init()  {
    name[0]   = char(0);
    classId   = 0;
    symNumber = 1;
    offset    = 0;
    M         = NULL;
  }

  void MolRefEntry::FreeMemory()  {
    if (M)  {
      delete M;
      M = NULL;
    }
  }

  int  MolRefEntry::readRefData ( mmdb::mmcif::PData mmCIF )  {
  int rc;
    if (!M)  M = new MolRefData();
    rc = M->readRefData ( mmCIF );
    if (!rc)  {
      strcpy ( name,M->name );
      classId   = M->classId;
      symNumber = M->symNumber;
    }
    return rc;
  }

  int  MolRefEntry::allocRefData ( ccp4srs::PMonomer monomer )  {
  int rc;
    if (!M)  M = new MolRefData();
    rc = M->allocRefData ( monomer );
    if (!rc)  {
      strcpy ( name,M->name );
      classId   = M->classId;
      symNumber = M->symNumber;
    }
    return rc;
  }

  void MolRefEntry::swap ( RMolRefEntry R )  {
  PMolRefData W;
  char         S[sizeof(name)];
  int          i;

    strcpy ( S,name );
    strcpy ( name,R.name );
    strcpy ( R.name,S );

    i = classId;    classId   = R.classId;    R.classId   = i;
    i = symNumber;  symNumber = R.symNumber;  R.symNumber = i;
    i = offset;     offset    = R.offset;     R.offset    = i;
    W = M;          M         = R.M;          R.M         = W;

  }

  void MolRefEntry::Copy ( RMolRefEntry R )  {
    strcpy ( name,R.name );
    classId   = R.classId;
    symNumber = R.symNumber;
    offset    = R.offset;
    M         = R.M;
  }

  void MolRefEntry::write ( mmdb::io::RFile f )  {
    f.WriteFile ( name,sizeof(name) );
    f.WriteInt  ( &classId   );
    f.WriteInt  ( &symNumber );
    f.WriteLong ( &offset    );
  }

  void MolRefEntry::read ( mmdb::io::RFile f )  {
    FreeMemory();
    f.ReadFile ( name,sizeof(name) );
    f.ReadInt  ( &classId   );
    f.ReadInt  ( &symNumber );
    f.ReadLong ( &offset    );
  }


  #define  STATE_Off     0
  #define  STATE_Index   1
  #define  STATE_Open    2
  #define  STATE_Closed  3

  MolRefIndex::MolRefIndex()  {
    InitMolRefIndex();
  }

  MolRefIndex::~MolRefIndex()  {
    closeData ();
    FreeMemory();
  }

  void MolRefIndex::InitMolRefIndex()  {
    nEntries  = 0;
    R         = NULL;
    MolRefDir = NULL;
    stateId   = STATE_Off;
    lastMRD   = NULL;
  }

  void MolRefIndex::FreeMemory()  {

    FreeIndex();

    if (MolRefDir)  {
      delete[] MolRefDir;
      MolRefDir = NULL;
    }

  }

  void MolRefIndex::FreeIndex()  {
  int i;

    if (R)  {
      for (i=0;i<nEntries;i++)
        R[i].FreeMemory();
      delete[] R;
      R = NULL;
    }
    nEntries = 0;
    lastMRD  = NULL;

  }

  void MolRefIndex::assign ( mmdb::pstr dirPath )  {
    closeData();
    FreeIndex();
    if (dirPath[strlen(dirPath)-1]=='/')
         mmdb::CreateCopy   ( MolRefDir,dirPath );
    else mmdb::CreateCopCat ( MolRefDir,dirPath,"/" );
    stateId = STATE_Off;
  }

  int MolRefIndex::loadIndex()  {
  mmdb::io::File f1;
  mmdb::pstr  S;
  int   rc;

    FreeIndex();

    rc = 0;

    if (MolRefDir)  {
      S = NULL;
      mmdb::CreateCopCat ( S,MolRefDir,molref_index_file );
      f1.assign ( S,false,true );
      if (f1.reset(true))  {
        read ( f1 );
        f1.shut();
        stateId = STATE_Index;
      } else  {
        stateId = STATE_Off;
        rc = -2;
      }
      if (S)  delete[] S;
    } else  {
      stateId = STATE_Off;
      rc = -1;
    }

    return rc;

  }

  int MolRefIndex::openData()  {
  mmdb::pstr S;
  int  rc;

    rc = 0;

    if (MolRefDir)  {
      S = NULL;
      mmdb::CreateCopCat ( S,MolRefDir,molref_data_file );
      f.assign ( S,false,true );
      if (f.reset(true))  {
        stateId = STATE_Open;
        lastMRD = NULL;
      } else  {
        FreeIndex();
        stateId = STATE_Off;
        rc = -2;
      }
      if (S)  delete[] S;
    } else  {
      FreeIndex();
      stateId = STATE_Off;
      rc = -1;
    }

    return rc;

  }

  void MolRefIndex::closeData()  {
    if (stateId==STATE_Open)  {
      f.shut();
      stateId = STATE_Closed;
    }
  }

  bool MolRefIndex::isDataOpen()  {
    return (stateId==STATE_Open);
  }


  PMolRefEntry  MolRefIndex::getMolRefEntry ( mmdb::pstr molName )  {
  mmdb::ResName resName;
  int           l1,l2,l,k;

    if (nEntries<=0)  return NULL;

    alignName ( molName,resName );

    l1 = 0;
    l2 = nEntries-1;
    do  {
      l = (l1+l2)/2;
      k = strcmp ( resName,R[l].name );
      if (k<0)      l2 = l;
      else if (k>0) l1 = l;
               else  break;
    } while (l1<l2-1);

    if (k!=0)  {
      if (l==l1)  l = l2;
            else  l = l1;
      k = strcmp ( resName,R[l].name );
      if (k!=0)
        return NULL;
    }

    return &(R[l]);

  }

  PMolRefData  MolRefIndex::getMolRefData ( mmdb::pstr molName )  {
  PMolRefEntry MRE;

    if (stateId!=STATE_Open)  return NULL;

    if (lastMRD)  {
      if (!strcmp(lastMRD->name,molName))  return lastMRD;
    }

    MRE = getMolRefEntry ( molName );
    if (!MRE)  return NULL;

    if (!MRE->M)  {
      f.seek ( MRE->offset );
      MRE->M = new MolRefData();
      MRE->M->read ( f );
    }

    lastMRD = MRE->M;

    return lastMRD;

  }

  PAtmRefData MolRefIndex::getAtmRefData ( mmdb::PAtom a )  {
  PMolRefData MRD;
    MRD = getMolRefData ( a->GetResName() );
    if (MRD)
      return MRD->getAtmRefData ( a->name );
    return NULL;
  }

  mmdb::realtype MolRefIndex::getAtomRadius ( mmdb::PAtom a )  {
  PAtmRefData ARD;
    ARD = getAtmRefData ( a );
    if (ARD) return ARD->radius;
    return mmdb::getVdWaalsRadius(a->GetElementName());
  }

  void MolRefIndex::allocMREntries ( int & nAlloc )  {
  PMolRefEntry R1;
  int           i;

    if (nEntries>=nAlloc)  {
      nAlloc = nEntries + 1000;
      R1 = new MolRefEntry[nAlloc];
      for (i=0;i<nEntries;i++)  {
        R1[i].Init();
        R1[i].Copy ( R[i] );
      }
      for (i=nEntries;i<nAlloc;i++)
        R1[i].Init();
      if (R)  delete[] R;
      R = R1;
    }

  }


  mmdb::realtype const defaultASP[mmdb::nElementNames] = {
  // == " H",  "HE",
        0.001, 0.001,
  // == "LI", "BE", " B", " C",  " N",  " O", " F", "NE",
        0.05, 0.05, 0.05, 0.016, -0.02, 0.02, -0.02, 0.001,
  // == "NA", "MG", "AL", "SI", " P",  " S", "CL",   "AR",
        0.05, 0.1, 0.05, 0.05, -0.02, 0.04, -0.001, 0.001,
  // == " K", "CA",
        0.05, 0.15,
  // ==        "SC", "TI", " V", "CR", "MN", "FE",
               0.05, 0.05, 0.05, 0.05, 0.10, 0.10,
  // ==        "CO", "NI", "CU", "ZN",
               0.05, 0.05, 0.05, 0.50,
  // ==              "GA", "GE", "AS", "SE", "BR",   "KR",
                     0.05, 0.05, 0.04, 0.04, -0.001, 0.001,
  // ==  "RB", "SR",
         0.05, 0.05,
  // ==        " Y", "ZR", "NB", "MO", "TC", "RU",
               0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
  // ==        "RH", "PD", "AG", "CD",
               0.05, 0.05, 0.05, 0.05,
  // ==              "IN", "SN", "SB", "TE", " I",   "XE",
                     0.05, 0.05, 0.05, 0.05, -0.001, 0.001,
  // ==  "CS", "BA",
         0.05, 0.05,
  // ==        "LA", "CE", "PR", "ND", "PM", "SM", "EU",
               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
  // ==        "GD", "TB", "DY", "HO", "ER", "TM", "YB",
               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
  // ==              "LU", "HF", "TA", " W", "RE", "OS",
                     0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
  // ==              "IR", "PT", "AU", "HG",
                     0.05, 0.05, 0.05, 0.20,
  // ==                    "TL", "PB", "BI", "PO", "AT", "RN",
                           0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
  // ==  "FR", "RA",
         0.05, 0.05,
  // ==        "AC", "TH", "PA", " U", "NP", "PU", "AM",
               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
  // ==        "CM", "BK", "CF", "ES", "FM", "MD", "NO",
               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
  // ==              "LR", "RF", "DB", "SG", "BH", "HS",
                     0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
  // ==              "MT", "UN", "UU", "UB",
                     0.05, 0.05, 0.05, 0.05,
  // ==                          "UQ",       "UH",       "UO",
                                 0.05,       0.05,       0.05,
  // ==  " D", "AN"
         0.001, 0.001
  };


  mmdb::realtype  getDefaultASP ( const mmdb::pstr element )  {
  int type;
    type = mmdb::getElementNo ( element );
    if (type!=mmdb::ELEMENT_UNKNOWN)
         return defaultASP[type-1];
    else return 0.0;
  }

  int  MolRefIndex::fillMRData ( RMolRefEntry      MRE,
                                 ccp4srs::PMonomer monomer,
                                 int               nMR0 )  {
  mmdb::Atom    a1,a2;
  mmdb::realtype radius,asp;
  int      c[nASPs];
  int      rc,i,j,k,nrad,nasp,aspId;

    aspId = 0;  // only to keep compiler happy

    rc = MRE.allocRefData ( monomer );
    if (rc)  return rc;

    for (i=0;i<MRE.M->nAtoms;i++)  {
      strcpy ( a1.name,MRE.M->atom[i].name );
      strcpy ( a1.element,"  "  );
      a1.MakePDBAtomName();
      for (j=0;j<nASPs;j++)
        c[j] = 0;
      radius = 0.0;
      asp    = 0.0;
      nrad   = 0;
      nasp   = 0;
      for (j=0;(j<nMR0) && (nrad>=0);j++)
        for (k=0;(k<R[j].M->nAtoms) && (nrad>=0);k++)
          if (!strcmp(MRE.M->atom[i].name,R[j].M->atom[k].name))  {
            // exact name match - use as hypothesis
            nrad   = -1;  // signal "exact name match"
            radius = R[j].M->atom[k].radius;
            asp    = R[j].M->atom[k].asp;
            aspId  = R[j].M->atom[k].aspId;
          } else  {
            strcpy ( a2.name,R[j].M->atom[k].name );
            strcpy ( a2.element,"  " );
            a2.MakePDBAtomName();
            if (!strcmp(a1.element,a2.element))  {
              radius += R[j].M->atom[k].radius;
              nrad++;
              aspId = R[j].M->atom[k].aspId;
              if (aspId==-2) aspId = 3;
              if (aspId>0)
                c[aspId]++;
              else if (aspId==0)  {
                asp  += R[j].M->atom[k].asp;
                nasp++;
              }
            }
          }
      if (nrad<0)  {
        MRE.M->atom[i].radius = radius;
        MRE.M->atom[i].asp    = asp;
        MRE.M->atom[i].aspId  = aspId;
      } else if (nrad>0)  {
        MRE.M->atom[i].radius = radius/nrad;
        k = -1;
        for (j=0;j<nASPs;j++)
          if (c[j]>k)  {
            k = c[j];
            aspId = j;
          }
        if (k>nasp)  {
          MRE.M->atom[i].asp   = 0.0;
          MRE.M->atom[i].aspId = aspId;
        } else  {
          MRE.M->atom[i].asp   = asp/nasp;
          MRE.M->atom[i].aspId = 0;
        }
      } else  {
        // take defaults from tables
        MRE.M->atom[i].radius = mmdb::getVdWaalsRadius ( a1.element );
        MRE.M->atom[i].asp    = getDefaultASP ( a1.element );
        MRE.M->atom[i].aspId  = 0;
      }
    }

    MRE.M->SortAtoms();

    return 0;

  }


  int MolRefIndex::addSRSEntries ( mmdb::pstr srsDir, int & nAlloc )  {
  ccp4srs::PManager SRS;
  mmdb::io::PFile   structFile;
  ccp4srs::PMonomer monomer;
  mmdb::ResName     resName;
  int               nSRSEntries,nMR0,structNo,i,j,rc;

    SRS = new ccp4srs::Manager();
    if (SRS->loadIndex(srsDir)!=ccp4srs::CCP4SRS_Ok)  {
      delete SRS;
      printf ( " ***** cannot open SRS at\n"
               "       '%s'\n",srsDir );
      return -1;
    }

    structFile = SRS->getStructFile();
    if (!structFile)  {
      delete SRS;
      printf ( " ***** cannot open SRS structure file at\n"
               "       '%s'\n",srsDir );
      return -2;
    }

    nSRSEntries = SRS->n_entries();
    printf ( " ... total %i structures in SRS\n",nSRSEntries );
    nMR0 = nEntries;

    for (structNo=0;structNo<nSRSEntries;structNo++)  {
      monomer = SRS->getMonomer ( structNo,structFile );
      if (monomer)  {
        alignName ( monomer->ID(),resName );
        j = -1;
        for (i=0;(i<nEntries) && (j<0);i++)
          if (!strcmp(resName,R[i].name))  j = i;
        if (j<0)  {
          // new ligand
          allocMREntries ( nAlloc );
          rc = fillMRData ( R[nEntries],monomer,nMR0 );
          if (!rc)  nEntries++;
              else  printf ( "  entry '%s' not added (%i)\n",
                             resName,rc );
        }
        delete monomer;
      }
    }

    structFile->shut();
    delete structFile;
    delete SRS;

    return 0;

  }

  int MolRefIndex::makeMolReference ( mmdb::pstr mmCIFFile,
                                      mmdb::pstr outDir,
                                      mmdb::pstr srsDir )  {
  mmdb::mmcif::PData  mmCIF;
  char                S[1000];
  char                M[1000];
  int                 lcount,rc,nAlloc,i,fNo;
  bool                fSeries;

    closeData();
    FreeIndex();

    fNo = 0;
    fSeries = (mmdb::FirstOccurence(mmCIFFile,'%')!=NULL);

    mmCIF = new mmdb::mmcif::Data();

    do {

      if (fSeries)  {
        sprintf  ( S,mmCIFFile,fNo );
        f.assign ( S,true );
      } else
        f.assign ( mmCIFFile,true );

      if (f.reset(true))  {

        lcount = 0;
        nAlloc = 0;
        rc     = mmdb::mmcif::CIFRC_Ok;
        while (!f.FileEnd())  {
          rc = mmCIF->ReadMMCIFData ( f,S,lcount );
          if (rc==mmdb::mmcif::CIFRC_Ok)  {
            allocMREntries ( nAlloc );
            rc = R[nEntries].readRefData ( mmCIF );
            if (!rc)  nEntries++;
          } else if (rc!=mmdb::mmcif::CIFRC_NoDataLine)
            printf ( " \n *** mmCIF read error:\n"
                     "   %s\n"
                     "   offending line #%i:\n"
                     "   %s\n",
                     mmdb::mmcif::GetCIFMessage(M,rc),lcount,S );
        }

        f.shut();
        rc = 0;
        fNo++;

      } else if (fNo<=0)  {
        printf ( "\n *** error: cannot open file\n"
                 "     '%s'\n",f.FileName() );
        rc = -1;
      } else
        rc = 1;

      if (!fSeries)  rc = 2;

    } while (!rc);

    delete mmCIF;

    if (rc<0)  return -1;

    printf ( " --- total %i input entries in the reference\n",nEntries );

    if (srsDir)
      addSRSEntries ( srsDir,nAlloc );

    SortIndex    ();
    CheckRepeats ();
    CheckASPRange();

    printf ( " --- total %i entries in the reference\n",nEntries );

    strcpy ( S,outDir );
    if (S[strlen(S)-1]!='/')  strcat ( S,"/" );
    strcat ( S,molref_data_file );

    f.assign ( S,false,true );
    if (f.rewrite())  {
      for (i=0;i<nEntries;i++)  {
        R[i].offset = f.Position();
        R[i].M->write ( f );
      }
      f.shut();
      strcpy ( S,outDir );
      if (S[strlen(S)-1]!='/')  strcat ( S,"/" );
      strcat ( S,molref_index_file );
      f.assign ( S,false,true );
      if (f.rewrite())  {
        write ( f );
        f.shut();
        rc = 0;
      } else  {
        printf ( "\n *** cannot write file\n"
                 "    '%s'\n",S );
        rc = -3;
      }
    } else  {
      printf ( "\n *** cannot write file\n"
               "    '%s'\n",S );
      rc = -2;
    }

    for (i=nEntries;i<nAlloc;i++)
      R[i].FreeMemory();

    return rc;

  }

  void MolRefIndex::SortIndex()  {
  int i,j;

    for (i=0;i<nEntries;i++)
      for (j=i+1;j<nEntries;j++)
        if (strcmp(R[i].name,R[j].name)>0)
          R[i].swap ( R[j] );

  }

  void MolRefIndex::CheckRepeats()  {
  int i,j;

    for (i=0;i<nEntries;i++)  {
      for (j=i+1;j<nEntries;j++)
        if (!strcmp(R[i].name,R[j].name))
          printf ( " *** compound %s is entered more than once.\n",
                   R[i].name );
      if (R[i].M)  R[i].M->CheckRepeats();
    }

  }

  void MolRefIndex::CheckASPRange()  {
  int i;
    for (i=0;i<nEntries;i++)
      if (R[i].M)  R[i].M->CheckASPRange();
  }

  void MolRefIndex::listEntries()  {
  int rc,i;

    if (!MolRefDir)  {
      printf ( " *** reference not initialized\n" );
      return;
    }

    rc = loadIndex();
    if (rc)  {
      printf ( " *** cannot read the reference index (%i)\n",rc );
      return;
    }

    printf ( "\n  Total %i entries in the reference.\n\n",nEntries );
    for (i=0;i<nEntries;i++)
      printf ( " %5i. '%s'  %2i %5i\n",i+1,R[i].name,R[i].classId,
                                       R[i].symNumber );

  }

  void MolRefIndex::write ( mmdb::io::RFile f1 )  {
  int i;

    f1.WriteInt ( &nEntries );
    for (i=0;i<nEntries;i++)
      R[i].write ( f1 );

  }

  void MolRefIndex::read ( mmdb::io::RFile f1 )  {
  int i;

    FreeIndex();

    f1.ReadInt ( &nEntries );
    if (nEntries>0)  {
      R = new MolRefEntry[nEntries];
      for (i=0;i<nEntries;i++)  {
        R[i].Init ();
        R[i].read ( f1 );
      }
    }

  }

  void MolRefIndex::rewriteMolRef()  {
  mmdb::io::File f1;
  mmdb::pstr     S;
  int            i;

    loadIndex();
    if (!isDataOpen())
      openData();

    for (i=0;i<nEntries;i++)  {
      f.seek ( R[i].offset );
      R[i].M = new MolRefData();
      R[i].M->read ( f );
    }

  mmdb::__modify4();

    S = NULL;
    mmdb::CreateCopCat ( S,MolRefDir,molref_data_file,".new" );
    printf ( " S='%s'\n",S );
    f1.assign ( S,false,true );
    f1.rewrite();
    for (i=0;i<nEntries;i++)  {
      R[i].offset = f1.Position();
      R[i].M->write ( f1 );
    }
    f1.shut ();

    mmdb::CreateCopCat ( S,MolRefDir,molref_index_file,".new" );
    printf ( " S='%s'\n",S );
    f1.assign ( S,false,true );
    f1.rewrite();
    write ( f1 );
    f1.shut();
    delete[] S;

  }

}  // namespace pisa
