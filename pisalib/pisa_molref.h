// $Id: pisa_molref.h $
// =================================================================
//
//    19.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_molref <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::AtmRefData
//       ~~~~~~~~~  pisa::MolRefData
//                  pisa::MolRefIndex
//
//  (C) E. Krissinel, 2006-2013
//
// =================================================================
//

#ifndef __PISA_MolRef__
#define __PISA_MolRef__

#include "ccp4srs/ccp4srs_monomer.h"

namespace pisa  {

  // =========================  AtmRefData  =========================

  DefineStructure(AtmRefData);

  struct AtmRefData  {

    char           name[5],old_name[5];
    mmdb::realtype radius,refASA,asp;
    int            aspId;

    int  readRefData  ( mmdb::mmcif::PLoop mmCIFLoop, int aNo );
    void allocRefData ( mmdb::pstr atname, mmdb::pstr old_atname );
    void swap         ( RAtmRefData a );

    void write ( mmdb::io::RFile f );
    void read  ( mmdb::io::RFile f );

  };


  // =========================  MolRefData  =========================

  enum  MOLREF_CLASS  {
    MRCLASS_Unknown       = 0,
    MRCLASS_Aminoacid     = 1,
    MRCLASS_Nucleotide    = 2,
    MRCLASS_ModAminoacid  = 3,
    MRCLASS_ModNucleotide = 4,
    MRCLASS_Ligand        = 5,
    MRCLASS_Solvent       = 6
  };

  DefineClass(MolRefData);

  class MolRefData  {

    public :
      char     name[4];
      char  refName[4];
      int      classId;
      int    symNumber;
      int       nAtoms;
      PAtmRefData atom;

      MolRefData ();
      ~MolRefData();

      int  readRefData  ( mmdb::mmcif::PData mmCIF   );
      int  allocRefData ( ccp4srs::PMonomer  monomer );
      void SortAtoms    ();

      PAtmRefData getAtmRefData ( mmdb::pstr atmName );
      void CheckRepeats ();
      void CheckASPRange();

      void write ( mmdb::io::RFile f );
      void read  ( mmdb::io::RFile f );

    protected :
      void InitMolRefData();
      void FreeMemory    ();

  };


  // =========================  MolRefIndex  =======================

  DefineStructure(MolRefEntry);

  struct MolRefEntry  {

    char   name[4];
    int    classId;
    int  symNumber;
    long    offset;
    PMolRefData  M;

    void Init      ();
    void FreeMemory();

    int  readRefData  ( mmdb::mmcif::PData mmCIF   );
    int  allocRefData ( ccp4srs::PMonomer  monomer );

    void swap  ( RMolRefEntry R );
    void Copy  ( RMolRefEntry R );

    void write ( mmdb::io::RFile f );
    void read  ( mmdb::io::RFile f );

  };

  DefineClass(MolRefIndex);

  extern mmdb::cpstr molref_index_file;
  extern mmdb::cpstr molref_data_file;

  class MolRefIndex  {

    public :
      int   nEntries;
      PMolRefEntry R;

      MolRefIndex ();
      ~MolRefIndex();

      void  assign   ( mmdb::pstr dirPath );
      int   loadIndex();
      int   openData ();
      void  closeData();

      bool  isDataOpen();

      PMolRefEntry getMolRefEntry ( mmdb::pstr molName );
      PMolRefData  getMolRefData  ( mmdb::pstr molName );

      PAtmRefData    getAtmRefData ( mmdb::PAtom a );
      mmdb::realtype getAtomRadius ( mmdb::PAtom a );

      int  makeMolReference ( mmdb::pstr mmCIFFile, mmdb::pstr outDir,
                              mmdb::pstr SBaseDir );
      void listEntries();

      void rewriteMolRef();  // service function

    protected :
      mmdb::io::File f;
      PMolRefData    lastMRD;
      mmdb::pstr     MolRefDir;
      int            stateId;

      void InitMolRefIndex ();
      void FreeMemory      ();
      void FreeIndex       ();

      void allocMREntries  ( int & nAlloc );
      int  fillMRData      ( RMolRefEntry   MRE,
                             ccp4srs::PMonomer monomer,
                             int nMR0 );
      int  addSRSEntries   ( mmdb::pstr srsDir, int & nAlloc );
      void SortIndex       ();
      void CheckRepeats    ();
      void CheckASPRange   ();

      void write ( mmdb::io::RFile f1 );
      void read  ( mmdb::io::RFile f1 );

  };

}  // namespace pisa

#endif
