// $Id: pisa_engine.h $
// =================================================================
//
//    29.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  asm_engine <interface>
//       ~~~~~~~~~
//  **** Project :  Protein interfaces
//       ~~~~~~~~~
//  **** Classes :  pisa::Assembler
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2004-2013
//
// =================================================================
//

#ifndef __PISA_Engine__
#define __PISA_Engine__

#include "time.h"

#include "pisa_assemblies.h"

namespace pisa  {

  // =========================  Assembler  ==========================

  DefineClass(Assembler);

  typedef mmdb::vect3 * pvect3;
  typedef pvect3      * ppvect3;
  typedef ppvect3     * pppvect3;

  class Assembler  {

    public :

      ///  selection rules for fit matrix
      enum SELECTION_RULE  {
        SELRULE_DeltaG = 0,
        SELRULE_BSA    = 1
      };

      Assembler ();
      ~Assembler();

      void Reset            ();

      void SetDissThreshold ( mmdb::realtype threshold );
      void SetLigandKey     ( int            ligandKey );
      void SetFullSearch    ( bool           fSearch   );

      void Set_as_is_Key ( int as_is_Key );

      ASSMB_RC MakePAGraph    ( mmdb::PManager MMDBManager,
                                PDomains       domains,
                                PInterfaces    PI );

      ASSMB_RC CalcAssemblies ();
      ASSMB_RC CalcAssemblies ( mmdb::PManager MMDBManager,
                                PDomains       domains,
                                PInterfaces    PI );
      void GetAssemblies      ( RPAssemblies   assemblies );

      ASSMB_RC AnalyseComplex ( mmdb::PManager MMDBManager,
                                PDomains       domains,
                                PInterfaces    PI );
      void GetComplex         ( RPAssembly     assembly );

      int  GetRC            ();

    protected :
      SymNumber      SNum;        //!< symmetry number calculator
      PDomains       D;           //!< pointers retained after makePAGraph(),
      mmdb::PManager MMDB;        //!<       never allocated nor disposed
      PPInterface    interface;   //!< allocated and resorted
      PPMonomer      M;           //!< molecules
      pppvect3       ucv;         //!< unit cell translations
      PPMultimer     U;           //!< vector of multimers
      PPIntfRecStack recStack;    //!< recursion stack
      PPMultimerSet  multSet;     //!< sets of resulting multimers
      PMultimer      Complex;     //!< multimer for single-complex analysis
      mmdb::mat44    rom;         //!< orthogonalisation matrix
      mmdb::omatrix  parMon;      //!< matrix of parallel monomers
      mmdb::ivector  eType;       //!< nterface enumeration types
      int            ligKey;      //!< ligand proccessing key
      int            asisKey;     //!< assembly calculation key 
      ASSMB_RC       retcode;     //!< return code
      int            nSymOps;     //!< number of symmetry operations
      int            nInterfaces; //!< calculated number of interfaces
      int            nInterfaces0;//!< original number of interfaces
      int            nMonomers;   //!< number of monomers in unit cell
      int            nCellLayers; //!< number of unit cell layers
      int            nETypes;     //!< number of enumeration types
      int            nMultimers;  //!< current number of multimers
      int            nMultSets;   //!< number of resulting sets
      bool           Crystal;     //!< True if crystal is given
      bool           fullSearch;  //!< True if full search is requested

      ///   enumeration parameters to control the enumeration length
      ///   for large structures
      mmdb::realtype dissThreshold;  //!< dissociation threshold
      int            maxStableSize;
      int            fixLigandSize;  //!< max. size of ligands to fix compulsory
      int            fixLigandNum;   //!< maximal number of non-fixed ligands
      int            noImproveCnt;
      int            noImproveLimit; //!< enumeration limit for large structures
      int            runStatus;
      int            timeLimit;      //!< CPU time limit in seconds

      PPPPBrick      MB;
      mmdb::realtype brick_size;
      mmdb::realtype x1_brick,x2_brick, y1_brick,y2_brick, z1_brick,z2_brick;
      int            nx_bricks,ny_bricks,nz_bricks;

      void InitAssembler      ();
      void DeletePAGraph      ();
      void DeleteMonomers     ( PPMonomer  & Mon,  int & nM, int & nMA );
      void DeleteMultimers    ( PPMultimer & Mult, int & nM, int & nMA );
      void DeleteComplex      ();
      void DeleteRecStack     ();
      void DeleteMultimerSets ( PPMultimerSet & mSet, int & nM, int & nMA );
      void RemoveBricks       ();

      ASSMB_RC MakeInterfaces ( PInterfaces PI, PPDomain Domain );
      bool traceMonomerChain  ( PMonomer M1, int engType,
                                int i1, int j1, int k1,
                                PMonRef L, int cnt );
      void CalcEnumerationTypes ();
      ASSMB_RC MakeCrystalMonomers();
      ASSMB_RC MakeComplexMonomers();
      int  FixInterface         ( int iNo, mmdb::ovector fixed, PMonRef G,
                                  mmdb::ovector mmol );
      void FixLigands           ();
      void MakeBricks           ();
      void getBricks            ( PPBrick brick, mmdb::mat44 & uct,
                                  PDomain dom );
      void CalcCrystalInterfaces();
      void CalcComplexInterfaces();
      void CalcParallelMonomers ();
      void SeedMultimers        ();
      int  makeInterfaceFlags   ( mmdb::ivector iflag );
      int  MakeRecStack         ();

      int  NoteMultimer         ( int r );
      bool EngageInterface      ( int p, int r  );
      void RollBack             ( int r );
      void EnumerateInterfaces  ( int r, int i0 );
      void CalcEqIds            ();

    private :
      clock_t startCPUclock;
      int     nMonAlloc,nMultAlloc,nRSAlloc,nMultSetAlloc;

  };


  // =================================================================

  extern void SetMaxNofETypes ( int maxNETypes );

}  // namespace pisa

#endif
