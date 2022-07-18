// $Id: pisa_domains.h $
// =================================================================
//
//    01.11.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_domains <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Domains
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2007-2013
//
// =================================================================
//

#ifndef __PISA_Domains__
#define __PISA_Domains__

#include "pisa_domain.h"

namespace pisa  {

  //  ========================  Domains  =============================

  DefineStructure(DAlign);

  struct DAlign  {
    mmdb::realtype Qscore,rmsd,seqId;
    mmdb::mat44    TMatrix;
    int            Nalign;
    bool           equal;
    void  Init     ();
    void  calcEqual();
    void  Copy  ( RDAlign DA, bool direct );
    void  write ( mmdb::io::RFile f );
    void  read  ( mmdb::io::RFile f );
  };

  DefineClass(Brick)
  typedef PPBrick  * PPPBrick;
  typedef PPPBrick * PPPPBrick;

  class Brick  {
    public :
      mmdb::ivector n;
      int           nObjects;
      Brick ();
      ~Brick();
      void AddObject ( int objectNo );
    protected :
      int nAlloc;
  };

  DefineStreamFunctions(Domains);

  class Domains : public mmdb::io::Stream  {

    public :
      PPDomain       domain; //!< array of monomeric units ("domains")
      PPDAlign       DA;     //!< symmetric, only D[i][j>i] are allocated
      ssm::PPGraph   G;      //!< SSM structral graphs
      mmdb::realtype Dmax;   //!< maximal domain radius for bricking
      int            nDomains;     //!< total number of domains
      int            nDTypes;      //!< total number of domain types
      int            nNCSOps;      //!< number of NCS operations in ASU
      int            nNCSParents;  //!< number of NCS paren domains

      Domains ();
      Domains ( mmdb::io::RPStream Object );
      ~Domains();

      void  FreeMemory();

      void  calcDimensions   ( mmdb::PManager MMDB );
      void  MakeDomains      ( mmdb::PManager MMDB,
                               PMolRefIndex   MolRef,
                               mmdb::cpstr    AgentsRef,
                               mmdb::pstr   & warnUnk,
                               mmdb::pstr   & warnExcl );

      PROSURF_RC CalcSurfaces ( mmdb::PManager MMDB,
                                PProSurf       ProSurf,
                                PMolRefIndex   MolRef );

      void  AllocateAlignments();
      void  AlignSimilar     ( mmdb::PManager MMDB );
      void  makeSSGraphs     ( mmdb::PManager MMDB,
                               mmdb::ivector & selHnd );
      void  completeSSGraphs ( mmdb::PManager MMDB,
                               mmdb::ivector & selHnd );
      void  AlignAll         ( mmdb::PManager MMDB,
                               int startNo, int endNo );
      void  CalcTypes        ();
      void  SortDomains      ();

      int   getNnotAligned   ();
      void  writeAlignments  ( mmdb::io::RFile f );

      bool  isAligned     ( int domain1, int domain2 );
      bool  isEquivalent  ( int domain1, int domain2 );
      mmdb::realtype getQscore     ( int domain1, int domain2 );
      mmdb::realtype getRMSD       ( int domain1, int domain2 );
      mmdb::realtype getSeqId      ( int domain1, int domain2 );
      int      getNalign     ( int domain1, int domain2 );
      void     getTMatrix    ( mmdb::mat44 & TMatrix, int domain1,
                               int domain2 );
      mmdb::pstr getDomainRange      ( mmdb::pstr S, int domainNo,
                                       int fieldLen );
      int        getDomainClass      ( int domainNo );
      mmdb::pstr getChainID          ( mmdb::pstr S, int domainNo );
      mmdb::pstr getDomainRange_html ( mmdb::pstr S, int domainNo );
      void       getNofDomains ( int & nAADomains, int & nNADomains,
                                 int & nLigands );
      PDomain  getDomain     ( int domainNo );  //!< 0..nDomains-1
      int      getNCSParent  ( int domainNo );
      int      getNAtoms     ( int domainNo );
      void     getLigandList ( mmdb::PResName & ligName,
                               int & nLigNames );

      /// Returns the occurence number of specified monomer type in the
      /// crystal.
      int getMonTypeOccurence ( int monType );

      bool  getAssembleStatus ( mmdb::pstr ligandName );
      bool  setAssembleStatus ( mmdb::pstr ligandName, bool assemble );
      bool  assembleAll       ();
      void  getExclLigandList ( mmdb::PResName & ligName,
                                int & nLigNames );

      /// ligList must contain comma-separated list of ligand names
      /// to be excluded from oligomeric state analysis. The list
      /// must start and finish with commas, and contain no spaces.
      ///   Example:  ",SO4,Zn,"
      void  excludeLigands     ( mmdb::pstr ligList );
      void  setLigandsAssemble ( bool On            );

      void  RemoveBricks  ();
      void  MakeBricks    ();
      void  getBricks     ( PPBrick brick, PDomain D,
                            mmdb::mat44 & TMatrix );

      void write ( mmdb::io::RFile f );
      void read  ( mmdb::io::RFile f );

    protected :
      PPPPBrick      DB;
      mmdb::realtype brick_size;
      mmdb::realtype x1_brick,x2_brick;
      mmdb::realtype y1_brick,y2_brick;
      mmdb::realtype z1_brick,z2_brick;
      int            nx_bricks,ny_bricks,nz_bricks;

      void InitDomains();
      void FreeAligns ();
      void FreeGraphs ();
      void allocateDomains ( int & nAlloc );
      void excludeLigands  ( int symFactor, int maxNofLigands,
                             mmdb::pstr & warnExcl );
      void randomizeReferenceFrames();
      void correctReferenceFrames  ();

  };

  extern bool Mat4Compare ( mmdb::mat44 & t1, mmdb::mat44 & t2   );
  extern bool isMat4Rot   ( mmdb::mat44 & t , mmdb::realtype eps );
  extern void SetMaxNofLigands ( int maxNLigandsUC, int maxNLigandsASU );
  extern bool CompareReferenceFrames ( RTFrame F1, RTFrame F2,
                                       mmdb::realtype eps );
  extern bool ParallelFrames         ( RTFrame F1, RTFrame F2,
                                       mmdb::realtype eps );


}  // namespace pisa

#endif
