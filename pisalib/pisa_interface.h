// $Id: pisa_interface.h $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_interface <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Interface
//       ~~~~~~~~~  pisa::Interfaces
//
//  (C) E. Krissinel, 2007-2014
//
// =================================================================
//

#ifndef __PISA_Interface__
#define __PISA_Interface__

#include "mmdb2/mmdb_xml_.h"
#include "ccp4srs/ccp4srs_manager.h"

#include "pisa_domains.h"
#include "pisa_rcsbdat.h"


namespace pisa  {

  // ========================  Interface  ===========================

  DefineStructure(IBond);

  struct IBond  {
    int     serNum1,serNum2; //!< serial numbers of bonded atoms
    int     res1,res2;       //!< residue positions (0,1,...)
    mmdb::realtype dist;
    void SetAtomPair ( ccp4srs::RAtomPair AP );
    void SetContact  ( mmdb::RContact     C  );
    void Copy        ( RIBond             IB );
    void write       ( mmdb::io::RFile    f  );
    void read        ( mmdb::io::RFile    f  );
  };

  DefineStructure(SREffect);

  struct SREffect  {
    mmdb::realtype effect;  //!< single residue effect, kcal/mol
    int            resNo;   //!< residue number (position in the chain)
    void SetEffect ( mmdb::realtype eff, int resN );
    void Copy      ( RSREffect sre     );
    void write     ( mmdb::io::RFile f );
    void read      ( mmdb::io::RFile f );
  };

  extern void MakeDFrame ( RDFrame dframe, PDomain D1, PDomain D2,
                           mmdb::mat44 & TMatrix );
  extern void MakeDFrame ( RDFrame dframe,
                           PDomain D1, mmdb::mat44 & TM1,
                           PDomain D2, mmdb::mat44 & TM2 );

  DefineClass(Interface);
  DefineStreamFunctions(Interface);

  class Interface : public mmdb::io::Stream  {

    public :
      int            id;                //!< interface serial number
      int            type;              //!< equivalence type 1,2,...
      int            domain1,domain2;   //!< 0,1,...
      DOMAIN_CLASS   dclass1,dclass2;   //!< DCLASS_XX
      int            symOpNo;           //!< 0..nSymOps-1
      int            rcsb_symop;        //!< rcsb symop serial number
      int            cell_i,cell_j,cell_k;
      int            nOcc;              //!< occurence number in ASU
      mmdb::mat44    TMatrix;           //!< to be applied to domain2
      DFrame         dframe;            //!< interframe distances
      int            selHndInt1 ,selHndInt2;
      int            nIntAtoms1 ,nIntAtoms2;
      int            nRes1      ,nRes2; //!< number of residues in domains
      int            nIntRes1   ,nIntRes2;
      mmdb::realtype intArea1   ,intArea2   ,intArea;
      mmdb::realtype intDeltaG1 ,intDeltaG2 ,intDeltaG;
      mmdb::realtype aveDeltaG1 ,aveDeltaG2 ,aveDeltaG;
      mmdb::realtype probDeltaG1,probDeltaG2,probDeltaG;
      mmdb::realtype PValue1    ,PValue2    ,PValue;
      mmdb::realtype DeltaSAS1[nASPs];  //!< partial atom surface
      mmdb::realtype DeltaSAS2[nASPs];  ///             area changes
      mmdb::realtype stabEn;            //!< stabilization energy
      mmdb::realtype css;               //!< complexation significance score
      mmdb::pstr     symOp;
      PIBond         HBond,SBridge,DSBond,CovBond,OthBond;
      int            nHBonds,nSBridges,nDSBonds,nCovBonds,nOthBonds;
      mmdb::rvector  bsa1   ,bsa2;      //!< buried surface area per residue
      mmdb::rvector  solvEn1,solvEn2;   //!< solvation energy changes per residue
      PSREffect      sre_stab1,sre_stab2,sre_destab1,sre_destab2;
      int            nStab1,nStab2,nDestab1,nDestab2;
      bool           overlap;           //!< flag "overlap suspicion"
      bool           casual;            //!< flag "casual interface"
      bool           Xrel;              //!< flag "crystallographically related"
      bool           fixedLigand;       //!< flag "fixed ligand interface"

      Interface ();
      Interface ( mmdb::io::RPStream Object );
      ~Interface();

      void   set_asisKey_xml_int ( AS_IS_KEY as_is_Proc_int );
      void   Reset      ();
      void   SetContact ( int   dom1, int  dom2, PDomains D,
                          int nSymOp, mmdb::mat44 & TM,
                          int  icell, int jcell, int kcell );

      bool  checkTMatrix ( mmdb::mat44 & TM );
      bool  checkDFrame  ( RDFrame df, bool direct );
      mmdb::realtype getStabEn();

      PROSURF_RC calcInterface ( mmdb::PManager    MMDB,
                                 PDomain           D1,
                                 PDomain           D2,
                                 PProSurf          proSurf,
                                 PMolRefIndex      molRef,
                                 ccp4srs::PManager SRS );

      bool  isSimilar (
                   PDomain A, mmdb::mat44 & Ta, mmdb::realtype rmsdA,
                   PDomain B, mmdb::mat44 & Tb, mmdb::realtype rmsdB,
                   mmdb::mat44 & Tm );

      int  isEquivalent ( PInterface Interface, PDomains Domains );

      mmdb::xml::PXMLObject getXML ( mmdb::PManager MMDB,
                                     PPDomain Domain, int as_is_param );

      void  Copy  ( PInterface interface );
      void  write ( mmdb::io::RFile f );
      void  read  ( mmdb::io::RFile f );

    protected :

      AS_IS_KEY  asisKey ;
      void  InitInterface      ();
      void  FreeMemory         ();
      void  deleteHSDBonds     ();
      void  deleteCovBonds     ();
      void  deleteOthBonds     ();
      void  deleteSREffects    ();

      void  checkOverlap       ( mmdb::PManager MMDB,
                                 mmdb::PPAtom   atom1, int nat1,
                                 mmdb::PPAtom   atom2, int nat2 );
    //GLD: added function to calculate other contacts that are not h-bonds,s-bridges, cov-bonds and disulfide bonds
      void  CalcOtherContacts ( mmdb::PManager MMDB,
                                 mmdb::PPAtom atom1, int nat1,
                                 mmdb::PPAtom atom2, int nat2);
    
      void  calcChemBonds      ( mmdb::PManager MMDB,
                                 ccp4srs::PManager SRS);

      void  calcResBSA         ( mmdb::PManager  MMDB,
                                 int             selHndA,
                                 mmdb::rvector & bsa,
                                 mmdb::rvector & solvEn,
                                 int           & nRes,
                                 mmdb::rvector   intSAS,
                                 mmdb::rvector   atomSE );

      void  calcSREffs         ( mmdb::rvector   solvEn,
                                 mmdb::rvector   asa,
                                 mmdb::rvector   bsa,
                                 int             nRes,
                                 int             domNo,
                                 RPSREffect      sre_stab,
                                 int           & nStab,
                                 RPSREffect      sre_destab,
                                 int           & nDestab );

      mmdb::xml::PXMLObject getStructXML (
                                 mmdb::PManager MMDB,
                                 PPDomain       domain,
                                 int            structNo,
                                 int            dclass,
                                 int            symop,
                                 int            icell,
                                 int            jcell,
                                 int            kcell,
				 int  	        asis_param,
                                 mmdb::mat44  & t,
                                 int            nIntAtoms,
                                 int            nIntRes,
                                 mmdb::realtype intArea,
                                 mmdb::realtype intDeltaG,
                                 mmdb::realtype PValue,
                                 mmdb::cpstr    symOpn,
                                 mmdb::rvector  bsa,
                                 mmdb::rvector  solvEn );

      mmdb::xml::PXMLObject makeBondTXML ( PIBond Bond, int nBonds,
                                 mmdb::cpstr Tag, mmdb::PManager MMDB );

      void  writeMaxFrameDiff  ( mmdb::mat44 & df );

  };


  // ========================  Interfaces  ==========================

  DefineClass(Interfaces);
  DefineStreamFunctions(Interfaces);

  class Interfaces : public mmdb::QuickSort  {

    public :
      bool rcsb_symops;  // True if rcsb symops have been assigned

      Interfaces ();
      Interfaces ( mmdb::io::RPStream Object );
      ~Interfaces();

      void  Reset();

      void  AddContact ( PDomains D, int  domain1, int domain2,
                         int nSymOpNo, mmdb::mat44 & TMatrix,
                         int cell_i, int  cell_j, int cell_k );

      void  AddInterface ( PInterface Interface );
      RESULT_CODE CalcContacts ( mmdb::PManager MMDB, PDomains D,
                                 int nCellLayers );

      void  DeleteDummyInterfaces();
      void  CalcInterfaceTypes   ( PDomains  domains  );
      void  Sort                 ( ISORT_KEY sortmode );

      int   AnalyseInterfaces    ( PDomains Domains );

      inline int getNofInterfaces()  { return nInterfaces; }
      inline int getNofITypes    ()  { return nITypes;     }

      PInterface  getInterface ( int intfNo );  // 0,1,...n-1
      inline PPInterface getInterfaces() { return PI; }

      int   getSimilarInterface (
        PDomain A, mmdb::mat44 & Ta, int domain1, mmdb::realtype rmsdA,
        PDomain B, mmdb::mat44 & Tb, int domain2, mmdb::realtype rmsdB,
        mmdb::mat44 & Tm, PDomains Domains );

      void  assignRCSBSymOps ( PRCSBData rcsbData );

      mmdb::xml::PXMLObject getXML ( mmdb::PManager MMDB,
                                     PPDomain Domain, int as_is_param );

      void  write ( mmdb::io::RFile f );
      void  read  ( mmdb::io::RFile f );

    protected :
      PPInterface      PI;
      INTERFACE_STATUS intStatus; // interface status
      int              nInterfaces;
      int              nITypes;   // number of equivalent interface types
      ISORT_KEY        smode;     // sorting mode

      void  InitInterfaces ();
      void  FreeMemory     ();
      void  checkAllocation();

      // sorting virtuals
      int  Compare ( int i, int j );
      void Swap    ( int i, int j );

    private :
      int nAlloc;

  };

  extern void  SetASPFitSE ( bool fitMode );

  extern void  writeBondFactors   ( mmdb::io::RFile f );
  extern void  readBondFactors    ( mmdb::io::RFile f );
  extern void  SetBondFactors     ( mmdb::realtype HBEn,
                                    mmdb::realtype SBEn,
                                    mmdb::realtype DSEn,
                                    mmdb::realtype epsf );
  extern void  GetBondFactors     ( mmdb::realtype & HBEn,
                                    mmdb::realtype & SBEn,
                                    mmdb::realtype & DSEn,
                                    mmdb::realtype & epsf );
  extern void  makeBondFactorsCIF ( mmdb::mmcif::PStruct mmCIFStruct );
  extern void  getBondFactorsCIF  ( mmdb::mmcif::PStruct mmCIFStruct );
  extern mmdb::PAtom GetAtom ( mmdb::PManager MMDB, int  serNum );
  extern void  addXML0       ( mmdb::xml::PXMLObject xml,
                               mmdb::cpstr tag,
                               mmdb::realtype V );
  extern bool isFixed ( int domainId, PPInterface I, int nInterfaces );


}  // namespace pisa

#endif

