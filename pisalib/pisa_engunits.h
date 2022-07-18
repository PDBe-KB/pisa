// $Id: pisa_asmunits.h $
// =================================================================
//
//    19.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  AsmUnits <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::Monomer
//       ~~~~~~~~~  pisa::Multimer
//                  pisa::IntfRecStack
//                  pisa::MultimerSet
//
//  (C) E. Krissinel 2007-2013
//
// =================================================================
//

#ifndef __PISA_AsmUnits__
#define __PISA_AsmUnits__

#include "pisa_domain.h"
#include "pisa_interface.h"
#include "pisa_symnumber.h"

namespace pisa  {

  //  ==================  Assembly Unit Classes  =====================

  enum INTERFACE_FLAG  {
    INTF_Undefined = -5,
    INTF_Off       = -4,
    INTF_Void      = -3,
    INTF_Open      = -2,
    INTF_Fixed     = -1
  };

  extern void writeEntropyFactors ( mmdb::io::RFile f );
  extern void readEntropyFactors  ( mmdb::io::RFile f );

  extern void makeEntropyFactorsCIF ( mmdb::mmcif::PStruct mmCIFStruct );
  extern void getEntropyFactorsCIF  ( mmdb::mmcif::PStruct mmCIFStruct );

  extern void SetEntropyFactors ( mmdb::realtype ent_const,
                                  mmdb::realtype ent_rbody,
                                  mmdb::realtype ent_surf );

  extern void GetEntropyFactors ( mmdb::realtype & ent_const,
                                  mmdb::realtype & ent_rbody,
                                  mmdb::realtype & ent_surf );

  extern bool repeated_assignment;
  extern bool multimer_overflow;


  DefineClass(Monomer);
  DefineClass(Multimer);


  //  =========================  MonRef  ============================

  DefineStructure(MonRef);

  struct MonRef  {
    PMonomer M;      //!< interface-linked monomer (not allocated)
    int      i,j,k;  //!< relative cell shift
    bool     c;      //!< used for specifying internal interfaces
    void  Init();
    bool  isInternal ( RMonRef L );
    void  GetTMatrix ( mmdb::mat44 & T, mmdb::mat44 & rom );
  };

  //  =========================  Monomer  ===========================

  #define _max_n_int  10

  class Monomer  {

    public :
      int         id;      //!< monomer serial number in ASU
      int         ncsParent; //!< monomer ncs parent = domain->ncsParent
      int         ix;      //!< monomer serial number in unit cell
      int         mx;      //!< monomer serial number in multimer
      int         id2;     //!< auxiliary id for temporary use
      int         symOpNo; //!< symmetry operation (ASU) number
      int         dclass;  //!< protein, dna/rna or ligand
      mmdb::mat44 uct;     //!< transformation to ASU in cell 333
      PMonRef     L[_max_n_int]; //!< interface type link vectors:
                           //!<   macromolecular monomeric units may make
                           //!<   the same interface twice, one- and
                           //!<   two-atom ligands may make more identical
                           //!<   interfaces
      PMultimer U;         //!< reference to multimer

      Monomer ( int monId, int monNCSParent, int monDClass,
                mmdb::mat44 & TMatrix, int nInterfaces );
      ~Monomer();

      void   AddInterface ( PMonomer Mon, int interfaceNo,
                            int i, int j, int k, int nInterfaces );

      bool  isParallel ( PMonomer Mon );

      void  getMassCenter ( PDomains         D,
                            mmdb::realtype & x,
                            mmdb::realtype & y,
                            mmdb::realtype & z );


  };


  //  =========================  Multimer  ==========================

  DefineStructure(DIRecursion);

  struct DIRecursion  {
    PDomains      domains;
    PPDomain      D;
    PPInterface   interface;
    PSymNumber    SNum;
    mmdb::vect3 * mc;
    mmdb::pmat44  TMatrix;
    mmdb::ivector eType;
    int           nInterfaces;
    int           nMMInt;  // number of m'molecular interfaces
    int           nit;     // number of engaged significant interfaces

    void Init ( PDomains        Doms,
                PPInterface     Ints,
                PSymNumber      symNum,
                mmdb::ivector   intType,
                int             nInts,
                PMultimer       U,
                int             nETypes,
                mmdb::mat44   & rom,
                mmdb::ivector & itp );

    void FreeMemory();

  };


  class Multimer  {

    public :
      PMonRef        R;           //!< list of monomers in assembly
      int            mSize;       //!< multimer size (number of monomers)
      mmdb::realtype dissEn;      //!< maximal free energy of dissociation
      mmdb::realtype entropy;     //!< entropy change at dissociation
      mmdb::realtype dissIntArea; //!< dissociation interface area
      mmdb::ivector  dissMult;    //!< dissociation subunits ids
      mmdb::realtype dissEn0;     //!< ground-level free energy of dissociation
      mmdb::realtype entropy0;    //!< ground-level entropy change at dissociation
      mmdb::realtype asa;         //!< accessible surface area
      mmdb::realtype bsa;         //!< buried surface area
      mmdb::realtype seGain;      //!< solvation energy gain
      int            nUC;         //!< number of these assemblies in unit cells
      int            type;        //!< assembly type id in assembly set
      bool           stable;      //!< True if multimer is stable

      Multimer ( int maxSize );
      ~Multimer();

      void  FreeMemory();
      void  MakeMultimer      ( PMonomer Monomer );
      void  MakeMultimer      ( PPMonomer Monomer, int nMonomers );
      void  MakeBoundMultimer ( PPMonomer Monomer, int nMonomers );
      void  SetMultimerReferences();

      int   getMaxEType    ( mmdb::realtype  _entropy,
                             mmdb::realtype  dissThreshold,
                             int             r,
                             mmdb::ivector   eType,
                             mmdb::ivector   iflag,
                             PPInterface     interface,
                             int             nInterfaces );

      void  CalcInternalInterfaces ( int nInterfaces);

      void  SetMultimer    ( PMultimer     U,
                             int           nInterfaces,
                             PPInterface   interface,
                             mmdb::ivector eType );

      bool  isMacroMolecule();
      int   getMMSize      ();

      bool  isComplemenatry ( PMultimer U );

      bool  isIdentical     ( PMultimer U,  PDomains Domains,
                              mmdb::mat44 &  rom,  bool withLigands );
      bool  isInternal      ( int  intfType, mmdb::ivector eType,
                              int nInterfaces );

      void  CalcProperties ( PDomains      D,
                             PPInterface   interface,
                             mmdb::ivector eType,
                             int           nInterfaces,
                             int           nETypes,
                             PSymNumber    SNum,
                             mmdb::mat44 & rom );

      void  CalcTMatrices  ( mmdb::rpmat44 TMatrix, mmdb::mat44 & rom );

      void  Dissociate     ( PDomains      D,
                             PPInterface   interface,
                             mmdb::ivector eType,
                             int           nInterfaces,
                             int           nETypes,
                             PSymNumber    SNum,
                             mmdb::mat44 & rom );

      int   EngageInterface ( int interfaceNo, mmdb::omatrix  parMon,
                              mmdb::realtype intArea );

      void  getEntropyTerms ( mmdb::realtype & const_ent,
                              mmdb::realtype & rbody_ent,
                              mmdb::realtype & surf_ent );

      void GetDissProperties ( PPInterface      interface,
                               int              nInterfaces,
                               mmdb::rvector    DeltaSAS, // [0..nASPs-1]
                               mmdb::realtype & DeltaG,
                               int            & nHBonds,
                               int            & nSBridges );

    protected :

      void calcRBETerms  ( int              n,
                           mmdb::realtype & logW,
                           mmdb::realtype & logJ,
                           RDIRecursion     P );

      void calcRBEntropy ( int nDiss, mmdb::realtype & entropy_change,
                           RDIRecursion P );

      void dissInterface ( int              typeNo,
                           mmdb::ivector    itp,
                           RDIRecursion     P,
                           mmdb::realtype & dissEnX,
                           mmdb::realtype & entropyX,
                           mmdb::realtype & dissIntAreaX,
                           mmdb::ivector    dissMultX );

  };


  // ========================  IntfRecStack  =======================

  DefineClass(IntfRecStack);

  class IntfRecStack  {

    public :
      mmdb::ivector iflag;   //!< interface engagement flags
      mmdb::ivector mSize;   //!< sizes of the multimers
      PPMultimer        U;   //!< list of multimers
      int      nMultimers;   //!< number of multimers

      IntfRecStack ( int nETypes, int maxNofMultimers );
      ~IntfRecStack();

  };


  //  ========================  MultimerSet  ========================

  DefineClass(MultimerSet);

  class MultimerSet  {

    public :
      mmdb::ivector intf;  //!< vector of engaged interfaces
      PPMultimer       U;  //!< list of multimers obtained
      int     nMultimers;  //!< number of multimers obtained
      int  maxStableSize;  //!< maximal size of stable multimer
      bool     equiv_all;  //!< True if all multimers are equivalent
      bool    stable_all;

      MultimerSet ();
      ~MultimerSet();

      bool isEqual ( PIntfRecStack recStack, int nInterfaces );
      int     getMaxEType ( int            r,
                            mmdb::realtype dissThreshold,
                            mmdb::ivector  eType,
                            PPInterface    interface,
                            int            nInterfaces );

      void    makeSet     ( PIntfRecStack recStack,
                            int           nETypes,
                            int           nInterfaces,
                            int           nMols,
                            PPInterface   interface,
                            PDomains      D,
                            mmdb::mat44 & rom,
                            mmdb::ivector eType );

      void CalcProperties ( PDomains      D,
                            PPInterface   interface,
                            int           nInterfaces,
                            mmdb::ivector eType,
                            int           nETypes,
                            PSymNumber    SNum,
                            mmdb::mat44 & rom );


      int getMaxMultimerSize();

    protected :
      void FreeMemory();

  };

}  // namespace pisa

#endif
