// $Id: chem_equilibrium.h$
// =================================================================
//
//    14.08.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  chem_equilibrium <interface>
//       ~~~~~~~~~
//  **** Project :  Protein interfaces
//       ~~~~~~~~~
//  **** Classes :  chem::Species
//       ~~~~~~~~~  chem::Product
//                  chem::Equilibrium
//
//  (C) E. Krissinel, 2013
//
// =================================================================
//

#ifndef __CHEM_Equlibrium__
#define __CHEM_Equlibrium__

namespace chem  {

  class   Species;
  class   Product;

  typedef Species * PSpecies;
  typedef Product * PProduct;


  // =========================  Equilibrium  =========================

  class Equilibrium  {

    public:

      Equilibrium();
      virtual ~Equilibrium();

      void clear    ();
      void allocate ( int nOfSpecies, int nOfProducts );

      void addSpecies         ( int productNo, int speciesNo, int n=1 );
      void setKdiss           ( int productNo, long double Kdiss );
      void setLnKdiss         ( int productNo, double lnKdiss );
      void setStoichiometry   ( int speciesNo, double alpha );
      void setIterationLimit  ( int iterLimit );
      void prepareSolution    ();
      bool equilibrate        ( double totalMass,
                                double lnPrecision=1.0e-6 );
      long double Kdiss       ( int productNo );
      int         productSize ( int productNo );
      long double speciesConc ( int speciesNo );
      long double productConc ( int productNo );

    protected:
      PProduct *products;
      PSpecies *species;
      int       nProducts;
      int       nSpecies;
      int       iterationLimit;

    private:
      int nProdAlloc;
      int nSpecAlloc;

      long double getMass ( int speciesNo );
      bool        bracket ( int speciesNo, double totalMass,
                            double deltaLnC );

  };


}


#endif
