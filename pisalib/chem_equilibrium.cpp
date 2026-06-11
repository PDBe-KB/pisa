// $Id: chem_equilibrium.cpp$
// =================================================================
//
//    15.03.19   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  chem_equilibrium <implementation>
//       ~~~~~~~~~
//  **** Project :  Protein interfaces
//       ~~~~~~~~~
//  **** Classes :  chem::Species
//       ~~~~~~~~~  chem::Product
//                  chem::Equilibrium
//
//  (C) E. Krissinel, 2013-2019
//
// =================================================================
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "chem_equilibrium.h"

namespace chem  {

  // ===========================  Species  ===========================

  class Species  {

    public:
      double alpha;  //!< stoichiometry coefficient
      double   lnC;  //!< logarithm concentration
      bool    used;  //!< temporary field for equilibration procedure

      Species();
      virtual ~Species();

      void clear();

  };


  // ===========================  Product  ===========================

  class Product  {

    public:
      int         *species;  //!< species numbers in the product
      long double  Kd;       //!< dissociation constant
      long double  lnKd;     //!< logarithm dissociation constant
      int          size;     //!< total number of species in the product

      Product ( int nSpecies );
      virtual ~Product();

      void clear();

  };

}


// ===========================  Species  ===========================

chem::Species::Species()  {
  alpha = 1.0;
  lnC   = 0.0;
}

chem::Species::~Species()  {
  clear();
}

void chem::Species::clear()  {
}


// ===========================  Product  ===========================

chem::Product::Product ( int nSpecies )  {
int i;
  species = new int[nSpecies];
  for (i=0;i<nSpecies;i++)
    species[i] = 0;
  size = 0;
}

chem::Product::~Product()  {
  clear();
}

void chem::Product::clear()  {
  if (species)
    delete[] species;
}


// =========================  Equilibrium  =========================

chem::Equilibrium::Equilibrium()  {
  products       = NULL;
  species        = NULL;
  nProducts      = 0;
  nSpecies       = 0;
  nProdAlloc     = 0;
  nSpecAlloc     = 0;
  iterationLimit = 50;
}

chem::Equilibrium::~Equilibrium()  {
  clear();
}

void chem::Equilibrium::clear()  {
int i;

  for (i=0;i<nProdAlloc;i++)
    delete products[i];
  for (i=0;i<nSpecAlloc;i++)
    delete species[i];

  products   = NULL;
  species    = NULL;
  nProducts  = 0;
  nSpecies   = 0;
  nProdAlloc = 0;
  nSpecAlloc = 0;

}

void chem::Equilibrium::setIterationLimit ( int iterLimit )  {
  iterationLimit = iterLimit;
}

void chem::Equilibrium::allocate ( int nOfSpecies, int nOfProducts ) {
int i;

  clear();

  nSpecAlloc = nOfSpecies;
  nSpecies   = nOfSpecies;
  nProdAlloc = nOfProducts;
  nProducts  = nOfProducts;

  products = new PProduct[nProdAlloc];
  species  = new PSpecies[nSpecAlloc];

  for (i=0;i<nProdAlloc;i++)
    products[i] = new Product ( nSpecAlloc );

  for (i=0;i<nSpecAlloc;i++)
    species[i] = new Species();

}

void chem::Equilibrium::addSpecies ( int productNo, int speciesNo,
                                     int n )  {
  products[productNo]->species[speciesNo] += n;
  products[productNo]->size += n;
}

void chem::Equilibrium::setKdiss ( int productNo, long double Kdiss )  {
  products[productNo]->Kd   = Kdiss;
  products[productNo]->lnKd = logl(Kdiss);
}

void chem::Equilibrium::setLnKdiss ( int productNo, double lnKdiss )  {
  products[productNo]->Kd   = expl ( lnKdiss );
  products[productNo]->lnKd = lnKdiss;
}


void chem::Equilibrium::setStoichiometry ( int speciesNo,
                                           double alpha )  {
  species[speciesNo]->alpha = alpha;
}

void chem::Equilibrium::prepareSolution()  {
//double s;
int    i;

//  s = 0.0;
//  for (i=0;i<nSpecies;i++)
//    s += species[i]->alpha;

  for (i=0;i<nSpecies;i++)  {
//    species[i]->alpha /= s;
    species[i]->lnC    = 1.0e-12*species[i]->alpha;
  }

}

long double chem::Equilibrium::getMass ( int speciesNo )  {
Product     *p;
int         *s;
long double  C;
double       lnC;
int          i,j,n;

  C = expl ( species[speciesNo]->lnC );
  for (i=0;i<nProducts;i++)  {
    p = products[i];
    s = p->species;
    n = s[speciesNo];
    if (n>0)  {
      lnC = -p->lnKd;
      for (j=0;j<nSpecies;j++)
        lnC += s[j]*species[j]->lnC;
      C += n*expl ( lnC );
    }
  }

  return C;

}

bool chem::Equilibrium::bracket ( int speciesNo, double totalMass,
                                  double deltaLnC )  {
// Returns true if already bracketed
long double mass1,M,dmass1,dmass2;
double      lnC1,lnC2;
int         n;

  M = static_cast<long double>(totalMass) * species[speciesNo]->alpha;

  lnC1  = species[speciesNo]->lnC;
  lnC2  = lnC1;
  mass1 = getMass ( speciesNo );
  n     = 0;

  dmass1 = M - mass1;
  dmass2 = dmass1;  // only to please compiler

  if (mass1>=M)  {

    while (mass1>=M)  {
      dmass2 = dmass1;
      lnC2   = lnC1;
      lnC1  -= deltaLnC;
      species[speciesNo]->lnC = lnC1;
      mass1  = getMass ( speciesNo );
      dmass1 = M - mass1;
      n++;
    }

  } else  {

    while (mass1<M)  {
      dmass2 = dmass1;
      lnC2   = lnC1;
      lnC1  += deltaLnC;
      species[speciesNo]->lnC = lnC1;
      mass1  = getMass ( speciesNo );
      dmass1 = M - mass1;
      n++;
    }

  }

//  in this case, this relaxation causes instability
//  species[speciesNo]->lnC = (lnC1+lnC2)/2.0;

//  if (dmass1<dmass2)  species[speciesNo]->lnC = lnC1;
//                else  species[speciesNo]->lnC = lnC2;

  if (dmass1!=dmass2)
    species[speciesNo]->lnC = lnC1 + (lnC2-lnC1)*dmass1/(dmass1-dmass2);
  else
    species[speciesNo]->lnC = lnC2;

  return (n<=1);

}


bool chem::Equilibrium::equilibrate ( double totalMass,
                                      double lnPrecision ) {
double deltaLnC, lnCi;
int    i,j,k,iter;
bool   b,bracketed, converged;

  iter     = 0;    // iteration counter
  deltaLnC = 0.1;  // initial bracket in log scale

  do {

    bracketed = true;
    converged = true;

    for (i=0;i<nSpecies;i++)
      species[i]->used = false;

    for (k=0;k<nSpecies;k++)  {

      lnCi = 1.0e10;
      i    = -1;
      for (j=0;j<nSpecies;j++)
        if ((!species[j]->used) && (species[j]->lnC<=lnCi))  {
          lnCi = species[j]->lnC;
          i    = j;
        }

      if (i>=0)  {
//        lnCi = species[i]->lnC;
        b = bracket ( i,totalMass,deltaLnC );
        if ((k>0) && (!b))
          bracketed = false;
        if (fabs(lnCi-species[i]->lnC)>lnPrecision)
          converged = false;
        species[i]->used = true;
      }

    }

    if (bracketed)
      deltaLnC /= sqrt(2.0);

    iter++;

  } while ((iter<iterationLimit) &&
           ((2.0*deltaLnC>lnPrecision) || (!converged)));

  return (iter<iterationLimit);

}

long double chem::Equilibrium::Kdiss ( int productNo )  {
  return products[productNo]->Kd;
}

int chem::Equilibrium::productSize ( int productNo )  {
  return products[productNo]->size;
}

long double chem::Equilibrium::speciesConc ( int speciesNo )  {
  return expl ( species[speciesNo]->lnC );
}

long double chem::Equilibrium::productConc ( int productNo )  {
Product *p;
int     *s;
double   lnC;
int      i;

  p   = products[productNo];
  s   = p->species;
  lnC = -p->lnKd;
  for (i=0;i<nSpecies;i++)
    lnC += s[i]*species[i]->lnC;

  return expl ( lnC );

}
