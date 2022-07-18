// $Id: pisa_symnumber.cpp $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_symnumber <implementation>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  pisa::SymNumber
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2007-2014
//
// =================================================================
//

#include "pisa_symnumber.h"
#include "pisa_engunits.h"
#include "pisa_assemblies.h"
#include "pisa_types.h"

#include "mmdb2/mmdb_math_linalg.h"

namespace pisa  {

  // =========================  SymNumber  ==========================

  SymNumber::SymNumber() {
    C   = NULL;
    F   = NULL;
    allocSize = 0;
    A   = NULL;
    U   = NULL;
    V   = NULL;
    W   = NULL;
    RV1 = NULL;
  }

  SymNumber::~SymNumber() {
    FreeMemory();
    mmdb::FreeVectorMemory ( RV1,1 );
    mmdb::FreeVectorMemory ( W  ,1 );
    mmdb::FreeMatrixMemory ( V  ,3,1,1 );
    mmdb::FreeMatrixMemory ( U  ,3,1,1 );
    mmdb::FreeMatrixMemory ( A  ,3,1,1 );
  }

  void SymNumber::FreeMemory()  {
    mmdb::FreeMatrixMemory ( C,allocSize,0,0 );
    if (F)  {
      delete[] F;
      F = NULL;
    }
    allocSize = 0;
  }

  void SymNumber::getStorage ( mmdb::rmatrix & S1, mmdb::rmatrix & S2,
                               mmdb::rvector & S3, mmdb::rvector & S4 ) {
    if (!A)  {
      mmdb::GetMatrixMemory ( A  ,3,3,1,1 );
      mmdb::GetMatrixMemory ( U  ,3,3,1,1 );
      mmdb::GetMatrixMemory ( V  ,3,3,1,1 );
      mmdb::GetVectorMemory ( W  ,3,1 );
      mmdb::GetVectorMemory ( RV1,3,1 );
    }
    S1 = A;
    S2 = U;
    S3 = W;
    S4 = RV1;
  }

  int SymNumber::getSymNumber ( PMultimer M, int n, PDomains D,
                                mmdb::mat44 * TMatrix )  {
  //   M       multimeric unit
  //   n       multimeric subunit number n, identified by index id2
  //           in monomeric units
  //   TMatrix transformation matrix for coordinates of ith monomeric
  //           unit in multimer M
  int  i,j,i1,j1,dId,dNo;
  bool B;

    dId = 0;  // only to keep compiler happy

    if (M->mSize<=1)
      return D->domain[M->R[0].M->id]->symNumber;

    if (M->mSize>allocSize)  {
      FreeMemory();
      allocSize = M->mSize;
      mmdb::GetMatrixMemory ( C,allocSize,allocSize,0,0 );
      F = new TFrame[allocSize];
    }

    // Fill the connectivity matrix C:
    //
    //    C[i][j] = N > 0 <=> successful trial number N
    //    C[i][j] = 0     <=> equivalent units
    //    C[i][j] = -1    <=> non-equivalent units
    //    C[i][j] = -2    <=> unsuccessful trial
    //
    i1 = 0;
    B  = true;
    for (i=0;(i<M->mSize) && B;i++)
      if ((M->R[i].M->dclass!=DCLASS_Ligand) && ((!n) ||
          (M->R[i].M->id2==n))) {
        dId = M->R[i].M->id;
        C[i1][i1] = 0;
        j1 = i1+1;
        B  = false;
        for (j=i+1;j<M->mSize;j++)
          if ((M->R[j].M->dclass!=DCLASS_Ligand) &&
              ((!n) || (M->R[j].M->id2==n)))  {
            if (D->domain[M->R[i].M->id]->type ==
                D->domain[M->R[j].M->id]->type)  {
              C[i1][j1] = 0;
              B = true;
            } else
              C[i1][j1] = -1;
            C[j1][i1] = C[i1][j1];
            j1++;
          }
        if (!B)
          for (j1=0;(j1<i1) && (!B);j1++)
            B = (C[i1][j1]==0);
        i1++;
      }

    if (i1==1)  // only one subunit is there
      return D->domain[dId]->symNumber;

    if (!B)  return 1;

    // Calculate coordinates frames F[i] of units transformed
    // to their actual orientation in the multimer
    i1 = 0;
    for (i=0;i<M->mSize;i++)
      if ((M->R[i].M->dclass!=DCLASS_Ligand) && ((!n) ||
          (M->R[i].M->id2==n)))  {
        dNo = D->domain[M->R[i].M->id]->ncsParent;
        D->domain[dNo]->getReferenceFrame ( F[i1],TMatrix[i] );
        i1++;
      }

    nUnits = i1;

    return calcSymNumber();

  }


  int SymNumber::getSymNumber ( PAssembly Asm, PDomains D )  {
  int  i,j,i1,j1,dNo;
  bool B;

    if (Asm->mmSize<=1)
      return D->domain[Asm->M[0]->id]->symNumber;

    if (Asm->mmSize>allocSize)  {
      FreeMemory();
      allocSize = Asm->mmSize;
      mmdb::GetMatrixMemory ( C,allocSize,allocSize,0,0 );
      F = new TFrame[allocSize];
    }

    nUnits = Asm->mmSize;

    // Fill the connectivity matrix C:
    //
    //    C[i][j] = N > 0 <=> successful trial number N
    //    C[i][j] = 0     <=> equivalent units
    //    C[i][j] = -1    <=> non-equivalent units
    //    C[i][j] = -2    <=> unsuccessful trial
    //
    B = true;
    i = 0;
    for (i1=0;(i1<Asm->asmSize) && B;i1++)
      if (D->domain[Asm->M[i1]->id]->dclass!=DCLASS_Ligand)  {
        C[i][i] = 0;
        B = false;
        j = i+1;
        for (j1=i1+1;j1<Asm->asmSize;j1++)
          if (D->domain[Asm->M[j1]->id]->dclass!=DCLASS_Ligand)  {
            if (D->domain[Asm->M[i1]->id]->type ==
                D->domain[Asm->M[j1]->id]->type)  {
              C[i][j] = 0;
              B = true;
            } else
              C[i][j] = -1;
            C[j][i] = C[i][j];
            j++;
          }
        if (!B)
          for (j=0;(j<i) && (!B);j++)
            B = (C[i][j]==0);
        i++;
      }

    if (!B)  return 1;

    // Calculate coordinates frames F[i] of units transformed
    // to their actual orientation in the multimer
    i = 0;
    for (i1=0;i1<Asm->asmSize;i1++)
      if (D->domain[Asm->M[i1]->id]->dclass!=DCLASS_Ligand)  {
        dNo = D->domain[Asm->M[i1]->id]->ncsParent;
        D->domain[dNo]->getReferenceFrame ( F[i],Asm->M[i1]->T );
        i++;
      }

    return calcSymNumber();

  }


  int  SymNumber::superposeFrames ( mmdb::mat44 & T,
                                    mmdb::ivector f1,
                                    mmdb::ivector f2,
                                    int nFrames )  {
  //  Calculates matrix T such that T*F[f1] -> F[f2]
  mmdb::realtype xc1,yc1,zc1, xc2,yc2,zc2, det,B;
  mmdb::vect3    vc1,vc2;
  int            i,j,k,n,i1,i2;

    //  1. Calculate mass centers

    xc1 = 0.0;
    yc1 = 0.0;
    zc1 = 0.0;
    xc2 = 0.0;
    yc2 = 0.0;
    zc2 = 0.0;

    for (i=0;i<nFrames;i++)  {
      i1 = f1[i];
      i2 = f2[i];
      for (j=0;j<frameLen;j++)  {
        xc1 += F[i1][j][0];
        yc1 += F[i1][j][1];
        zc1 += F[i1][j][2];
        xc2 += F[i2][j][0];
        yc2 += F[i2][j][1];
        zc2 += F[i2][j][2];
      }
    }

    k = frameLen*nFrames;
    xc1 /= k;
    yc1 /= k;
    zc1 /= k;
    xc2 /= k;
    yc2 /= k;
    zc2 /= k;

    //  2.  Calculate the correlation matrix

    if (!A)  {
      mmdb::GetMatrixMemory ( A  ,3,3,1,1 );
      mmdb::GetMatrixMemory ( U  ,3,3,1,1 );
      mmdb::GetMatrixMemory ( V  ,3,3,1,1 );
      mmdb::GetVectorMemory ( W  ,3,1 );
      mmdb::GetVectorMemory ( RV1,3,1 );
    }

    for (i=1;i<=3;i++)
      for (j=1;j<=3;j++)
        A[i][j] = 0.0;

    for (n=0;n<nFrames;n++)  {
      i1 = f1[n];
      i2 = f2[n];
      for (k=0;k<frameLen;k++)  {
        vc1[0] = F[i1][k][0] - xc1;
        vc1[1] = F[i1][k][1] - yc1;
        vc1[2] = F[i1][k][2] - zc1;
        vc2[0] = F[i2][k][0] - xc2;
        vc2[1] = F[i2][k][1] - yc2;
        vc2[2] = F[i2][k][2] - zc2;
        for (i=1;i<=3;i++)
          for (j=1;j<=3;j++)
            A[i][j] += vc1[j-1]*vc2[i-1];
      }
    }


    //  3. Calculate transformation matrix (to be applied to f1)

    det = A[1][1]*A[2][2]*A[3][3] +
          A[1][2]*A[2][3]*A[3][1] +
          A[2][1]*A[3][2]*A[1][3] -
          A[1][3]*A[2][2]*A[3][1] -
          A[1][1]*A[2][3]*A[3][2] -
          A[3][3]*A[1][2]*A[2][1];

    //  3.1 SV-decompose the correlation matrix

    mmdb::math::SVD ( 3,3,3,A,U,V,W,RV1,true,true,i );

    if (i!=0)   return i;  // SVD fail

    //  3.2 Check for parasite inversion and fix it if found

    if (det<=0.0)  {
      k = 0;
      B = mmdb::MaxReal;
      for (j=1;j<=3;j++)
        if (W[j]<B)  {
          B = W[j];
          k = j;
        }
      for (j=1;j<=3;j++)
        V[j][k] = -V[j][k];
    }


    //  3.3 Calculate rotational part of T

    for (j=1;j<=3;j++)
      for (k=1;k<=3;k++)  {
        B = 0.0;
        for (i=1;i<=3;i++)
          B += U[j][i]*V[k][i];
        T[j-1][k-1] = B;
      }


    //  3.4 Add translational part to T

    T[0][3] = xc2 - T[0][0]*xc1 - T[0][1]*yc1 - T[0][2]*zc1;
    T[1][3] = yc2 - T[1][0]*xc1 - T[1][1]*yc1 - T[1][2]*zc1;
    T[2][3] = zc2 - T[2][0]*xc1 - T[2][1]*yc1 - T[2][2]*zc1;

    return 0;

  }


  int SymNumber::calcSymNumber()  {
  mmdb::mat44    T;
  TFrame         tf;
  mmdb::ivector  f1,f2;
  mmdb::realtype dx,dy,dz, d,dmin, thres1,thres2;
  int            n,i,j,k,p,q, kmin, key;

    mmdb::GetVectorMemory ( f1,nUnits,0 );
    mmdb::GetVectorMemory ( f2,nUnits,0 );

    thres1 = 6.250*frameLen;
    thres2 = 225.0*frameLen;

    n = 1;  // symmetry number

    for (i=1;i<nUnits;i++)
      if (!C[0][i])  {
        //   There is a proper and yet untried rotation which
        // brings unit 0 into position of an equivalent unit i.
        // If this rotation brings all other units into positions
        // of their equivalents, then the rotation is counted
        // into the symmetry number.
        //   First find the rotation.
        f1[0] = 0;
        f2[0] = i;
        if (!superposeFrames(T,f1,f2,1))  {
          // now check that matrix T brings every other unit
          // into position of its equivalent
          key = 1;
          for (j=1;(j<nUnits) && (key>0);j++) {
            // get coordinate frame of transformed unit j in tf
            for (k=0;k<frameLen;k++)
              for (p=0;p<3;p++)  {
                tf[k][p] = T[p][3];
                for (q=0;q<3;q++)
                  tf[k][p] += T[p][q]*F[j][k][q];
              }
            // check that tf coincides with coordinate frame
            // of an equivalent unit
            dmin = mmdb::MaxReal;
            kmin = -1;
            for (k=0;k<nUnits;k++)
              if ((!C[j][k]) && (k!=j))  {
                d = 0.0;
                for (p=0;p<frameLen;p++)  {
                  dx = tf[p][0] - F[k][p][0];
                  dy = tf[p][1] - F[k][p][1];
                  dz = tf[p][2] - F[k][p][2];
                  d += dx*dx + dy*dy + dz*dz;
                }
                if (d<dmin)  {
                  dmin = d;
                  kmin = k;
                }
              }
            if (dmin>thres2)   key = 0;  // no solution
            else  {
              f1[j] = j;
              f2[j] = kmin;
              if (dmin>thres1)  key = 2;  // to be checked again
            }
          }
          if (key>1)  {
            //   This may be a self-superposing rotation, however
            // the superposition is not perfect when rotation matrix
            // is calculated using units 0 and i only. Check whether
            // a corrected matrix, taking into account all matched
            // frames, makes a better result, and if it does then
            // count it in (assign key=1)
            if (!superposeFrames(T,f1,f2,nUnits))  {
              key = 1;
              for (j=1;(j<nUnits) && (key>0);j++)  {
                // get coordinate frame of transformed unit j in tf
                for (k=0;k<frameLen;k++)
                  for (p=0;p<3;p++)  {
                    tf[k][p] = T[p][3];
                    for (q=0;q<3;q++)
                      tf[k][p] += T[p][q]*F[j][k][q];
                }
                // check that tf coincides with coordinate frame
                // of an equivalent unit
                k = f2[j];
                d = 0.0;
                for (p=0;p<frameLen;p++)  {
                  dx = tf[p][0] - F[k][p][0];
                  dy = tf[p][1] - F[k][p][1];
                  dz = tf[p][2] - F[k][p][2];
                  d += dx*dx + dy*dy + dz*dz;
                }
                if (d>thres1)  key = 0; // no solution
              }
            }
          }
          if (key==1)  n++;  // increase the symmetry number
        } else
          C[0][i] = -3;  // error - should never happen
      }

    mmdb::FreeVectorMemory ( f1,0 );
    mmdb::FreeVectorMemory ( f2,0 );

    return n;

  }

}  // namespace pisa
