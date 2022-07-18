// $Id: pisa_datastats.h $
// =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  ----------------------------------------------------------------
//
//  **** Module  :  pisa_datastats <interface>
//       ~~~~~~~~~
//  **** Project :  PISA
//       ~~~~~~~~~
//  **** Classes :  DataStats
//       ~~~~~~~~~
//
//  (C) E. Krissinel, 2013-2014
//
// =================================================================
//

#ifndef __PISA_DATASTATS__
#define __PISA_DATASTATS__

#include "mmdb2/mmdb_mattype.h"
#include "mmdb2/mmdb_io_file.h"


namespace pisa  {

  //  =====================  StatDistribution  =====================

  DefineClass(StatDistribution);

  class StatDistribution  {
    // distribution over one parameter

    public :

      StatDistribution ();
      ~StatDistribution();

      void clear();

      void addData   ( mmdb::realtype v, bool inAssembly );
      void makeParts ( int maxNParts    );

      mmdb::realtype getStats ( mmdb::realtype v );

      mmdb::realtype getMinValue();
      mmdb::realtype getMaxValue();

      void print();

      void write ( mmdb::io::File & f );
      void read  ( mmdb::io::File & f );

    protected:
      mmdb::shortreal *parts;
      mmdb::shortreal *prob;
      mmdb::ovector    in_assembly;
      int              nParts;

    private:
      int              nPartsAlloc;

  };


  //  =======================  StatDistSet  ========================

  DefineClass(StatDistSet);

  class StatDistSet  {
    // set of distributions at fixed scan parameter

    public :

      StatDistSet ();
      ~StatDistSet();

      void clear();

      void initSet ( int n_params, int param_no );

      void addData   ( mmdb::rvector p, bool inAssembly );
      void makeParts ( int   maxNParts );

      void getStats ( mmdb::rvector stats );

      inline int nParameters()  { return nParams; }
      inline PStatDistribution getDistribution ( int paramNo )
                                { return SD[paramNo]; }

      void write ( mmdb::io::File & f );
      void read  ( mmdb::io::File & f );

    protected:
      PPStatDistribution SD;
      int                nParams;

  };


  //  =========================  DataScan  =========================

  DefineClass(DataScan);

  class DataScan  {
    // scan over a given parameter

    public :

      DataScan ();
      ~DataScan();

      void clear();

      void initScan ( int n_parts, mmdb::rvector parts,
                      int nParams, int paramNo );

      int  getPartNo ( mmdb::realtype value );
      void addData   ( mmdb::rvector p, bool inAssembly );
      void makeParts ( int   maxNParts );

      inline PStatDistSet getDistSet ( int n ) { return distSet[n]; }

      inline int nPartitions()                 { return nParts;    }
      inline mmdb::shortreal *getPartitions()  { return partition; }

      void printPartition();

      void write ( mmdb::io::File & f );
      void read  ( mmdb::io::File & f );

    protected:
      PPStatDistSet    distSet;
      mmdb::shortreal *partition;
      int              nParts;
      int              scanParam;

  };


  //  ========================  DataStats  =========================

  DefineClass(DataStats);

  class DataStats  {
    // scan over all parameters

    public :

      enum PROC_KEY { Raw,Minus,MinusLog };

      DataStats ();
      ~DataStats();

      void clear();

      void initStats    ( int n_parameters );
      void setDataName  ( mmdb::cpstr name );
      void setParamID   ( int paramNo, mmdb::cpstr ID   );
      void setParamName ( int paramNo, mmdb::cpstr name );
      void initScan     ( int paramNo, mmdb::rvector partition,
                          int nPartitions );
      void setProcKey   ( int paramNo, PROC_KEY key );

      void addData      ( mmdb::rvector p, bool inAssembly );
      void makeParts    ( int   maxNParts );
      void addData      ( mmdb::rvector p, bool inAssembly,
                          RStatDistSet D );

      inline int nParameters()          { return nParams;  }
      inline mmdb::cpstr getDataName () { return dataName; }
      inline mmdb::cpstr getParamID  ( int paramNo )
                                        { return paramID[paramNo];   }
      inline mmdb::cpstr getParamName ( int paramNo )
                                        { return paramName[paramNo]; }
      inline mmdb::psvector getParamNames() { return paramName; }
      PDataScan getDataScan ( int paramNo );

      int  getStats ( int scanParamNo, mmdb::rvector stats,
                      PStatDistSet D );

      bool write ( mmdb::cpstr fileName );
      bool read  ( mmdb::cpstr fileName );

      void write ( mmdb::io::File & f );
      bool read  ( mmdb::io::File & f );

    protected:
      StatDistSet     SDS;        //!< unreferenced stats
      mmdb::pstr      dataName;
      PPDataScan      scan;
      mmdb::psvector  paramID;   //!< short name for radar plot
      mmdb::psvector  paramName; //!< long name for combobox etc
      mmdb::ivector   procKey;
      mmdb::rvector   buf;
      int             nParams;

      void convertParams ( mmdb::rvector p );

  };

}

#endif
