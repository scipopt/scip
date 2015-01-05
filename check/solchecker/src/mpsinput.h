/**
 * @file mpsinput.h
 * @brief MPS reader class
 *
 * @author Domenico Salvagnin
 */

#ifndef MPSINPUT_H
#define MPSINPUT_H

#include <string>
#include <stdexcept>
#include <stdio.h>
#include <zlib.h>

#define MPS_MAX_LINELEN  1024

class Model;

/**
 * @brief Exception class to store syntax error information
 */

class SyntaxError : public std::exception 
{
public:
   /**
    * Constructor: uses the line and the line number to build the error message
    */
   SyntaxError(const char* line, int lineno);
   ~SyntaxError() throw() {}
   const char* what() const throw();
protected:
   std::string msg; //< error message string
};

/**
 * @brief MPS Reader Class.
 * It reads an MPS file (possibly gzipped), both in the old fixed length format
 * and in the new "relaxed" format. It handles all basic extensions.
 * It does NOT handle SOS and quadratic stuff yet!
 */

class MpsInput
{
public:
   /** default constructor */
   MpsInput();
   /**
    * @brief reads an MPS file
    * @param _filename path to the MPS file (may be gzipped)
    * @param _model model to fill with the data read from file
    * @return true if the read was successful, false otherwise
    */
   bool readMps(const char* _filename, Model* _model);
protected:
   /**
    * Section of the MPS format we are into
    */
   enum MpsSection
   {
      MPS_NAME,
      MPS_OBJSEN,
      MPS_OBJNAME,
      MPS_ROWS,
      MPS_USERCUTS,
      MPS_LAZYCONS,
      MPS_COLUMNS,
      MPS_RHS,
      MPS_RANGES,
      MPS_BOUNDS,
      MPS_SOS,
      MPS_ENDATA
   };
   
   /** @name data */
   // *{
   MpsSection    section; //< current section of the MPS section
   Model*        model; //< model to fill
   FILE*         fp; //< FILE pointer for plain MPS files
   gzFile        gzfp; //< gzFile pointer for gzipped files
   bool          isZipped; //< are we reading from a gzipped file?
   int           lineno; //< current line number
   char          buf[MPS_MAX_LINELEN]; //< buffer to store and manipulate each line of input
   const char*   f0; //< field 0 of an MPS line
   const char*   f1; //< field 1 of an MPS line
   const char*   f2; //< field 2 of an MPS line
   const char*   f3; //< field 3 of an MPS line
   const char*   f4; //< field 4 of an MPS line
   const char*   f5; //< field 5 of an MPS line
   bool          isInteger; //< is the current variable in a integer block?
   bool          isNewFormat; //< is this a new format MPS file?
   bool          semiContWarning; //< was there any semicontinuous declaration?
   // *}

   /** @name helpers */
   // *{
   bool readLine();
   void insertName(const char* name, bool second);
   void readName();
   void readObjsen();
   void readObjname();
   void readRows();
   void readCols();
   void readRhs();
   void readRanges();
   void readBounds();
   void readSOS();
   // *}
};

#endif /* MPSINPUT_H */
