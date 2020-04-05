// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.


// ************************************************************************* //
//                            avtSlugCodeFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_SlugCode_FILE_FORMAT_H
#define AVT_SlugCode_FILE_FORMAT_H

#include <avtSTSDFileFormat.h>

#include <hdf5.h>

// ****************************************************************************
//  Class: avtSlugCodeFileFormat
//
//  Purpose:
//      Reads in SlugCode files as a plugin to VisIt.
//
//  Programmer: ylee -- generated by xml2avt
//  Creation:   Fri Feb 7 15:34:32 PST 2020
//
// ****************************************************************************

class avtSlugCodeFileFormat : public avtSTSDFileFormat
{
  public:
                       avtSlugCodeFileFormat(const char *filename);
    virtual           ~avtSlugCodeFileFormat();

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    // virtual void      *GetAuxiliaryData(const char *var, const char *type,
    //                                  void *args, DestructorFunction &);
    //

    //
    // These are used to declare what the current time and cycle are for the
    // file.  These should only be defined if the file format knows what the
    // time and/or cycle is.
    //
    // virtual int       GetCycle(void);
    // virtual double    GetTime(void);
    //

    virtual const char    *GetType(void)   { return "SlugCode"; };
    virtual void           FreeUpResources(void);
    virtual void           ActivateTimestep(void);

    virtual vtkDataSet    *GetMesh(const char *);
    virtual vtkDataArray  *GetVar(const char *);
    virtual vtkDataArray  *GetVectorVar(const char *);

  protected:
    void Initialize();
    void InitializeHDF5();
    void FinalizeHDF5();
    void GetSimInfo();

    // DATA MEMBERS
    bool                      initialized;
    static int                objcnt;
    std::string               filename;
    int                       dimension;
    hid_t                     file_id;
    std::vector<std::string>  varNames;

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);
};


#endif
