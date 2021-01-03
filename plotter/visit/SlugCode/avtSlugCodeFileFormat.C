// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

// ************************************************************************* //
//                            avtSlugCodeFileFormat.C                           //
// ************************************************************************* //

#include <avtSlugCodeFileFormat.h>

#include <string>

#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>

#include <avtDatabaseMetaData.h>

#include <DBOptionsAttributes.h>
#include <Expression.h>

#include <InvalidVariableException.h>
#include <InvalidDBTypeException.h>
#include <InvalidFilesException.h>

#include <DebugStream.h>

#include <hdf5.h>
#include <visit-hdf5.h>


using     std::string;


int avtSlugCodeFileFormat::objcnt = 0;

// ****************************************************************************
//  Method: avtSlugCodeFileFormat::ActivateTimestep
//
//  Purpose:
//      Tells the reader to activate the current time step.
//
//  Programmer: ylee
//
// ****************************************************************************
void
avtSlugCodeFileFormat::ActivateTimestep()
{
   Initialize();
}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::Initialize
//
//  Purpose:
//      It tries to open the file in HDF5 format and set 'initialized'
//      true if it success.
//
//  Programmer: ylee
//
// ****************************************************************************
void
avtSlugCodeFileFormat::Initialize()
{
    if(!initialized)
    {
        bool okay = false;
        // suppress "attempt to open an HDF5 file without H5F_CLOSE_SEMI" warnings
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fclose_degree(fapl, H5F_CLOSE_SEMI);
        file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);
        okay = true;
        // If your file format API could not open the file then throw // an exception.
        if (!okay)
        {
            EXCEPTION1(InvalidDBTypeException, "The file could not be opened");
        }
        initialized = true;
    }
}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::InitializeHDF5
//
//  Purpose:
//      It opens HDF5 interface.
//
//  Programmer: ylee
//
// ****************************************************************************
void
avtSlugCodeFileFormat::InitializeHDF5(void)
{
    H5open();
    debug1 << "HDF5 interface is opened" << endl;
}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::FinalizeHDF5
//
//  Purpose:
//      For enclosing HDF5 interface, it collects all HDF5-related garbages.
//      Implemented from FLASH plugin.
//
//  Programmer: ylee
//
// ****************************************************************************
void
avtSlugCodeFileFormat::FinalizeHDF5(void)
{
    H5garbage_collect();
    debug1 << "Garbage collected" << endl;
}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::GetSimInfo
//
//  Purpose:
//      It opens HDF5 file and determines dimensionality of the file.
//
//  Programmer: ylee
//
// ****************************************************************************
void
avtSlugCodeFileFormat::GetSimInfo()
{
    hid_t     fapl;

    // suppress "attempt to open an HDF5 file without H5F_CLOSE_SEMI" warnings
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_SEMI);
    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);
    if (file_id < 0)
    {
        EXCEPTION1(InvalidFilesException, filename.c_str());
    }

    // find dimensionality
    int check = H5Lexists(file_id, "/zCoord", H5P_DEFAULT);
    if (check)
    {
        dimension = 3;
    }
    else
    {
        dimension = 2;
    }

}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::GetCycle
//
//  Purpose:
//      It opens HDF5 file and returns current simulation step
//
//  Programmer: ylee
//
// ****************************************************************************
int
avtSlugCodeFileFormat::GetCycle()
{
    // open file temporarily
    // suppress "attempt to open an HDF5 file without H5F_CLOSE_SEMI" warnings
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_SEMI);
    hid_t tmp_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);
    if (tmp_file_id < 0)
        return INVALID_CYCLE;

    int cycle;

    hid_t dset_id = H5Dopen(tmp_file_id, "nStep", H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cycle);

    H5Dclose(dset_id);
    H5Fclose(tmp_file_id);
    return cycle;
}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::GetTime
//
//  Purpose:
//      It opens HDF5 file and returns current simulation time
//
//  Programmer: ylee
//
// ****************************************************************************
double
avtSlugCodeFileFormat::GetTime()
{
    // open file temporarily
    // suppress "attempt to open an HDF5 file without H5F_CLOSE_SEMI" warnings
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_SEMI);
    hid_t tmp_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);
    if (tmp_file_id < 0)
        return INVALID_TIME;

    double time;

    hid_t dset_id = H5Dopen(tmp_file_id, "time", H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &time);

    H5Dclose(dset_id);
    H5Fclose(tmp_file_id);
    return time;
}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat constructor
//
//  Programmer: ylee -- generated by xml2avt
//  Creation:   Fri Feb 7 15:34:32 PST 2020
//
// ****************************************************************************
avtSlugCodeFileFormat::avtSlugCodeFileFormat(const char *cfilename)
    : avtSTSDFileFormat(cfilename)
{
    // INITIALIZE DATA MEMBERS
    initialized = false;
    filename    = cfilename;
    file_id     = -1;
    dimension   = 0;

    // do HDF5 library initialization on consturction of first instance
    if (avtSlugCodeFileFormat::objcnt == 0)
        InitializeHDF5();
    avtSlugCodeFileFormat::objcnt++;
}


// ****************************************************************************
//  Function:  avtSlugCodeFileFormat::~avtSlugCodeFileFormat
//
//  Purpose:   Free up resources
//
//  Programmer:  ylee
//
// ****************************************************************************
avtSlugCodeFileFormat::~avtSlugCodeFileFormat()
{
    FreeUpResources();
}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::FreeUpResources
//
//  Purpose:
//      When VisIt is done focusing on a particular timestep, it asks that
//      timestep to free up any resources (memory, file descriptors) that
//      it has associated with it.  This method is the mechanism for doing
//      that.
//
//  Programmer: ylee -- generated by xml2avt
//  Creation:   Fri Feb 7 15:34:32 PST 2020
//
// ****************************************************************************
void
avtSlugCodeFileFormat::FreeUpResources(void)
{
    if (file_id >= 0)
    {
        H5Fclose(file_id);
        file_id = -1;
    }

    // handle HDF5 library termination on descrution of last instance
    avtSlugCodeFileFormat::objcnt--;
    if (avtSlugCodeFileFormat::objcnt == 0)
        FinalizeHDF5();
}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: ylee -- generated by xml2avt
//  Creation:   Fri Feb 7 15:34:32 PST 2020
//
// ****************************************************************************
void
avtSlugCodeFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md)
{

    GetSimInfo();

    // initialize mesh
    avtMeshMetaData *mesh = new avtMeshMetaData;
    mesh->name = "mesh";
    mesh->originalName = "mesh";

    mesh->meshType = AVT_RECTILINEAR_MESH;
    mesh->topologicalDimension = dimension;
    mesh->spatialDimension = dimension;
    mesh->blockOrigin = 0;

    md->Add(mesh);

    // add vriable names to mesh
    if (dimension > 2)  // 3D
        varNames = {"dens", "velx", "vely", "velz", "pres"};
    else                // 2D
        varNames = {"dens", "velx", "vely", "pres"};

    for (size_t v = 0 ; v < varNames.size(); v++)
    {
        AddScalarVarToMetaData(md, varNames[v], "mesh", AVT_ZONECENT);
    }

}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::GetMesh
//
//  Purpose:
//      Gets the mesh associated with this file.  The mesh is returned as a
//      derived type of vtkDataSet (ie vtkRectilinearGrid, vtkStructuredGrid,
//      vtkUnstructuredGrid, etc).
//
//  Arguments:
//      meshname    The name of the mesh of interest.  This can be ignored if
//                  there is only one mesh.
//
//  Programmer: ylee -- generated by xml2avt
//  Creation:   Fri Feb 7 15:34:32 PST 2020
//
// ****************************************************************************
vtkDataSet *
avtSlugCodeFileFormat::GetMesh(const char *meshname)
{
    GetSimInfo();

    int            rank;
    int            dims[3];

    hid_t          tmp_file_id, dset_id, dspace_id, fapl;
    hsize_t        count[1];
    vtkFloatArray  *coords[3];


    // open file temporarily
    // suppress "attempt to open an HDF5 file without H5F_CLOSE_SEMI" warnings
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_SEMI);
    tmp_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

    // xcoord
    // open dset and dspace
    dset_id = H5Dopen(tmp_file_id, "xCoord", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);

    // count number of nodes on axis
    H5Sget_simple_extent_dims(dspace_id, count, NULL);
    rank = count[0];
    dims[0] = rank+1;

    // read xcoord from HDF5
    float xcoord[rank], xcoord_ext[rank+1], dx;
    H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xcoord);
    // extend coordinate
    dx = abs(xcoord[1] - xcoord[0]);
    for (int i = 0; i < rank; i++)
    {
        xcoord_ext[i] = xcoord[i]-dx/2;
    }
    xcoord_ext[rank] = xcoord_ext[rank-1] + dx;

    coords[0] = vtkFloatArray::New();
    coords[0]->SetNumberOfTuples(rank+1);

    for (int i = 0; i < rank+1; i++)
    {
        coords[0]->SetComponent(i, 0, xcoord_ext[i]);
    }
    // close dset
    H5Dclose(dset_id);


    // ycoord
    // open dset and dspace
    dset_id = H5Dopen(tmp_file_id, "yCoord", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);

    // count number of nodes on axis
    H5Sget_simple_extent_dims(dspace_id, count, NULL);
    rank = count[0];
    dims[1] = rank+1;

    // read xcoord from HDF5
    float ycoord[rank], ycoord_ext[rank+1], dy;
    H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ycoord);
    // extend coordinate
    dy = abs(ycoord[1] - ycoord[0]);
    for (int i = 0; i < rank; i++)
    {
        ycoord_ext[i] = ycoord[i]-dy/2;
    }
    ycoord_ext[rank] = ycoord_ext[rank-1] + dy;

    coords[1] = vtkFloatArray::New();
    coords[1]->SetNumberOfTuples(rank+1);

    for (int i = 0; i < rank+1; i++)
    {
        coords[1]->SetComponent(i, 0, ycoord_ext[i]);
    }
    // close dset
    H5Dclose(dset_id);

    if (dimension > 2)
    {
        // zcoord
        // open dset and dspace
        dset_id = H5Dopen(tmp_file_id, "zCoord", H5P_DEFAULT);
        dspace_id = H5Dget_space(dset_id);

        // count number of nodes on axis
        H5Sget_simple_extent_dims(dspace_id, count, NULL);
        rank = count[0];
        dims[2] = rank+1;

        // read xcoord from HDF5
        float zcoord[rank], zcoord_ext[rank+1], dz;
        H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, zcoord);
        // extend coordinate
        dz = abs(zcoord[1] - zcoord[0]);
        for (int i = 0; i < rank; i++)
        {
            zcoord_ext[i] = zcoord[i]-dz/2;
        }
        zcoord_ext[rank] = zcoord_ext[rank-1] + dz;

        coords[2] = vtkFloatArray::New();
        coords[2]->SetNumberOfTuples(rank+1);

        for (int i = 0; i < rank+1; i++)
        {
            coords[2]->SetComponent(i, 0, zcoord_ext[i]);
        }
        // close dset
        H5Dclose(dset_id);
    }
    else
    {
        dims[2] = 1;
        coords[2] = vtkFloatArray::New();
        coords[2]->SetNumberOfTuples(1);
        coords[2]->SetComponent(0, 0, 0);
    }

    // making grid
    vtkRectilinearGrid  *rGrid = vtkRectilinearGrid::New();
    rGrid->SetDimensions(dims);
    rGrid->SetXCoordinates(coords[0]);
    coords[0]->Delete();
    rGrid->SetYCoordinates(coords[1]);
    coords[1]->Delete();
    rGrid->SetZCoordinates(coords[2]);
    coords[2]->Delete();


    // close temporarily file
    H5Fclose(tmp_file_id);

    return rGrid;
}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::GetVar
//
//  Purpose:
//      Gets a scalar variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      varname    The name of the variable requested.
//
//  Programmer: ylee -- generated by xml2avt
//  Creation:   Fri Feb 7 15:34:32 PST 2020
//
// ****************************************************************************
vtkDataArray *
avtSlugCodeFileFormat::GetVar(const char *varname)
{

    GetSimInfo();

    string vn_str = varname;
    int            pos, ndims, nvar, ntuples;
    hid_t          dset_id, dspace_id, memspace_id;

    // open all primitive variables
    dset_id = H5Dopen(file_id, "prim_vars", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);
    ndims = H5Sget_simple_extent_ndims(dspace_id);

    // count number of zones & variables
    // dims[nx][ny][nz][nvar]
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);
    nvar = dims[ndims-1];

    // count number of zones for each variable
    ntuples = 1;
    for (int i = 0; i < ndims-1; i++)
    {
        ntuples *= dims[i];
    }

    // determine variable name
    if (string(vn_str.c_str()) == "dens")
    {
        pos = 0;
    }
    else if (string(vn_str.c_str()) == "velx")
    {
        pos = 1;
    }
    else if (string(vn_str.c_str()) == "vely")
    {
        pos = 2;
    }
    else if (string(vn_str.c_str()) == "velz")
    {
        pos = 3;
    }
    else if (string(vn_str.c_str()) == "pres")
    {
        if (ndims > 3)  // 3D
            pos = 4;
        else            // 2D
            pos = 3;
    }
    else
    {
        EXCEPTION1(InvalidVariableException, vn_str.c_str());
    }

    // selecting hyperslab
    hsize_t start[ndims], stride[ndims], count[ndims];
    for (int i = 0; i < ndims; i++)
    {
        start[i] = 0;
        stride[i] = 1;
        count[i] = dims[i];
    }
    count[ndims-1] = 1;
    start[ndims-1] = pos;

    // select hyperslab & memspace
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, start, stride, count, NULL);
    memspace_id = H5Screate_simple(ndims, count, NULL);

    // init grid array
    vtkFloatArray * fa = vtkFloatArray::New();
    fa->SetNumberOfTuples(ntuples);
    float *data = fa->GetPointer(0);

    // read data
    H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace_id, dspace_id, H5P_DEFAULT, data);

    // close primitive variables
    H5Dclose(dset_id);


    return fa;
}


// ****************************************************************************
//  Method: avtSlugCodeFileFormat::GetVectorVar
//
//  Purpose:
//      Gets a vector variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      varname    The name of the variable requested.
//
//  Programmer: ylee -- generated by xml2avt
//  Creation:   Fri Feb 7 15:34:32 PST 2020
//
// ****************************************************************************
vtkDataArray *
avtSlugCodeFileFormat::GetVectorVar(const char *varname)
{
    // Not yet implemented
    return 0;
}
