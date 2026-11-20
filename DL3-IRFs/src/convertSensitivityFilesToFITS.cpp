/* \file convertSensitivityFilesToFITS.cpp
 *
 * convert a IRF ROOT file into a FITS file
 *
 */

#include "VDL3IRFs.h"
#include "TSystem.h"
#include <filesystem>

using namespace std;

/*
   try and guess the array name
*/
char* getArrayName( TFile *fData )
{
    if( !fData ) return (char*)"Unknown";

    string file_name = fData->GetName();
    if( file_name.find( "Paranal" ) != string::npos
       || file_name.find( "South" ) != string::npos )
    {
        return (char*)"CTAO Southern Array";
    }

    return (char*)"CTAO Northern Array";
}


int main( int argc, char* argv[] )
{
    if( argc == 1 || ( argc == 2 && string( argv[1] ) == "--help" ) || argc < 4 )
    {
        cout << endl;
        cout << "./convertSensitivityFilesToFITS <input root file (.root)> ";
        cout << "<output FITS file (.fits/.fits.gz)> ";
        cout << " <2D/3D> [noplots | --fit-output-dir DIR]" << endl;
        cout << endl;
        cout << " converts IRFs in ROOT format to GADF FITS format" << endl;
        cout << "     existing files are overwritten" << endl;
        cout << endl;
        cout << "\t <2D/3D> 3D: extrapolate 2D histograms to 3D (default)" << endl;
        cout << "\t noplots: do not write the PSF fit report (.csv) and plots (.png)" << endl;
        cout << "\t --fit-output-dir DIR: write PSF diagnostics to DIR" << endl;
        cout << "\t     default: <FITS-directory>_psf_fit_diagnostics, adjacent to the FITS directory" << endl;
        cout << endl;
        return argc == 1 || ( argc == 2 && string( argv[1] ) == "--help" ) ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    const bool fPSFFitPlots = !( argc == 5 && string( argv[4] ) == "noplots" );
    const bool customFitDir = argc == 6 && string( argv[4] ) == "--fit-output-dir"
                              && argv[5][0] != '\0';
    if( argc != 4 && fPSFFitPlots && !customFitDir )
    {
        cerr << "Expected noplots or --fit-output-dir DIR; see --help" << endl;
        return EXIT_FAILURE;
    }
    cout << "Converting IRFs in ROOT to GADF FITS format" << endl;
    cout << "-------------------------------------------" << endl;
    cout << endl;

    string fRootFile = argv[1];
    string fFitsFileName = argv[2];
    string f2D3D = argv[3];

    string fPSFFitBase;
    if( fPSFFitPlots )
    {
        const auto fitsPath = filesystem::absolute( fFitsFileName ).lexically_normal();
        const auto fitsDir = fitsPath.parent_path();
        const auto reportDir = customFitDir ? filesystem::path( argv[5] )
                               : fitsDir.parent_path() / ( fitsDir.filename().string() + "_psf_fit_diagnostics" );
        auto reportName = fitsPath.filename();
        if( reportName.extension() == ".gz" ) reportName = reportName.stem();
        if( reportName.extension() == ".fits" ) reportName = reportName.stem();

        error_code error;
        filesystem::create_directories( reportDir, error );
        if( error )
        {
            cerr << "Cannot create PSF diagnostics directory " << reportDir << ": " << error.message() << endl;
            return EXIT_FAILURE;
        }
        fPSFFitBase = ( reportDir / reportName ).string();
    }

    cout << "Reading root file from " << fRootFile << endl;
    cout << "Writing FITS file to " << fFitsFileName << endl;
    cout << endl;
    ostringstream o_ss;
    o_ss << "rm -f -v " << fFitsFileName;
    gSystem->Exec( o_ss.str().c_str() );

    TFile *fData = new TFile( fRootFile.c_str() );
    if( fData->IsZombie() )
    {
        cout << "file not found: " << fData->GetName() << endl;
        exit( EXIT_FAILURE );
    }

    VDL3IRFs a;
    if( !a.open_fits_file( fFitsFileName ) )
    {
        exit( EXIT_FAILURE );
    }

    a.write_fits_header(
             (char*)"CTA (MC prod5, v0.1.1)",
             (char*)getArrayName(fData) );

    // effective area
    cout << "Writing effective area" << endl;
    a.write_effarea(
         (TH2F*)fData->Get( "EffectiveAreaEtrueNoTheta2cut_offaxis" ),
         (char*)getArrayName(fData) );

    // psf
    TH3F* psf_table_hist = (TH3F*)fData->Get( "AngularPSF2DEtrue_offaxis" );
    if( psf_table_hist )
    {
        cout << "Writing gamma-ray point-spread function (Gaussian mixture)" << endl;
        if( !a.write_psf_gauss( psf_table_hist, (char*)getArrayName(fData),
                                fPSFFitBase, filesystem::path( fFitsFileName ).filename().string() ) )
        {
            cerr << "Error: Failed to write PSF_3GAUSS" << endl;
        }

        cout << "Writing gamma-ray point-spread function (3D table)" << endl;
        if( !a.write_psf_table( psf_table_hist, (char*)getArrayName(fData) ) )
        {
            cerr << "Error: Failed to write PSF_TABLE" << endl;
        }
    }
    else
    {
        cout << "Warning: AngularPSF2DEtrue_offaxis histogram not found, skipping PSF_3GAUSS and PSF_TABLE" << endl;
    }

    // edisp
    // note: note using migration matrix MigMatrixNoTheta2cut_offaxis
    //       with fine (500) binning on E_true axis
    cout << "Writing energy dispersion matrix" << endl;
    a.write_edisp(
         (TH3F*)fData->Get( "EestOverEtrueNoTheta2cut_offaxis" ),
         (char*)getArrayName(fData) );

    // background rates
    cout << "Writing background IRFs" << endl;
    if( f2D3D == "2D" )
    {
        a.write_background(
             (TH2F*)fData->Get( "BGRatePerSqDeg_offaxis" ),
             (char*)getArrayName(fData) );
    }
    else
    {
        a.write_background_3D_from_2d(
             (TH2F*)fData->Get( "BGRatePerSqDeg_offaxis" ),
             (char*)getArrayName(fData) );
    }

    // differential sensitivity
    cout << "Writing differential sensitivity IRFs" << endl;
    a.write_diffsens(
         (TH2F*)fData->Get( "DiffSens_offaxis" ),
         (char*)getArrayName(fData) );

    a.write_fits_file();

    return 0;
}
