/* \file convertSensitivityFilesToFITS.cpp
 *
 * convert a IRF ROOT file into a FITS file
 *
 */

#include "VDL3IRFs.h"
#include "TSystem.h"

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
    if( argc != 4 )
    {
        cout << endl;
        cout << "./convertSensitivityFilesToFITS <input root file (.root)> ";
        cout << "<output FITS file (.fits/.fits.gz)> ";
        cout << " <2D/3D>" << endl;
        cout << endl;
        cout << " converts IRFs in ROOT format to GADF FITS format" << endl;
        cout << "     existing files are overwritten" << endl;
        cout << endl;
        cout << "\t <2D/3D> 3D: extrapolate 2D histograms to 3D" << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    cout << "Converting IRFs in ROOT to GADF FITS format" << endl;
    cout << "-------------------------------------------" << endl;
    cout << endl;

    string fRootFile = argv[1];
    string fFitsFileName = argv[2];
    string f2D3D = argv[3];

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
             (char*)"CTA (MC prod5, v0.1)",
             (char*)getArrayName(fData) );

    // effective area
    cout << "Writing effective area" << endl;
    a.write_effarea(
         (TH2F*)fData->Get( "EffectiveAreaEtrueNoTheta2cut_offaxis" ) );

    // psf
    cout << "Writing gamma-ray point-spread function (Gaussian approximation)" << endl;
    a.write_psf_gauss( 
         (TH2F*)fData->Get( "AngResEtrue_offaxis" ) );

    /*
    cout << "Writing gamma-ray point-spread function (3D table)" << endl;
    a.write_psf_table(
         (TH3F*)fData->Get( "AngularPSF2DEtrue_offaxis" ) ); 
    */

    // edisp
    // note: note using migration matrix MigMatrixNoTheta2cut_offaxis
    //       with fine (500) binning on E_true axis
    cout << "Writing energy dispersion matrix" << endl;
    a.write_edisp( 
         (TH3F*)fData->Get( "EestOverEtrueNoTheta2cut_offaxis" ) );

    // background rates
    cout << "Writing background IRFs" << endl;
    if( f2D3D == "2D" )
    {
        a.write_background( 
             (TH2F*)fData->Get( "BGRatePerSqDeg_offaxis" ) );
    }
    else
    {
        a.write_background_3D_from_2d(
             (TH2F*)fData->Get( "BGRatePerSqDeg_offaxis" ) );
    }

    a.write_fits_file();

    return 0;
}
