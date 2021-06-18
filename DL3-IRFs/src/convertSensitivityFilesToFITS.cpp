/* \file convertSensitivityFilesToFITS.cpp
 *
 * convert a IRF ROOT file into a FITS file
 *
 */

#include "VDL3IRFs.h"
#include "TSystem.h"

using namespace std;

int main( int argc, char* argv[] )
{
    string fFitsFileName = "data/test.fits";
    gSystem->Exec( "rm -v -f data/test.fits" );
    string fRootFile = "./data/DESY.g20210610.V3.ID0NIM3LST3MST3SST4SCMST3.prod5-Paranal-20deg-sq10-LL.S-M6C5ax-14MSTs37SSTs-MSTF.180000s.root";
    TFile *fData = new TFile( fRootFile.c_str() );
    if( fData->IsZombie() )
    {
        cout << "file not found: " << fData->GetName() << endl;
        exit( EXIT_FAILURE );
    }

    VDL3IRFs a;
    a.open_fits_file( fFitsFileName );

    a.write_fits_header();

    // effective area
    cout << "Writing effective area" << endl;
    a.write_effarea(
         (TH2F*)fData->Get( "EffectiveAreaEtrueNoTheta2cut_offaxis" ) );

    // psf
    cout << "Writing background PSF Gauss" << endl;
    a.write_psf_gauss( 
         (TH2F*)fData->Get( "AngResEtrue_offaxis" ) );

    // edisp
    // note: note using migration matrix MigMatrixNoTheta2cut_offaxis
    //       with fine (500) binning on E_true axis
    cout << "Writing energy dispersion matrix" << endl;
    a.write_edisp( 
         (TH3F*)fData->Get( "EestOverEtrueNoTheta2cut_offaxis" ) );

    // background rates
    cout << "Writing background IRFs" << endl;
    a.write_background( 
         (TH2F*)fData->Get( "BGRatePerSqDeg_offaxis" ) );

    a.write_fits_file();

    return 0;
}
