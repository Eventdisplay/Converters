/* convert IRFS from ROOT to FITS
*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <fitsio.h>

#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"

using namespace std;

class VDL3IRFs
{
    private:

    fitsfile* fptr;

    vector< vector< float > > get_baseline_axes( TH1 *h );
    vector< float > calculate_norm_mev_background( TH1 *h );
    void normalise_pdf( TH3F *h );
    bool printerror( int status );
    bool write_fits_keyword( char*, char*, char* );
    bool write_fits_table_header( string irftype, char* instrument );
    bool write_table( vector< vector< float > > table );
    bool write_histo2D( TH2F *h,
                        string name,
                        char* col_name,
                        char* col_unit,
                        bool MEV_BACKGROUND_UNIT = false );


    public:

    VDL3IRFs();
   ~VDL3IRFs() {}

    bool open_fits_file( string fits_file_name );
    bool write_fits_header( char*, char* );
    bool write_background( TH2F *h, char* );
    bool write_background_3D_from_2d( TH2F *h, char* );
    bool write_effarea( TH2F *h, char* );
    bool write_edisp( TH3F *h, char* );
    bool write_diffsens( TH2F*, char* );
    bool write_psf_gauss( TH2F *h, char* );
    bool write_psf_table( TH3F *h, char* );
    bool write_fits_file();
};
