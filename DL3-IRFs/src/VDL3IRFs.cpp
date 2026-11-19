#include "VDL3IRFs.h"
#include <algorithm>
#include <cmath>
#include <utility>

VDL3IRFs::VDL3IRFs()
{
    fptr = 0;
}

/*
 * fits error reporting
*/
bool VDL3IRFs::printerror( int status )
{
    if( status )
    {
        fits_report_error( stderr, status );
    }
    return false;
}

/*
 *  create a new fits file
 */
bool VDL3IRFs::open_fits_file( string fits_file_name )
{
    int status = 0;
    if( fits_create_file( &fptr, fits_file_name.c_str(), &status ) )
    {
        return printerror( status );
    }
    if( fits_create_img( fptr,  DOUBLE_IMG, 0, 0, &status ) )
    {
        return printerror( status );
    }

    return true;
}

bool VDL3IRFs::write_fits_file()
{
    int status = 0;
    if( fits_close_file( fptr, &status ) )
    {
         return printerror( status );
    }
    return true;
}

/*
   write FITS header with
   basic information
*/
bool VDL3IRFs::write_fits_header(
             char *telescope,
             char *instrument )
{
   int status = 0;
   cout << "Instrument: " << instrument << endl;

   char author[] = "G.Maier";
   if( fits_update_key( fptr, TSTRING,
                        (char*)"AUTHOR",
                        author,
                        (char*)"Author",
                        &status ) )
   {
       return printerror( status );
   }
   if( fits_update_key( fptr, TSTRING,
                        (char*)"TELESCOP",
                        telescope,
                        (char*)"Telescope name",
                        &status ) )
   {
       return printerror( status );
   }
   if( fits_update_key( fptr, TSTRING,
                        (char*)"INSTRUME",
                        instrument,
                        (char*)"Instrument",
                        &status ) )
   {
       return printerror( status );
   }
   return true;
}

bool VDL3IRFs::write_fits_keyword(
                    char *key_name,
                    char *key_value,
                    char *key_comment )
{
   int status = 0;
   if( fits_update_key( fptr,
                        TSTRING,
                        key_name,
                        key_value,
                        key_comment,
                        &status ) )
   {
       return printerror( status );
   }
   return true;
}

bool VDL3IRFs::write_fits_table_header( string irftype, char *instrument )
{
   cout << "write_fits_table_header " << irftype << endl;

   write_fits_keyword( (char*)"TELESCOP",
                       (char*)"CTA",
                       (char*)"Name of telescope" );

   write_fits_keyword( (char*)"INSTRUME",
                       instrument,
                       (char*)"Instrument" );

   write_fits_keyword( (char*)"AUTHOR",
                       (char*)"G.Maier",
                       (char*)"Author" );
   write_fits_keyword( (char*)"ORIGIN",
                       (char*)"CTAO",
                       (char*)"Origin" );

   // date string
   char datestr[FLEN_KEYWORD];
   int timeref = 0;
   int status = 0;
   if( fits_get_system_time( datestr,
                             &timeref,
                             &status ) )
   {
       return printerror( status );
   }
   write_fits_keyword( (char*)"DATE",
                       datestr,
                       (char*)"File creation date (YYYY-MM-DDThh:mm:ss UTC)" );

   if( irftype == "DIFFSENS_2D" ) return true;

   write_fits_keyword( (char*)"HDUDOC",
                       (char*)"https://github.com/open-gamma-ray-astro/gamma-astro-data-formats",
                       (char*)"" );

   write_fits_keyword( (char*)"HDUVERS",
                       (char*)"0.3",
                       (char*)"HDU version" );

   write_fits_keyword( (char*)"HDUCLASS",
                       (char*)"GADF",
                       (char*)"HDUCLASS" );

   write_fits_keyword( (char*)"HDUCLAS1",
                       (char*)"RESPONSE",
                       (char*)"HDUCLAS1" );

   if( irftype == "PSF_3GAUSS" || irftype == "PSF_TABLE" )
   {
       write_fits_keyword( (char*)"HDUCLAS2",
                           (char*)"PSF",
                           (char*)"HDUCLAS2" );
   }
   else if( irftype.find( "BKG" ) != string::npos )
   {
       write_fits_keyword( (char*)"HDUCLAS2",
                           (char*)"BKG",
                           (char*)"HDUCLAS2" );
   }
   else if( irftype == "AEFF_2D" )
   {
       write_fits_keyword( (char*)"HDUCLAS2",
                           (char*)"EFF_AREA",
                           (char*)"HDUCLAS2" );
   }
   else if( irftype == "EDISP_2D" )
   {
       write_fits_keyword( (char*)"HDUCLAS2",
                           (char*)"EDISP",
                           (char*)"HDUCLAS2" );
   }

   write_fits_keyword( (char*)"HDUCLAS3",
                       (char*)"FULL-ENCLOSURE",
                       (char*)"HDUCLAS3" );

   write_fits_keyword( (char*)"HDUCLAS4",
                       (char*)irftype.c_str(),
                       (char*)"HDUCLAS4" );

   if( irftype.find( "BKG" ) != string::npos )
   {
       write_fits_keyword( (char*)"FOVALIGN",
//                           (char*)"ALTAZ",
                           (char*)"RADEC",
                           (char*)"FOV alignment" );
   }

   return true;
}

/*
   normalise rows of fixed E_true and off-axis bins
   to 1 (pdf)
*/
void VDL3IRFs::normalise_pdf( TH3F *h )
{
    if( !h ) return;

    // off-axis bins (fix)
    for( int j = 0; j < h->GetNbinsZ(); j++ )
    {
        // true energy axis bins (fix)
        for( int i = 0; i < h->GetNbinsX(); i++ )
        {
            // energy dispersion: normalise to one
            double sum = 0.;
            for( int k = 0; k < h->GetNbinsY(); k++ )
            {
                sum += h->GetBinContent( i+1, k+1, j+1 );
            }
            // require at least then entries
            // to be useful
            if( sum > 10. )
            {
                for( int k = 0; k < h->GetNbinsY(); k++ )
                {
                    h->SetBinContent( i+1, k+1, j+1,
                              h->GetBinContent(  i+1, k+1, j+1 ) / sum );
                }
            }
            // otherwise set all bins to zero
            else
            {
                for( int k = 0; k < h->GetNbinsY(); k++ )
                {
                    h->SetBinContent( i+1, k+1, j+1, 0. );
                }
            }
        }
    }
}

/*
 * Write energy dispersion
 * https://gamma-astro-data-formats.readthedocs.io/en/latest/irfs/full_enclosure/edisp/index.html
 *
*/
bool VDL3IRFs::write_edisp( TH3F *h, char *instrument )
{
   if( !h ) return false;

   // make sure that PDFs are normalised
   // (this is not guaranteed by Eventdisplay)
   normalise_pdf( h );

   int status = 0;

   const int nCol = 7;
   long nRows = h->GetNbinsX()
              * h->GetNbinsZ()
              * h->GetNbinsY();
   nRows = 0;
   char* tType[nCol] = { (char*)"ENERG_LO",
                      (char*)"ENERG_HI",
                      (char*)"MIGRA_LO",
                      (char*)"MIGRA_HI",
                      (char*)"THETA_LO",
                      (char*)"THETA_HI",
                      (char*)"MATRIX"};
   char* tUnit[nCol] = { (char*)"TeV",
                      (char*)"TeV",
                      (char*)"",
                      (char*)"",
                      (char*)"deg",
                      (char*)"deg",
                      (char*)""};
   // e_true axis is x-axis
   char x_form[10];
   sprintf( x_form, "%dE", h->GetNbinsX() );
   // e_mig = e_reco / e_true
   char y_form[10];
   sprintf( y_form, "%dE", h->GetNbinsY() );
   // offset angle axis is z-axis
   char z_form[10];
   sprintf( z_form, "%dE", h->GetNbinsZ() );
   // mig
   char m_form[10];
   sprintf( m_form, "%dE", h->GetNbinsX()
                         * h->GetNbinsZ()
                         * h->GetNbinsY() );
   char* tForm[nCol] = { &x_form[0],
                      &x_form[0],
                      &y_form[0],
                      &y_form[0],
                      &z_form[0],
                      &z_form[0],
                      &m_form[0] };
   ///////////////
   // create empty table
   if( fits_create_tbl( fptr,
                        BINARY_TBL,
                        nRows,
                        nCol,
                        tType,
                        tForm,
                        tUnit,
                        "ENERGY DISPERSION",
                        &status ) )
   {
       return printerror( status );
   }
   // set dimensions
   long int naxes[] = { h->GetNbinsX(),
                        h->GetNbinsY(),
                        h->GetNbinsZ() };
   if( fits_write_tdim( fptr,
                        7,
                        3,
                        naxes,
                        &status ) )
   {
      return printerror( status );
   }

   ///////////////
   // write data
   vector< vector< float > > table = get_baseline_axes( h );

   vector< float > data;
   for( int k = 0; k < h->GetNbinsZ(); k++ )
   {
       for( int j = 0; j < h->GetNbinsY(); j++ )
       {
           double width = h->GetYaxis()->GetBinWidth( j+1 );
           for( int i = 0; i < h->GetNbinsX(); i++ )
           {
               data.push_back( width > 0. ? h->GetBinContent( i+1, j+1, k+1 ) / width : 0. );
           }
       }
   }
   table.push_back( data );

   bool writing_success = write_table( table );
   if( writing_success )
   {
       write_fits_table_header( "EDISP_2D", instrument );
   }
   return writing_success;
}

/*
    PSF table
    https://gamma-astro-data-formats.readthedocs.io/en/latest/irfs/full_enclosure/psf/psf_table/index.html

*/
bool VDL3IRFs::write_psf_table( TH3F *h, char *instrument )
{
   if( !h ) return false;

   int status = 0;
   // make sure that PDfs are normalised
   // (this is not guaranteed by Eventdisplay)
   normalise_pdf( h );

   const int nCol = 7;
   long nRows = h->GetNbinsX()
              * h->GetNbinsZ()
              * h->GetNbinsY();
   nRows = 0;
   char* tType[nCol] = { (char*)"ENERG_LO",
                      (char*)"ENERG_HI",
                      (char*)"THETA_LO",
                      (char*)"THETA_HI",
                      (char*)"RAD_LO",
                      (char*)"RAD_HI",
                      (char*)"RPSF" };
   char* tUnit[nCol] = { (char*)"TeV",
                      (char*)"TeV",
                      (char*)"deg",
                      (char*)"deg",
                      (char*)"deg",
                      (char*)"deg",
                      (char*)"sr^-1" };
   // e_true axis is x-axis
   char x_form[10];
   sprintf( x_form, "%dE", h->GetNbinsX() );
   // Offset angle from source position (y-axis, RAD)
   char y_form[10];
   sprintf( y_form, "%dE", h->GetNbinsY() );
   // Field of view offset angle (z-axis, THETA)
   char z_form[10];
   sprintf( z_form, "%dE", h->GetNbinsZ() );
   // mig
   char m_form[10];
   sprintf( m_form, "%dE", h->GetNbinsX()
                         * h->GetNbinsZ()
                         * h->GetNbinsY() );
   char* tForm[nCol] = { &x_form[0],
                      &x_form[0],
                      &z_form[0],
                      &z_form[0],
                      &y_form[0],
                      &y_form[0],
                      &m_form[0] };

   ///////////////
   // create empty table
   if( fits_create_tbl( fptr,
                        BINARY_TBL,
                        nRows,
                        nCol,
                        tType,
                        tForm,
                        tUnit,
                        "PSF_TABLE",
                        &status ) )
   {
       return printerror( status );
   }
   // set dimensions
   long int naxes[] = { h->GetNbinsX(),
                        h->GetNbinsZ(),
                        h->GetNbinsY() };
   if( fits_write_tdim( fptr,
                        7,
                        3,
                        naxes,
                        &status ) )
   {
      return printerror( status );
   }
   ///////////////
   // write data
   vector< vector< float > > table = get_baseline_axes( h );
   std::swap( table[2], table[4] );
   std::swap( table[3], table[5] );

   vector< float > data;
   for( int j = 0; j < h->GetNbinsY(); j++ )
   {
       double radius_low = h->GetYaxis()->GetBinLowEdge( j+1 ) * TMath::DegToRad();
       double radius_high = h->GetYaxis()->GetBinUpEdge( j+1 ) * TMath::DegToRad();
       double solid_angle = 4. * TMath::Pi()
                            * sin( ( radius_high + radius_low ) / 2. )
                            * sin( ( radius_high - radius_low ) / 2. );
       for( int k = 0; k < h->GetNbinsZ(); k++ )
       {
           for( int i = 0; i < h->GetNbinsX(); i++ )
           {
               data.push_back( solid_angle > 0. ? h->GetBinContent( i+1, j+1, k+1 ) / solid_angle : 0. );
           }
       }
   }
   table.push_back( data );

   bool writing_success = write_table( table );
   if( writing_success )
   {
       write_fits_table_header( "PSF_TABLE", instrument );
   }
   return writing_success;
}

namespace
{

const int N_GAUSS = 3;

struct GaussMixture
{
    // component weights (sum to one), sigma in deg (ascending)
    double weight[N_GAUSS];
    double sigma[N_GAUSS];
    bool valid;
};

/*
 * Fit a mixture of 2D Gaussians to a binned radial distribution
 *
 * Binned maximum-likelihood fit using the EM algorithm, with the radial bin
 * probabilities of a 2D Gaussian (Rayleigh distribution) integrated
 * over each bin: exp(-r_lo^2/2s^2) - exp(-r_hi^2/2s^2).
 *
 * r_lo, r_hi:  bin edges (deg)
 * prob:        bin content (any normalisation)
 * r68:         68% containment radius, used for the starting values
 * sigma_ratio: null for free sigmas; otherwise sigma_k / sigma_1 of the
 *              components, fixing the shape of the mixture so that only the
 *              scale (sigma_1) and the weights are fitted
 * n_eff_out:   effective number of entries of the distribution
 */
GaussMixture fit_gauss_mixture( const vector< double >& r_lo,
                                const vector< double >& r_hi,
                                const vector< double >& prob,
                                double r68,
                                const double* sigma_ratio,
                                double* n_eff_out )
{
    GaussMixture m;
    m.valid = false;
    for( int k = 0; k < N_GAUSS; k++ )
    {
        m.weight[k] = 0.;
        m.sigma[k] = 0.;
    }

    // only filled bins contribute to the multinomial likelihood
    vector< double > a, b, w;
    double sum = 0.;
    double sum2 = 0.;
    double min_width = 1.e99;
    for( unsigned int i = 0; i < prob.size(); i++ )
    {
        if( prob[i] > 0. )
        {
            a.push_back( r_lo[i] );
            b.push_back( r_hi[i] );
            w.push_back( prob[i] );
            sum += prob[i];
            sum2 += prob[i] * prob[i];
            min_width = min( min_width, r_hi[i] - r_lo[i] );
        }
    }
    if( w.empty() || sum <= 0. || r68 <= 0. ) return m;
    // effective number of entries (bin contents can be weighted)
    const double n_eff = sum * sum / sum2;
    if( n_eff_out ) *n_eff_out = n_eff;
    for( unsigned int i = 0; i < w.size(); i++ ) w[i] /= sum;

    // sigma below 1/2 bin width cannot be constrained by the data
    const double sigma_min = 0.5 * min_width;
    const double sigma_max = 10. * ( b.back() );

    const unsigned int n = w.size();
    vector< double > gamma_et[N_GAUSS];
    vector< double > gamma[N_GAUSS];
    for( int k = 0; k < N_GAUSS; k++ )
    {
        gamma_et[k].assign( n, 0. );
        gamma[k].assign( n, 0. );
    }

    // EM iterations from the given starting values (updated in place),
    // returns the log likelihood
    auto run_em = [&]( double* sigma, double* weight ) -> double
    {
    double last_ll = -1.e99;
    for( int iter = 0; iter < 20000; iter++ )
    {
        // E step
        double ll = 0.;
        for( unsigned int i = 0; i < n; i++ )
        {
            double q[N_GAUSS];
            double et[N_GAUSS];
            double total = 0.;
            for( int k = 0; k < N_GAUSS; k++ )
            {
                const double s2 = 2. * sigma[k] * sigma[k];
                const double ta = a[i] * a[i] / s2;
                const double tb = b[i] * b[i] / s2;
                const double ea = exp( -ta );
                const double eb = exp( -tb );
                const double p = ea - eb;
                if( p > 1.e-300 )
                {
                    // mean of t = r^2/2sigma^2 for an exponential truncated to [ta,tb]
                    et[k] = ( ( ta + 1. ) * ea - ( tb + 1. ) * eb ) / p;
                }
                else
                {
                    et[k] = 0.5 * ( ta + tb );
                }
                q[k] = weight[k] * p;
                total += q[k];
            }
            if( total > 1.e-300 )
            {
                ll += w[i] * log( total );
                for( int k = 0; k < N_GAUSS; k++ )
                {
                    gamma[k][i] = w[i] * q[k] / total;
                    gamma_et[k][i] = gamma[k][i] * et[k];
                }
            }
            else
            {
                // bin far outside all components: assign to the widest one
                for( int k = 0; k < N_GAUSS; k++ )
                {
                    gamma[k][i] = ( k == N_GAUSS - 1 ) ? w[i] : 0.;
                    gamma_et[k][i] = gamma[k][i] * et[k];
                }
            }
        }

        // M step
        double g_all = 0.;
        double ge_all = 0.;
        for( int k = 0; k < N_GAUSS; k++ )
        {
            double g = 0.;
            double ge = 0.;
            for( unsigned int i = 0; i < n; i++ )
            {
                g += gamma[k][i];
                ge += gamma_et[k][i];
            }
            weight[k] = g;
            g_all += g;
            ge_all += ge;
            if( !sigma_ratio && g > 1.e-12 )
            {
                sigma[k] *= sqrt( ge / g );
                sigma[k] = min( max( sigma[k], sigma_min ), sigma_max );
            }
        }
        if( sigma_ratio && g_all > 0. )
        {
            // common scale of the components with fixed ratios
            const double s = min( max( sigma[0] * sqrt( ge_all / g_all ), sigma_min ), sigma_max );
            for( int k = 0; k < N_GAUSS; k++ ) sigma[k] = s * sigma_ratio[k];
        }

        if( iter > 0 && fabs( ll - last_ll ) < 1.e-13 * fabs( ll ) ) break;
        last_ll = ll;
    }
    return last_ll;
    };

    // the likelihood of a mixture has local maxima: start from Gaussians
    // scaled from the Gaussian-equivalent of R68 and keep the best result
    const double sigma_r68 = r68 * 0.6624305;
    const double start_scale[][N_GAUSS] = { { 0.5, 1.0, 2.5 },
                                            { 0.3, 0.8, 1.5 },
                                            { 0.6, 1.5, 4.0 } };
    const double start_weight[N_GAUSS] = { 0.5, 0.35, 0.15 };
    double sigma[N_GAUSS] = { 0., 0., 0. };
    double weight[N_GAUSS] = { 0., 0., 0. };
    double best_ll = -1.e99;
    for( const auto& scales : start_scale )
    {
        double s[N_GAUSS];
        double wt[N_GAUSS];
        for( int k = 0; k < N_GAUSS; k++ )
        {
            // fixed shape: components scaled to start with a Gaussian of R68
            s[k] = sigma_ratio ? sigma_r68 * scales[1] * sigma_ratio[k] / sigma_ratio[1]
                               : sigma_r68 * scales[k];
            s[k] = min( max( s[k], sigma_min ), sigma_max );
            wt[k] = start_weight[k];
        }
        const double ll = run_em( s, wt );
        if( ll > best_ll )
        {
            best_ll = ll;
            for( int k = 0; k < N_GAUSS; k++ )
            {
                sigma[k] = s[k];
                weight[k] = wt[k];
            }
        }
    }

    // order by increasing sigma
    int order[N_GAUSS] = { 0, 1, 2 };
    sort( order, order + N_GAUSS,
          [&]( int x, int y ) { return sigma[x] < sigma[y]; } );
    for( int k = 0; k < N_GAUSS; k++ )
    {
        m.sigma[k] = sigma[order[k]];
        m.weight[k] = weight[order[k]];
    }
    m.valid = ( m.weight[0] > 0. );
    return m;
}

/*
 * containment radius of a binned radial distribution (linear interpolation)
 */
double containment_radius( const vector< double >& r_hi,
                           const vector< double >& prob,
                           double containment )
{
    double sum = 0.;
    for( unsigned int i = 0; i < prob.size(); i++ ) sum += prob[i];
    if( sum <= 0. ) return 0.;
    double cum = 0.;
    double r_last = 0.;
    for( unsigned int i = 0; i < prob.size(); i++ )
    {
        const double next = cum + prob[i] / sum;
        if( next >= containment )
        {
            const double f = prob[i] > 0. ? ( containment - cum ) / ( next - cum ) : 0.;
            return r_last + f * ( r_hi[i] - r_last );
        }
        cum = next;
        r_last = r_hi[i];
    }
    return r_last;
}

}

/*
 *  Multi-Gauss mixture model
 *  https://gamma-astro-data-formats.readthedocs.io/en/latest/irfs/full_enclosure/psf/psf_3gauss/index.html#
 *
 *  Fits a mixture of three 2D Gaussians to the radial PSF distribution in each
 *  true-energy and off-axis bin of the PSF histogram (same input as PSF_TABLE):
 *
 *  P(r) = SCALE * ( exp(-r^2/2 SIGMA_1^2) + AMPL_2 exp(-r^2/2 SIGMA_2^2)
 *                                         + AMPL_3 exp(-r^2/2 SIGMA_3^2) )
 *
 *  with SCALE = 1 / ( 2 pi ( SIGMA_1^2 + AMPL_2 SIGMA_2^2 + AMPL_3 SIGMA_3^2 ) )
 *  in sr^-1 (sigma in rad), such that P integrates to one.
 *
 *  h: PSF histogram (x: log10 E_true, y: radius [deg], z: off-axis angle [deg])
*/
bool VDL3IRFs::write_psf_gauss( TH3F *h, char *instrument )
{
   if( !h ) return false;

   int status = 0;
   const int nE = h->GetNbinsX();
   const int nT = h->GetNbinsZ();
   const int nR = h->GetNbinsY();
   const int nCol = 10;
   long nRows = 0;
   char* tType[nCol] = { (char*)"ENERG_LO",
                      (char*)"ENERG_HI",
                      (char*)"THETA_LO",
                      (char*)"THETA_HI",
                      (char*)"SCALE",
                      (char*)"SIGMA_1",
                      (char*)"SIGMA_2",
                      (char*)"SIGMA_3",
                      (char*)"AMPL_2",
                      (char*)"AMPL_3" };
   char* tUnit[nCol] = { (char*)"TeV",
                      (char*)"TeV",
                      (char*)"deg",
                      (char*)"deg",
                      (char*)"sr^(-1)",
                      (char*)"deg",
                      (char*)"deg",
                      (char*)"deg",
                      (char*)"",
                      (char*)"" };
   char x_form[10];
   sprintf( x_form, "%dE", nE );
   char y_form[10];
   sprintf( y_form, "%dE", nT );
   char z_form[10];
   sprintf( z_form, "%dE", nE * nT );
   char* tForm[nCol] = { &x_form[0],
                      &x_form[0],
                      &y_form[0],
                      &y_form[0],
                      &z_form[0],
                      &z_form[0],
                      &z_form[0],
                      &z_form[0],
                      &z_form[0],
                      &z_form[0] };
   ///////////////
   // create empty table
   if( fits_create_tbl( fptr,
                        BINARY_TBL,
                        nRows,
                        nCol,
                        tType,
                        tForm,
                        tUnit,
                        "POINT SPREAD FUNCTION",
                        &status ) )
   {
       return printerror( status );
   }
   // set dimensions
   long int naxes[] = { nE, nT };
   int dim_colnum[] = { 5, 6, 7, 8, 9, 10 };
   for( int i = 0; i < 6; i++ )
   {
       if( fits_write_tdim( fptr,
                            dim_colnum[i],
                            2,
                            naxes,
                            &status ) )
      {
          return printerror( status );
      }
   }
   ///////////////
   // write data
   vector< vector< float > > table;
   vector< float > xedge_low, xedge_hig, zedge_low, zedge_hig;
   for( int i = 0; i < nE; i++ )
   {
       xedge_low.push_back( TMath::Power( 10., h->GetXaxis()->GetBinLowEdge( i+1 ) ) );
       xedge_hig.push_back( TMath::Power( 10., h->GetXaxis()->GetBinUpEdge( i+1 ) ) );
   }
   for( int k = 0; k < nT; k++ )
   {
       zedge_low.push_back( h->GetZaxis()->GetBinLowEdge( k+1 ) );
       zedge_hig.push_back( h->GetZaxis()->GetBinUpEdge( k+1 ) );
   }
   table.push_back( xedge_low );
   table.push_back( xedge_hig );
   table.push_back( zedge_low );
   table.push_back( zedge_hig );

   vector< double > r_lo, r_hi;
   for( int j = 0; j < nR; j++ )
   {
       r_lo.push_back( h->GetYaxis()->GetBinLowEdge( j+1 ) );
       r_hi.push_back( h->GetYaxis()->GetBinUpEdge( j+1 ) );
   }

   // read radial distributions (same statistics requirement as for PSF_TABLE)
   vector< vector< double > > profile( nE * nT );
   for( int k = 0; k < nT; k++ )
   {
       for( int i = 0; i < nE; i++ )
       {
           vector< double > prob( nR, 0. );
           double sum = 0.;
           for( int j = 0; j < nR; j++ )
           {
               prob[j] = h->GetBinContent( i+1, j+1, k+1 );
               sum += prob[j];
           }
           if( sum > 10. ) profile[k*nE+i] = prob;
       }
   }

   // Fitting every bin independently gives mixtures with unrelated components
   // for noisy bins. As gammapy interpolates the parameters between bins,
   // this results in wrong PSFs at energies between two bins.
   // The shape (ratios of the sigmas) of the mixture is therefore determined
   // from free fits, smoothed over neighbouring energy bins, and kept fixed
   // when fitting the scale and the weights of each bin.
   double n_eff_max = 0.;
   vector< double > n_eff( nE * nT, 0. );
   vector< GaussMixture > free_fit( nE * nT );
   for( unsigned int p = 0; p < profile.size(); p++ )
   {
       free_fit[p].valid = false;
       if( profile[p].empty() ) continue;
       free_fit[p] = fit_gauss_mixture( r_lo, r_hi, profile[p],
                                        containment_radius( r_hi, profile[p], 0.68 ),
                                        0, &n_eff[p] );
       n_eff_max = max( n_eff_max, n_eff[p] );
   }
   // log of sigma_k/sigma_1 and weight of each bin in the smoothing
   vector< double > log_ratio[N_GAUSS];
   vector< double > ratio_weight( nE * nT, 0. );
   for( int c = 0; c < N_GAUSS; c++ ) log_ratio[c].assign( nE * nT, 0. );
   for( unsigned int p = 0; p < profile.size(); p++ )
   {
       const GaussMixture& f = free_fit[p];
       // ignore bins where a component is not used
       if( f.valid && f.weight[1] > 0.02 && f.weight[2] > 0.005 )
       {
           ratio_weight[p] = min( n_eff[p], 0.1 * n_eff_max );
           for( int c = 0; c < N_GAUSS; c++ ) log_ratio[c][p] = log( f.sigma[c] / f.sigma[0] );
       }
   }
   const double kernel_width = 1.0;
   int n_fitted = 0;
   int n_failed = 0;
   vector< float > scale, sigma[N_GAUSS], ampl[N_GAUSS];
   // off-axis angle bins outer, true energy bins inner (as in the FITS table)
   for( int k = 0; k < nT; k++ )
   {
       for( int i = 0; i < nE; i++ )
       {
           const vector< double >& prob = profile[k*nE+i];
           GaussMixture m;
           m.valid = false;
           if( !prob.empty() )
           {
               double sigma_ratio[N_GAUSS] = { 1., 2., 4. };
               double w_sum = 0.;
               double l_sum[N_GAUSS] = { 0., 0., 0. };
               for( int j = 0; j < nE; j++ )
               {
                   const double wj = ratio_weight[k*nE+j]
                                   * exp( -0.5 * ( i - j ) * ( i - j ) / ( kernel_width * kernel_width ) );
                   w_sum += wj;
                   for( int c = 0; c < N_GAUSS; c++ ) l_sum[c] += wj * log_ratio[c][k*nE+j];
               }
               if( w_sum > 0. )
               {
                   for( int c = 0; c < N_GAUSS; c++ ) sigma_ratio[c] = exp( l_sum[c] / w_sum );
               }
               m = fit_gauss_mixture( r_lo, r_hi, prob,
                                      containment_radius( r_hi, prob, 0.68 ),
                                      sigma_ratio, 0 );
               if( m.valid ) n_fitted++;
               else          n_failed++;
           }
           if( m.valid )
           {
               const double s1 = m.sigma[0] * TMath::DegToRad();
               // component k: weight_k / ( 2 pi sigma_k^2 ) = SCALE * AMPL_k
               const double c1 = m.weight[0] / ( 2. * TMath::Pi() * s1 * s1 );
               scale.push_back( c1 );
               sigma[0].push_back( m.sigma[0] );
               ampl[0].push_back( 1. );
               for( int c = 1; c < N_GAUSS; c++ )
               {
                   const double sc = m.sigma[c] * TMath::DegToRad();
                   const double ck = m.weight[c] / ( 2. * TMath::Pi() * sc * sc );
                   sigma[c].push_back( m.sigma[c] );
                   ampl[c].push_back( ck / c1 );
               }
           }
           else
           {
               scale.push_back( 0. );
               for( int c = 0; c < N_GAUSS; c++ ) sigma[c].push_back( 0. );
               ampl[0].push_back( 0. );
               for( int c = 1; c < N_GAUSS; c++ ) ampl[c].push_back( 0. );
           }
       }
   }
   cout << "\tPSF_3GAUSS: fitted " << n_fitted << " of " << nE * nT
        << " (E_true, off-axis) bins";
   if( n_failed > 0 ) cout << ", " << n_failed << " failed";
   cout << endl;

   table.push_back( scale );
   for( int c = 0; c < N_GAUSS; c++ ) table.push_back( sigma[c] );
   // AMPL_1 is fixed to one and not stored
   table.push_back( ampl[1] );
   table.push_back( ampl[2] );

   bool writing_success = write_table( table );
   if( writing_success )
   {
       write_fits_table_header( "PSF_3GAUSS", instrument );
   }
   return writing_success;
}

/*
    normalisation factor from TeV to MeV conversion

*/
vector< float > VDL3IRFs::calculate_norm_mev_background( TH1 *h )
{
   if( !h )
   {
        vector< float > norm_mev_background;
        return norm_mev_background;
   }
   vector< float > norm_mev_background( h->GetNbinsX(), 1. );
   float dE = 1.;
   for( unsigned int i = 0; i < norm_mev_background.size(); i++ )
   {
      norm_mev_background[i] /= TMath::DegToRad() * TMath::DegToRad();
      // first two axes are E low and high
      // TeV --> MeV
      dE  = TMath::Power( 10., h->GetXaxis()->GetBinUpEdge( i+1 ) );
      dE -= TMath::Power( 10., h->GetXaxis()->GetBinLowEdge( i+1 ) );
      dE *= 1.e6;
      norm_mev_background[i] /= dE;
   }
   return norm_mev_background;
}

// data

/*
 * background IRF
 *
 * expand 2D background IRF to 3D
 *
 * - https://github.com/open-gamma-ray-astro/gamma-astro-data-formats/issues/153
 */
bool VDL3IRFs::write_background_3D_from_2d( TH2F* h, char *instrument )
{
   int status = 0;

   // oversampling multiplier
   int oversample_mult = 5;

   const int nCol = 7;
   long nRows = h->GetNbinsX()
              * h->GetNbinsY() * oversample_mult
              * h->GetNbinsY() * oversample_mult;;
   nRows = 0;
   char* tType[nCol] = { (char*)"ENERG_LO",
                      (char*)"ENERG_HI",
                      (char*)"DETX_LO",
                      (char*)"DETX_HI",
                      (char*)"DETY_LO",
                      (char*)"DETY_HI",
                      (char*)"BKG"};

   char* tUnit[nCol] = { (char*)"TeV",
                      (char*)"TeV",
                      (char*)"deg",
                      (char*)"deg",
                      (char*)"deg",
                      (char*)"deg",
                      (char*)"s^-1 MeV^-1 sr^-1"};
   // e_true axis is x-axis
   char x_form[10];
   sprintf( x_form, "%dE", h->GetNbinsX() );
   // FOV coordinates (4x)
   char y_form[10];
   sprintf( y_form, "%dE", h->GetNbinsY() * oversample_mult * 2 );
   //
   char z_form[10];
   sprintf( z_form, "%dE", h->GetNbinsX()
                         * h->GetNbinsY() * oversample_mult * 2
                         * h->GetNbinsY() * oversample_mult * 2);
   char* tForm[nCol] = { &x_form[0],
                      &x_form[0],
                      &y_form[0],
                      &y_form[0],
                      &y_form[0],
                      &y_form[0],
                      &z_form[0] };
   ///////////////
   // create empty table
   if( fits_create_tbl( fptr,
                        BINARY_TBL,
                        nRows,
                        nCol,
                        tType,
                        tForm,
                        tUnit,
                        "BACKGROUND",
                        &status ) )
   {
       return printerror( status );
   }
   // set dimensions
   long int naxes[] = { h->GetNbinsX(),
                        h->GetNbinsY() * oversample_mult * 2,
                        h->GetNbinsY() * oversample_mult * 2 };
   if( fits_write_tdim( fptr,
                        7,
                        3,
                        naxes,
                        &status ) )
   {
      return printerror( status );
   }

   ///////////////
   // write data
   vector< vector< float > > table;

   // xaxis
   // (always expected to be energy axis in log10)
   vector< float > xedge_low;
   vector< float > xedge_hig;
   for( int i = 0; i < h->GetNbinsX(); i++ )
   {
       xedge_low.push_back( TMath::Power( 10.,
                                         h->GetXaxis()->GetBinLowEdge( i+1 ) ) );
       xedge_hig.push_back( TMath::Power( 10.,
                                         h->GetXaxis()->GetBinUpEdge( i+1 ) ) );
   }
   table.push_back( xedge_low );
   table.push_back( xedge_hig );
   // FOV axis
   vector< float > fov_low;
   vector< float > fov_hig;
   float d_xy = h->GetYaxis()->GetXmax() / (h->GetNbinsY() * oversample_mult);
   float y_min = -1.*h->GetYaxis()->GetXmax();
   for( int i = 0; i < 2.*h->GetNbinsY() * oversample_mult; i++ )
   {
       fov_low.push_back( y_min + i * d_xy );
       fov_hig.push_back( y_min + (i+1) * d_xy );
   }
   table.push_back( fov_low );
   table.push_back( fov_hig );
   table.push_back( fov_low );
   table.push_back( fov_hig );

   vector< float > data;

   vector< float > norm_mev_background = calculate_norm_mev_background( h );

   float dr = 0.;
   for( int k = 0; k < fov_low.size(); k++ )
   {
       for( int j = 0; j < fov_low.size(); j++ )
       {
           dr = sqrt( (fov_low[k]+0.5*d_xy)*(fov_low[k]+0.5*d_xy)
                     +(fov_low[j]+0.5*d_xy)*(fov_low[j]+0.5*d_xy) );

           for( int i = 0; i < h->GetNbinsX(); i++ )
           {
               if( dr < h->GetYaxis()->GetXmax() )
               {
                   data.push_back(
                        h->Interpolate(
                            h->GetXaxis()->GetBinCenter( i+1 ),
                            dr ) );
               }
               else
               {
                   data.push_back( 0. );
               }
               if( i < norm_mev_background.size() )
               {
                   data.back() *= norm_mev_background[i];
               }
           }
       }
   }
   table.push_back( data );

   bool writing_success = write_table( table );

   if( writing_success )
   {
       write_fits_table_header( "BKG_3D", instrument );
   }
   return writing_success;
}

/*
 * background IRF
 * - https://github.com/open-gamma-ray-astro/gamma-astro-data-formats/issues/153
 */
bool VDL3IRFs::write_background( TH2F *h, char *instrument )
{
   bool writing_success = write_histo2D( h,
                      "BACKGROUND",
                      (char*)"BKG",
                      (char*)"s^-1 MeV^-1 sr^-1",
                      true );
   if( writing_success )
   {
       write_fits_table_header( "BKG_2D", instrument );
   }
   return writing_success;
}

/*
 * effective area IRF
 * - https://gamma-astro-data-formats.readthedocs.io/en/latest/irfs/full_enclosure/aeff/index.html
 */
bool VDL3IRFs::write_effarea( TH2F *h, char *instrument )
{
   bool writing_success = write_histo2D( h,
                      "EFFECTIVE AREA",
                      (char*)"EFFAREA",
                      (char*)"m**2",
                      false,
                      true );
   if( writing_success )
   {
       write_fits_table_header( "AEFF_2D", instrument );
   }
   return writing_success;
}

/*
 * differential sensitivity
 *
 */
bool VDL3IRFs::write_diffsens( TH2F *h, char *instrument )
{
   bool writing_success = write_histo2D( h,
                      "DIFFERENTIAL SENSITIVITY",
                      (char*)"DIFFSENS",
                      (char*)"erg cm^-2 s^-1",
                      false );
   if( writing_success )
   {
       write_fits_table_header( "DIFFSENS_2D", instrument );
   }
   return writing_success;
}

/*
   write a 2D histogram to FITS
*/
bool VDL3IRFs::write_histo2D( TH2F *h,
                               string name,
                               char* col_name,
                               char* col_unit,
                               bool MEV_BACKGROUND_UNIT,
                               bool include_uncertainty )
{
   if( !h ) return false;

   int status = 0;
   const int nCol = include_uncertainty ? 6 : 5;
   const int maxCol = 6;  // Maximum possible columns (with uncertainty)
   long nRows = h->GetNbinsX() * h->GetNbinsY();
   nRows = 0;
   
   // Column names - create error column name using std::string
   string err_col_name;
   if( include_uncertainty )
   {
       err_col_name = string(col_name) + "_ERR";
   }
   
   // Arrays sized for maximum columns; fits_create_tbl uses only first nCol entries
   char* tType[maxCol] = { (char*)"ENERG_LO",
                           (char*)"ENERG_HI",
                           (char*)"THETA_LO",
                           (char*)"THETA_HI",
                           (char*)col_name,
                           include_uncertainty ? (char*)err_col_name.c_str() : (char*)"" };
   
   char* tUnit[maxCol] = { (char*)"TeV",
                           (char*)"TeV",
                           (char*)"deg",
                           (char*)"deg",
                           (char*)col_unit,
                           (char*)col_unit };
   char x_form[10];
   sprintf( x_form, "%dE", h->GetNbinsX() );
   char y_form[10];
   sprintf( y_form, "%dE", h->GetNbinsY() );
   char z_form[10];
   sprintf( z_form, "%dE", h->GetNbinsX()*h->GetNbinsY() );
   char* tForm[maxCol] = { &x_form[0],
                           &x_form[0],
                           &y_form[0],
                           &y_form[0],
                           &z_form[0],
                           &z_form[0] };

   ///////////////
   // create empty table
   if( fits_create_tbl( fptr,
                        BINARY_TBL,
                        nRows,
                        nCol,
                        tType,
                        tForm,
                        tUnit,
                        name.c_str() ,
                        &status ) )
   {
       return printerror( status );
   }
   // set dimensions
   long int naxes[] = { h->GetNbinsX(),  h->GetNbinsY() };
   if( fits_write_tdim( fptr,
                        5,
                        2,
                        naxes,
                        &status ) )
   {
      return printerror( status );
   }
   // set dimensions for uncertainty column if present
   if( include_uncertainty )
   {
       if( fits_write_tdim( fptr,
                            6,
                            2,
                            naxes,
                            &status ) )
       {
          return printerror( status );
       }
   }
   ///////////////
   // write data
   vector< vector< float > > table = get_baseline_axes( h );
   // coordinate conversion from
   // 1/deg^2/s --> 1/sr^2/s
   vector< float > norm_mev_background( h->GetNbinsX(), 1. );
   if( MEV_BACKGROUND_UNIT )
   {
      norm_mev_background = calculate_norm_mev_background( h );
   }

   // data
   vector< float > data;
   vector< float > errors;
   for( int j = 0; j < h->GetNbinsY(); j++ )
   {
      for( int i = 0; i < h->GetNbinsX(); i++ )
      {
          // background normalisation
          // (all vector sizes checked)
          if( MEV_BACKGROUND_UNIT )
          {
              data.push_back( h->GetBinContent( i+1, j+1 )
                             * norm_mev_background[i] );
              if( include_uncertainty )
              {
                  // Get the symmetric bin error and scale it with the background normalisation
                  float error = h->GetBinError( i+1, j+1 );
                  errors.push_back( error * norm_mev_background[i] );
              }
          }
          else
          {
              data.push_back( h->GetBinContent( i+1, j+1 ) );
              if( include_uncertainty )
              {
                  errors.push_back( h->GetBinError( i+1, j+1 ) );
              }
          }

      }
   }
   table.push_back( data );
   
   // Add uncertainty column if requested
   if( include_uncertainty )
   {
       table.push_back( errors );
   }

   return write_table( table );
}

bool VDL3IRFs::write_table( vector< vector< float > > table )
{
   int status = 0;

   for( unsigned int i = 0; i < table.size(); i++ )
   {
       cout << "\t writing column " << i+1 << "\t" << table[i].size() << endl;
       if( fits_write_col( fptr,
                       TFLOAT,
                       i+1,
                       1,
                       1,
                       table[i].size(),
                       &table[i][0],
                       &status ) )
       {
           return printerror( status );
       }
   }

   return true;
}

/*
  get baseline axes with
  upper or lower bin edge of a histogram axis
*/
vector<vector<float>> VDL3IRFs::get_baseline_axes( TH1 *h )
{
   vector< vector< float > > table;
   if( !h ) return table;

   // xaxis
   // (always expected to be energy axis in log10)
   vector< float > xedge_low;
   vector< float > xedge_hig;
   for( int i = 0; i < h->GetNbinsX(); i++ )
   {
       xedge_low.push_back( TMath::Power( 10.,
                                         h->GetXaxis()->GetBinLowEdge( i+1 ) ) );
       xedge_hig.push_back( TMath::Power( 10.,
                                         h->GetXaxis()->GetBinUpEdge( i+1 ) ) );
   }
   table.push_back( xedge_low );
   table.push_back( xedge_hig);
   // yaxis
   vector< float > yedge_low;
   vector< float > yedge_hig;
   for( int i = 0; i < h->GetNbinsY(); i++ )
   {
       yedge_low.push_back( h->GetYaxis()->GetBinLowEdge( i+1 ) );
       yedge_hig.push_back( h->GetYaxis()->GetBinUpEdge( i+1 ) );
   }
   table.push_back( yedge_low );
   table.push_back( yedge_hig);

   if( h->GetDimension() == 3 )
   {
       // zaxis
       vector< float > zedge_low;
       vector< float > zedge_hig;
       for( int i = 0; i < h->GetNbinsZ(); i++ )
       {
           zedge_low.push_back( h->GetZaxis()->GetBinLowEdge( i+1 ) );
           zedge_hig.push_back( h->GetZaxis()->GetBinUpEdge( i+1 ) );
       }
       table.push_back( zedge_low );
       table.push_back( zedge_hig);
   }

   return table;
}
