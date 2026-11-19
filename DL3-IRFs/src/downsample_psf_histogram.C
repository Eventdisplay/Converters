// Reduce the statistics of the PSF histogram AngularPSF2DEtrue_offaxis (Poisson resampling)
// for testing the PSF_3GAUSS fits, e.g.:
//   root -l -b -q 'downsample_psf_histogram.C("in.root","out.root",0.1)'
// Usage: downsample_psf_histogram( input file, output file, fraction of events to keep )
void downsample_psf_histogram( const char* in, const char* out, double f )
{
    gSystem->CopyFile( in, out, kTRUE );
    TFile file( out, "UPDATE" );
    TH3F* h = (TH3F*)file.Get( "AngularPSF2DEtrue_offaxis" );
    TRandom3 rnd( 12345 );
    for( int k = 1; k <= h->GetNbinsZ(); k++ )
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
        // weight of a single event: smallest filled bin of the radial profile
        double w = 1.e99;
        for( int j = 1; j <= h->GetNbinsY(); j++ )
        {
            const double c = h->GetBinContent( i, j, k );
            if( c > 0. ) w = TMath::Min( w, c );
        }
        for( int j = 1; j <= h->GetNbinsY(); j++ )
        {
            const double c = h->GetBinContent( i, j, k );
            if( c > 0. )
            {
                const double n = rnd.Poisson( c / w * f );
                h->SetBinContent( i, j, k, n * w / f );
                h->SetBinError( i, j, k, sqrt( n ) * w / f );
            }
        }
    }
    h->Write( "", TObject::kOverwrite );
    file.Close();
}
