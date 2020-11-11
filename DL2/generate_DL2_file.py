# Script to convert Eventdisplay DL2a output to FITS
#
# - expected Eventdisplay output including gamma/hadron cuts
#
#  Script by T.Hassan
#
import click
import logging

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
cuts_info = '''Cut level of the events to be included.
0: Events before applying gamma/hadron and direction cuts.
1: Events surviving gamma/hadron separation cut and not direction cut.
2: Events surviving gamma/hadron separation and direction cuts. [DEFAULT]
'''

# Definition of TtypeID within eventDisplay:
TELESCOPE_TYPES = {
    138704810: 'LST',
    10408618: 'MST-FlashCam',
    10608418: 'MST-NectarCam',
    201409917: 'SST',
    909924: 'SST-DC',
}

REQUIRED_NTTYPE = {
    'lapalma': 2,
    'paranal': 3,
}

TELESCOP = {
    'lapalma': 'CTA-N.BL-4LSTs15MSTs',
    'paranal': 'CTA-S.BL-4LSTs25MSTs70SSTs',
}

ALTITUDE = {
    'lapalma': 2158,
    'paranal': 2147,
}

# Common Header Keys
HDUCLASS = 'GADF'
HDUDOC = 'https://gamma-astro-data-formats.readthedocs.io'
HDUVERS = '0.2'



@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('filename', type=click.Path(exists=True, dir_okay=False))
@click.option(
    '-c', '--cut-level',
    default=2, type=click.IntRange(min=0, max=2),
    help=cuts_info
)
@click.option('--debug', '-d', is_flag=True)
@click.option('-o', '--output')
@click.argument('site', type=click.Choice(['paranal', 'lapalma']))
def cli(filename, cut_level, debug, output, site):
    """
    Command line tool for converting ED root files to DL2 fits files
    """

    if output is None:
        click.secho("No output file specified.", fg='yellow')
        click.secho("We will use the same filename, changing the extension to fits.", fg='yellow')
        output = filename.replace(".root", ".fits.gz")

    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.debug("Importing dependencies.")
    import uproot
    from astropy import units as u
    from astropy.io import fits
    from astropy.table import Table
    import numpy as np

    logging.debug("Opening Eventdisplay ROOT file and extracting content.")
    particle_file = uproot.open(filename)
    if 'hEmcUW' in particle_file:
        mc_energy_hist = particle_file['hEmcUW']
    else:
         mc_energy_hist = particle_file['hEmc']
    bin_content, bin_edges = mc_energy_hist.numpy()
    cuts = particle_file["fEventTreeCuts"]
    data = particle_file["data"]

    # Identify the telescope types within the 'NImages_Ttype' array:
    tel_types = [TELESCOPE_TYPES[t] for t in data.array("TtypeID")[0]]
    print('File contains the following telescope types: {}'.format(tel_types))

    cut_class = cuts.array('CutClass')
    # Cut 1: Events surviving gamma/hadron separation and direction cuts:
    mask_gamma_like_and_direction = cut_class == 5

    # Cut 2: Events surviving gamma/hadron separation cut and not direction cut
    mask_gamma_like_no_direction = cut_class == 0

    # Cut 0: Events before gamma/hadron and direction cuts (classes 0, 5 and 7)
    mask_before_cuts = mask_gamma_like_no_direction | mask_gamma_like_and_direction
    mask_before_cuts = mask_before_cuts | (cut_class == 7)

    if cut_level == 0:
        data_mask = mask_before_cuts
    elif cut_level == 1:
        data_mask = mask_gamma_like_no_direction
    elif cut_level == 2:
        data_mask = mask_gamma_like_and_direction


    logging.info(f'Surviving events: {np.count_nonzero(data_mask)}')

    # Remove events with NTtype!=2 in case of La Palma, and NTtype!=3 for Paranal.
    required_nttype = REQUIRED_NTTYPE[site]
    data_mask = data_mask & (data.array("NTtype") == required_nttype)

    # columns readable without transformation
    EVENTS_COLUMNS = {
        'OBS_ID': ('runNumber', None),
        'EVENT_ID': ('eventNumber', None),
        'MC_AZ': ('MCaz', u.deg),
        'AZ': ('Az', u.deg),
        'MC_ENERGY': ('MCe0', u.TeV),
        'ENERGY': ('ErecS', u.TeV),
        'MULTIP': ('NImages', None),
    }

    events = Table()
    for fits_key, (ed_key, unit) in EVENTS_COLUMNS.items():
        if unit is not None:
            events[fits_key] = u.Quantity(
                data.array(ed_key)[data_mask], unit, copy=False
            )
        else:
            events[fits_key] = data.array(ed_key)[data_mask]

    events['ALT'] = u.Quantity(
        90 - data.array("Ze")[data_mask], u.deg, copy=False
    )
    events['MC_ALT'] = u.Quantity(
        90 - data.array("MCze")[data_mask], u.deg, copy=False
    )
    events['GH_MVA'] = cuts.array('MVA')[data_mask]
    events['PNT_ALT'] = u.Quantity(
        data.array("ArrayPointing_Elevation")[data_mask], u.deg, copy=False
    )
    events['PNT_AZ'] = u.Quantity(
        data.array("ArrayPointing_Azimuth")[data_mask], u.deg, copy=False
    )

    logging.debug("Creating HDUs to be contained within the fits file.")

    # Create primary HDU:
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['TELESCOP'] = TELESCOP[site], 'Telescope and array codename'
    primary_hdu.header['COMMENT'] = "FITS (Flexible Image Transport System) format is defined in 'Astronomy"
    primary_hdu.header['COMMENT'] = "and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"

    # Create HDU
    header = fits.Header()
    header['HDUCLASS'] = HDUCLASS, 'This FITS file follows the GADF data format'
    header['HDUDOC'] = HDUDOC
    header['HDUVERS'] =  HDUVERS, 'Specification version'
    header['HDUCLAS1'] = 'EVENTS', 'Primary extension class'
    header['TELESCOP'] = TELESCOP[site]
    header['CREATOR'] = 'Eventdisplay prod5-v06'
    header['ORIGIN'] = 'CTA', 'Data generated by G.Maier'
    header['EUNIT'] = 'TeV', 'energy unit'
    header['ALTITUDE'] = ALTITUDE[site], 'Altitude of array center [m]'

    events_hdu = fits.BinTableHDU(events, header=header, name='EVENTS')

    # Store the histogram of simulated events vs MC_ENERGY
    energy_low = u.Quantity(10**bin_edges[:-1], u.TeV, copy=False)
    energy_high = u.Quantity(10**bin_edges[1:], u.TeV, copy=False)
    num_events = bin_content

    simulated_events = Table()
    simulated_events['MC_ENERG_LO'] = energy_low
    simulated_events['MC_ENERG_HI'] = energy_high
    simulated_events['EVENTS'] = num_events

    # Store the histogram of simulated events:
    hdu2 = fits.BinTableHDU(data=simulated_events, name='SIMULATED EVENTS')
    hdu2.header.set(
        'LO_THRES', energy_low.min().to_value(u.TeV),
        'Low energy threshold of validity [TeV]'
    )
    hdu2.header.set(
        'HI_THRES', energy_high.max().to_value(u.TeV),
        'High energy threshold of validity [TeV]'
    )
    hdu2.header.set('CREF3', '(MC_ENERG_LO:MC_ENERG_HI)', '')

    run_header = particle_file['MC_runheader']
    run_header = {
        k.replace('_5f_', '_'): [getattr(run_header, '_' + k)]
        for k in run_header._fields
    }
    run_header_hdu = fits.BinTableHDU(data=Table(run_header), name='RUNHEADER')

    logging.debug("Generating fits file.")

    # Generate output fits file:
    hdulist = fits.HDUList([primary_hdu, events_hdu, hdu2, run_header_hdu])
    hdulist.writeto(output, overwrite=True)


if __name__ == '__main__':
    cli()
