# Script to convert Eventdisplay DL2a output to FITS
# 
# - expected Eventdisplay output including gamma/hadron cuts
#
#  Script by T.Hassan
#
import click
import logging

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
cuts_info = 'Cut level of the events to be included.'
cuts_info += ' 0: Events before applying gamma/hadron and direction cuts.'
cuts_info += ' 1: Events surviving gamma/hadron separation cut and not direction cut.'
cuts_info += ' 2: Events surviving gamma/hadron separation and direction cuts. [DEFAULT]'

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('filename', metavar='<filename>', type=click.Path(exists=True))
@click.option('--cut_level', '-c', nargs=1, type=click.INT, help=cuts_info)
@click.option('--debug', '-d', is_flag=True)
@click.argument('output', metavar='<output>')
@click.argument('site', metavar='<site>')
def cli(filename, cut_level, debug, output, site):
    """Command line tool for converting ED root files to "DL2-like" fits files

    \b
    For the moment, just one mode is implemented:



    Note: One one mode can be used at a time.
    """
    if len(filename) == 0:
        click.echo(cli.get_help(click.Context(cli)))
        raise click.Abort()
    if len(filename) > 0 and len(output) == 0:
        import string
        click.secho("No output file specified.", fg='yellow')
        click.secho("We will use same filename, changing the extension to fits.", fg='yellow')
        output = string.replace(".root", ".fits")
    if cut_level is None:
        # If cut level is not set, set the default:
        click.secho("No cut level specified.", fg='yellow')
        click.secho("We will use the default: events surviving gamma/hadron separation and direction cuts.",
                    fg='yellow')
        cut_level = 2
    if int(cut_level) < 0 or int(cut_level) > 2:
        click.echo(cli.get_help(click.Context(cli)))
        click.secho("The cut level requested is not implemented.", fg='yellow')
        raise click.Abort()
    if len(site) == 0 or (site.lower().find("lapalma")<0 and site.lower().find("paranal")<0):
        click.echo(cli.get_help(click.Context(cli)))
        click.secho("Sites implemented: lapalma, paranal")
        raise click.Abort()

    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.debug("Importing dependencies.")
    import uproot
    import numpy as np
    from astropy.coordinates.angle_utilities import angular_separation
    from astropy import units as u
    from astropy.io import fits

    logging.debug("Opening ED ROOT file and extracting content.")
    particle_file = uproot.open(filename)
    mc_energy_hist = particle_file['hEmc;1']
    bin_content, bin_edges = mc_energy_hist.numpy()
    cuts = particle_file["fEventTreeCuts"]
    data = particle_file["data"]

    # Cut 1: Events surviving gamma/hadron separation and direction cuts:
    mask_gamma_like_and_direction = cuts.array('CutClass') == 5

    # Cut 2: Events surviving gamma/hadron separation cut and not direction cut:
    mask_gamma_like_no_direction = cuts.array('CutClass') == 0

    # Cut 0: Events before gamma/hadron and direction cuts (classes 0, 5 and 7):
    mask_before_cuts = np.logical_or(mask_gamma_like_no_direction, cuts.array('CutClass') == 5)
    mask_before_cuts = np.logical_or(mask_before_cuts, cuts.array('CutClass') == 7)

    if cut_level is 0:
        data_mask = mask_before_cuts
    elif cut_level is 1:
        data_mask = mask_gamma_like_no_direction
    elif cut_level is 2:
        data_mask = mask_gamma_like_and_direction

    obs_id = data.array("runNumber")[data_mask]         # obs_id = tb.Int16Col(dflt=-1, pos=0)
    event_id = data.array("eventNumber")[data_mask]     # event_id = tb.Int32Col(dflt=-1, pos=1)
    NTels_trig = data.array("NTrig")[data_mask]         # NTels_trig = tb.Int16Col(dflt=0, pos=2)
    NTels_reco = data.array("NImages")[data_mask]       # NTels_reco = tb.Int16Col(dflt=0, pos=3)
    if site.lower().find("paranal") >= 0:
        NTels_reco_lst = [images_type[2] for images_type in data.array("NImages_Ttype")[data_mask]]
                                                            # (2)  NTels_reco_lst = tb.Int16Col(dflt=0, pos=4)
        NTels_reco_mst = [images_type[1] for images_type in data.array("NImages_Ttype")[data_mask]]
                                                            # (1)  NTels_reco_mst = tb.Int16Col(dflt=0, pos=5)
        NTels_reco_sst = [images_type[0] for images_type in data.array("NImages_Ttype")[data_mask]]
                                                            # (0)  NTels_reco_sst = tb.Int16Col(dflt=0, pos=6)
    else:
        NTels_reco_lst = [images_type[1] for images_type in data.array("NImages_Ttype")[data_mask]]
                                                            # (2)  NTels_reco_lst = tb.Int16Col(dflt=0, pos=4)
        NTels_reco_mst = [images_type[0] for images_type in data.array("NImages_Ttype")[data_mask]]
                                                            # (1)  NTels_reco_mst = tb.Int16Col(dflt=0, pos=5)
    mc_energy = data.array("MCe0")[data_mask]  # mc_energy = tb.Float32Col(dflt=np.nan, pos=7)
    reco_energy = data.array("ErecS")[data_mask]  # reco_energy = tb.Float32Col(dflt=np.nan, pos=8)
    mc_alt = 90 - data.array("MCze")[data_mask]
    mc_az = data.array("MCaz")[data_mask]
    reco_alt = 90 - data.array("Ze")[data_mask]  # reco_alt = tb.Float32Col(dflt=np.nan, pos=9)
    reco_az = data.array("Az")[data_mask]  # reco_az = tb.Float32Col(dflt=np.nan, pos=10)

    pointing_elevation = data.array("ArrayPointing_Elevation")[data_mask]
    pointing_azimuth = data.array("ArrayPointing_Azimuth")[data_mask]

    # Angular separation bewteen the center of the camera and the reco direction.
    offset = angular_separation(  # offset = tb.Float32Col(dflt=np.nan, pos=11)
        pointing_azimuth * u.deg,  # az
        pointing_elevation * u.deg,  # alt
        reco_az * u.deg,
        reco_alt * u.deg,
    )

    xi = angular_separation(
        mc_az * u.deg, mc_alt * u.deg, reco_az * u.deg, reco_alt * u.deg
    )

    h_max = data.array("EmissionHeight")[data_mask]  # h_max = tb.Float32Col(dflt=np.nan, pos=18)
    reco_core_x = data.array("Xcore")[data_mask]  # reco_core_x = tb.Float32Col(dflt=np.nan, pos=19)
    reco_core_y = data.array("Ycore")[data_mask]  # reco_core_y = tb.Float32Col(dflt=np.nan, pos=20)
    mc_core_x = data.array("MCxcore")[data_mask]  # mc_core_x = tb.Float32Col(dflt=np.nan, pos=21)
    mc_core_y = data.array("MCycore")[data_mask]  # mc_core_y = tb.Float32Col(dflt=np.nan, pos=22)

    logging.debug("Filling event list dictionary.")
    evt_dict = dict()
    # Filling Event List
    evt_dict['EVENT_ID'] = event_id
    evt_dict['OBS_ID'] = obs_id
    # evt_dict['TIME'] = timeArr
    # evt_dict['RA'] = raArr
    # evt_dict['DEC'] = decArr
    evt_dict['MC_ALT'] = mc_alt
    evt_dict['MC_AZ'] = mc_az
    evt_dict['MC_ENERGY'] = mc_energy
    evt_dict['ALT'] = reco_alt
    evt_dict['AZ'] = reco_az
    evt_dict['ENERGY'] = reco_energy
    evt_dict['MULTIP'] = NTels_reco

    # Filling Header info
    evt_dict['ALT_PNT'] = pointing_elevation
    evt_dict['AZ_PNT'] = pointing_azimuth
    evt_dict['COREX'] = reco_core_x
    evt_dict['COREY'] = reco_core_y
    # evt_dict['COREX'] = reco_core_x
    # evt_dict['COREY'] = reco_core_y
    # evt_dict['TELLIST'] = produce_tel_list(telConfig)
    # evt_dict['N_TELS'] = len(telConfig['TelID'])
    # evt_dict['GEOLON'] = REFERENCE_LON
    # evt_dict['GEOLAT'] = REFERENCE_LAT
    if site.lower().find("paranal") >= 0:
        evt_dict['ALTITUDE'] = 2150
        evt_dict['B_FIELD'] = "23.1177 microT (-22.7127,0)"
    else:
        evt_dict['ALTITUDE'] = 2147
        evt_dict['B_FIELD'] = "38.728 microT (37.9449,0)"

    logging.debug("Creating HDUs to be contained within the fits file.")
    # Create primary HDU:
    hdu0 = fits.PrimaryHDU()
    hdu0.header.set('TELESCOP',
                    'CTA-S.3HB9-FD',
                    'Telescope and array codename.')
    # hdu0.header.set('LICENSE ',
    #                 '',
    #                 '')
    hdu0.header['COMMENT'] = "FITS (Flexible Image Transport System) format is defined in 'Astronomy"
    hdu0.header['COMMENT'] = "and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"

    # Create event list HDU:
    # Columns to be saved
    columns = [fits.Column(name='OBS_ID', format='1K', array=evt_dict['OBS_ID']),
               fits.Column(name='EVENT_ID', format='1K', array=evt_dict['EVENT_ID']),
               #            fits.Column(name='TIME', format='1D', array=evt_dict['TIME'], unit="s"),
               fits.Column(name='MC_ALT', format='1E', array=evt_dict['MC_ALT'], unit="deg"),
               fits.Column(name='MC_AZ', format='1E', array=evt_dict['MC_AZ'], unit="deg"),
               fits.Column(name='MC_ENERGY', format='1E', array=evt_dict['MC_ENERGY'], unit="TeV"),
               fits.Column(name='ALT', format='1E', array=evt_dict['ALT'], unit="deg"),
               fits.Column(name='AZ', format='1E', array=evt_dict['AZ'], unit="deg"),
               fits.Column(name='ENERGY', format='1E', array=evt_dict['ENERGY'], unit="TeV"),
               fits.Column(name='MULTIP', format='1J', array=evt_dict['MULTIP'])
               ]

    # Create HDU
    hdu1 = fits.BinTableHDU.from_columns(columns)
    hdu1.name = "EVENTS"

    # For filling HDU headers
    HDUCLASS = 'GADF'
    HDUDOC = 'https://gamma-astro-data-formats.readthedocs.io'
    HDUVERS = '0.2'
    RADECSYS = 'FK5'

    hdu1.header.set('HDUCLASS', HDUCLASS, 'This FITS file follows the GADF data format')
    hdu1.header.set('HDUDOC', HDUDOC)
    hdu1.header.set('HDUVERS', HDUVERS, 'Specification version')
    hdu1.header.set('HDUCLAS1', 'EVENTS', 'Primary extension class')

    # Filling Header
    hdu1.header.set('CREATOR', 'Eventdisplay g500')
    hdu1.header.set('ORIGIN', 'CTA Collaboration', 'Data generated by G.Maier')
    if site.lower().find("paranal") >= 0:
        hdu1.header.set('TELESCOP', 'CTA-S.3HB9-FD')
    else:
        hdu1.header.set('TELESCOP', 'CTA-Nb.3AL4-BN15')

    # hdu1.header.set('RA_PNT  ', evt_dict['RA_PNT'], 'observation position RA [deg]')
    # hdu1.header.set('DEC_PNT ', evt_dict['DEC_PNT'], 'observation position DEC [deg]')
    hdu1.header.set('ALT_PNT ', pointing_elevation[0], 'Average altitude of pointing [deg]')
    hdu1.header.set('AZ_PNT  ', pointing_azimuth[0], 'Average azimuth of pointing [deg]')

    # get the list of telescopes that participate in the event
    # hdu1.header.set('TELLIST',
    #                 evt_dict['TELLIST'],
    #                 'comma-separated list of tel IDs')
    # hdu1.header.set('N_TELS',evt_dict['N_TELS'],
    #                 'number of telescopes in event list')

    hdu1.header.set('EUNIT   ', 'TeV', 'energy unit')
    # hdu1.header.set('GEOLON  ', VTS_REFERENCE_LON, 'longitude of array center [deg]')
    # hdu1.header.set('GEOLAT  ', VTS_REFERENCE_LAT, 'latitude of array center [deg]')
    hdu1.header.set('ALTITUDE', evt_dict['ALTITUDE'], 'Altitude of array center [m]')

    # Store the histogram of simulated events vs MC_ENERGY
    energy_low = np.power(10, bin_edges[0:-1])
    energy_high = np.power(10, bin_edges[1:len(bin_edges)])
    num_events = bin_content

    x = np.array([(energy_low, energy_high, num_events)],
                 dtype=[('MC_ENERG_LO', '>f4', np.shape(energy_low)),
                        ('MC_ENERG_HI', '>f4', np.shape(energy_high)),
                        ('EVENTS', '>f4', np.shape(num_events))])

    # Store the histogram of simulated events:
    hdu2 = fits.BinTableHDU(data=x)
    hdu2.name = "SIMULATED EVENTS"
    hdu2.header.set('TUNIT1 ', 'TeV', "")
    hdu2.header.set('TUNIT2 ', 'TeV', "")
    hdu2.header.set('TUNIT3 ', '', "")
    hdu2.header.set('LO_THRES', np.power(10, min(bin_edges)),
                    'Low energy threshold of validity [TeV]')
    hdu2.header.set('HI_THRES', np.power(10, max(bin_edges)),
                    'High energy threshold of validity [TeV]')
    hdu2.header.set('CREF3', '(MC_ENERG_LO:MC_ENERG_HI)', '')

    logging.debug("Generating fits file.")
    # Generate output fits file:
    hdus = list()
    hdus.append(hdu0)
    hdus.append(hdu1)
    hdus.append(hdu2)
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(output, overwrite=True)


if __name__ == '__main__':
    cli()



