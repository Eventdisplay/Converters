# Script to convert Eventdisplay DL2a output to FITS
# following the open data formats for gamma-ray
# astronomy
# (https://gamma-astro-data-formats.readthedocs.io/en/latest/)
#
# - expect Eventdisplay output including gamma/hadron cuts
#
#  Script started by T.Hassan
#
import logging

import click
import numpy as np
import uproot
from astropy import units as u
from astropy.io import fits
from astropy.table import Table

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
cuts_info = """Cut level of the events to be included.
0: Events before applying gamma/hadron and direction cuts.
1: Events surviving gamma/hadron separation cut and not direction cut.
2: Events surviving gamma/hadron separation and direction cuts. [DEFAULT]
"""

# Common Header Keys
HDUCLASS = "GADF"
HDUDOC = "https://gamma-astro-data-formats.readthedocs.io"
HDUVERS = "0.2"


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("filenames", type=click.Path(exists=True, dir_okay=False), nargs=-1)
@click.option("-c", "--cut-level", default=2, type=click.IntRange(min=0, max=2), help=cuts_info)
@click.option("--debug", "-d", is_flag=True)
@click.option("-o", "--output")
@click.option("-l", "--layout", default="CTA-layout", help="layout name")
@click.option(
    "-t",
    "--event_type",
    nargs=1,
    type=click.INT,
    default=-1,
    help="Event type (in case they were computed)",
)
def cli(filenames, cut_level, debug, output, layout, event_type):
    """Convert Eventdisplay root files to DL2 fits files"""

    if output is None:
        click.secho("No output file specified.", fg="yellow")
        click.secho("Output file name will match input name (fits extension)", fg="yellow")
        if event_type > 0:
            output = filenames[0].replace(
                ".eff-0.root", "_event_type_{}.fits.gz".format(event_type)
            )
        else:
            output = filenames[0].replace(".eff-0.root", ".fits.gz")

    # If event_type is larger than 0, extract event-wise event types:
    event_types = None
    if event_type > 0:
        event_types = np.loadtxt(filenames[0].replace(".eff-0.root", ".txt"), dtype=float)
        event_types = event_types.astype(int)
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.debug("Opening Eventdisplay ROOT file and extracting content.")
    logging.info(f"Reading from {filenames}")
    particle_file = uproot.open(filenames[0])
    bin_content, bin_edges = particle_file["hEmcUW"].to_numpy(flow=True)
    # get min/max energy from MC run header
    # (as the under and overflow bin edges are -/+ inf)
    run_header = {
        k: [v]
        for k, v in particle_file["MC_runheader"].all_members.items()
        if k.find("@") != 0 and k != "fName" and k != "fTitle"
    }
    bin_edges[0] = np.log10(run_header["E_range"][0][0])
    bin_edges[-1] = np.log10(run_header["E_range"][0][1])

    site_altitude = uproot.open(filenames[0])["MC_runheader"].member("obsheight")

    # Combine a branch of several files:
    cut_class = uproot.lazy([{f: "DL2EventTree"} for f in filenames], "CutClass")
    cut_class = cut_class["CutClass"]
    # cut_class = particle_file['DL2EventTree/CutClass'].array(library='np')
    # Cut 1: Events surviving gamma/hadron separation and direction cuts:
    mask_gamma_like_and_direction = cut_class == 5

    # Cut 2: Events surviving gamma/hadron separation cut and not direction cut
    mask_gamma_like_no_direction = cut_class == 0

    # Cut 0: Events before gamma/hadron and direction cuts (classes 0, 5 and 7)
    mask_before_cuts = mask_gamma_like_no_direction | mask_gamma_like_and_direction
    mask_before_cuts = mask_before_cuts | (cut_class == 7)

    if cut_level == 0 or event_type > 0:
        data_mask = mask_before_cuts
    elif cut_level == 1:
        data_mask = mask_gamma_like_no_direction
    elif cut_level == 2:
        data_mask = mask_gamma_like_and_direction

    data_mask = data_mask.to_numpy()

    logging.info(f"Total number of events read from root file: {len(data_mask)}")
    logging.info(f"Total of selected events: {np.count_nonzero(data_mask)}")

    if event_type > 0:
        logging.info(f"Number of events within event_type file: {len(event_types)}")
        ratio = np.sum(event_types == -1) / len(event_types)
        logging.info(
            f"Ratio of training to simulated events: {ratio}"
        )
        logging.info(f"len(data_mask): {len(data_mask)}")
        data_mask[data_mask] = event_types == event_type
        if np.sum(data_mask) == 0:
            raise ValueError("No events to export of event types == {}".format(event_type))
        logging.info(f"Surviving events: {np.count_nonzero(data_mask)} (event type {event_type})")

    else:
        logging.info(f"Surviving events: {np.count_nonzero(data_mask)} (cut level {cut_level})")
        ratio = 0

    # columns readable without transformation
    EVENTS_COLUMNS = {
        "MC_AZ": ("MCaz", u.deg),
        "MC_ALT": ("MCel", u.deg),
        "AZ": ("az", u.deg),
        "ALT": ("el", u.deg),
        "MC_ENERGY": ("MCe0", u.TeV),
        "ENERGY": ("erec", u.TeV),
        "MULTIP": ("nimages", None),
        "PNT_AZ": ("ArrayPointing_Azimuth", u.deg),
        "PNT_ALT": ("ArrayPointing_Elevation", u.deg),
        "GH_MVA": ("MVA", None),
    }

    events = Table()
    for fits_key, (ed_key, unit) in EVENTS_COLUMNS.items():
        if unit is not None:
            events[fits_key] = u.Quantity(
                uproot.lazy([{f: "DL2EventTree"} for f in filenames], ed_key)[ed_key][data_mask],
                # particle_file['DL2EventTree/'+ed_key].array(library='np'),
                unit,
                copy=False,
            )
        else:
            events[fits_key] = uproot.lazy([{f: "DL2EventTree"} for f in filenames], ed_key)[
                ed_key
            ][data_mask]
            # events[fits_key] = \
            # particle_file['DL2EventTree/'+ed_key].array(library='np')[data_mask]

    # fake event numbers
    events["EVENT_ID"] = np.arange(len(events["MC_AZ"]))
    events["OBS_ID"] = uproot.lazy([{f: "DL2EventTree"} for f in filenames], "runNumber")[
        "runNumber"
    ][data_mask]
    # events['OBS_ID'] = particle_file['DL2EventTree/runNumber'].array(library='np')[data_mask]
    logging.debug("Creating HDUs to be contained within the fits file.")

    # Create primary HDU:
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header["TELESCOP"] = layout, "Telescope and array codename"
    primary_hdu.header[
        "COMMENT"
    ] = "FITS (Flexible Image Transport System) format is defined in 'Astronomy"
    primary_hdu.header[
        "COMMENT"
    ] = "and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"

    # Create HDU
    header = fits.Header()
    header["HDUCLASS"] = HDUCLASS, "This FITS file follows the GADF data format"
    header["HDUDOC"] = HDUDOC
    header["HDUVERS"] = HDUVERS, "Specification version"
    header["HDUCLAS1"] = "EVENTS", "Primary extension class"
    header["TELESCOP"] = layout
    header["CREATOR"] = "Eventdisplay prod5"
    header["ORIGIN"] = "CTA", "Data generated by G.Maier"
    header["EUNIT"] = "TeV", "energy unit"
    header["ALTITUDE"] = site_altitude, "Altitude of array center [m]"

    events_hdu = fits.BinTableHDU(events, header=header, name="EVENTS")

    # Store the histogram of simulated events vs MC_ENERGY
    energy_low = u.Quantity(10 ** bin_edges[:-1], u.TeV, copy=False)
    energy_high = u.Quantity(10 ** bin_edges[1:], u.TeV, copy=False)
    # If there are events used in the event-types training, that proportion must be
    # excluded from the number of events in the simulation, so the IRFs can be
    # then computed properly
    num_events = bin_content * (1 - ratio)

    simulated_events = Table()
    simulated_events["MC_ENERG_LO"] = energy_low
    simulated_events["MC_ENERG_HI"] = energy_high
    simulated_events["EVENTS"] = num_events

    # Store the histogram of simulated events:
    hdu2 = fits.BinTableHDU(data=simulated_events, name="SIMULATED EVENTS")
    hdu2.header.set(
        "LO_THRES", energy_low.min().to_value(u.TeV), "Low energy threshold of validity [TeV]"
    )
    hdu2.header.set(
        "HI_THRES", energy_high.max().to_value(u.TeV), "High energy threshold of validity [TeV]"
    )
    hdu2.header.set("CREF3", "(MC_ENERG_LO:MC_ENERG_HI)", "")

    run_header_hdu = fits.BinTableHDU(data=Table(run_header), name="RUNHEADER")

    logging.debug("Generating fits file.")

    # Generate output fits file:
    hdulist = fits.HDUList([primary_hdu, events_hdu, hdu2, run_header_hdu])
    logging.info(f"Writing to {output}")
    hdulist.writeto(output, overwrite=True)


if __name__ == "__main__":
    cli()
