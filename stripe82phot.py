from __future__ import print_function, division, absolute_import

""" These functions are designed to run on the PDAC server with access to
the Stripe82 data
"""

import os
import functools32
import tempfile
import numpy as np
from lsst.daf.persistence import Butler
import lsst.afw.table as afwTable
from lsst.afw.geom import Angle, degrees
from astropy.table import Table
import astropy.units as u
from forcedPhotExternalCatalog import ForcedPhotExternalCatalogTask

# Define this output Butler (needed only for schema)
if os.path.exists('/home/shupe/work/forcephot/output'):
    out_butler = Butler('/home/shupe/work/forcephot/output')
elif os.path.exists('/hydra/workarea/forcephot/output'):
    out_butler = Butler('/hydra/workarea/forcephot/output')


def make_refcat(ra, dec):
    """
    Make a reference catalog for forced photometry

    Parameters:
    -----------
    ra : sequence or array
        Right Ascension in decimal degrees
    dec : sequence or array
        Declination in decimal degrees

    Returns:
    --------
    src_cat : lsst.afw.table.tableLib.SourceCatalog
        Source catalog for the forced photometry task
    """
    schema = out_butler.get('src_schema', immediate=True).schema
    mapper = afwTable.SchemaMapper(schema)
    mapper.addMinimalSchema(schema)
    newSchema = mapper.getOutputSchema()
    src_cat = afwTable.SourceCatalog(newSchema)
    for row in zip(ra,dec):
        record = src_cat.addNew()
        record.set('coord_ra', Angle(row[0]*degrees))
        record.set('coord_dec', Angle(row[1]*degrees))
    return(src_cat)


def conv_afwtable_astropy(afwtable):
    """
    Convert an afwTable to an Astropy table

    Parameters:
    -----------
    afwtable : lsst.afw.table table
        input table

    Returns:
    --------
    atab : astropy.table.Table
        table in Astropy format
    """
    with tempfile.NamedTemporaryFile() as tf:
        afwtable.writeFits(tf.name)
        tf.flush()
        tf.seek(0)
        atab = Table.read(tf.name, hdu=1)
    return(atab)


def parse_phot_table(afwTable, convert=True):
    """
    Parse a forced photometry output table, optionally converting to Astropy format

    Parameters:
    -----------
    afwTable : lsst.afw.table table or astropy.table.Table
        Output of forced photometry task
    convert : bool
        Convert input to Astropy table? Default is True
        Only use this if input is already an Astropy table

    Returns:
    --------
    tab : astropy.table.Table
        converted table with metadata propagated to columns and magnitudes added
    """
    if convert:
        tab = conv_afwtable_astropy(afwTable)
    else:
        tab = afwTable
    tab['run'] = tab.meta['RUN']
    tab['camcol'] = tab.meta['CAMCOL']
    tab['field'] = tab.meta['FIELD']
    tab['filterName'] = tab.meta['FILTER']
    tab['psfMag'] = -2.5*np.log10(tab['base_PsfFlux_flux']/tab.meta['FLUXM0'])
    tab['psfMagErr'] = -2.5*np.log10(1.- + (tab['base_PsfFlux_fluxSigma']
                                     /tab['base_PsfFlux_flux']))
    tab['psfMag'].unit = u.mag
    tab['psfMagErr'].unit = u.mag
    del tab.meta['RUN']
    del tab.meta['CAMCOL']
    del tab.meta['FIELD']
    del tab.meta['FILTER']
    del tab.meta['FLUXM0']
    del tab.meta['FLUXM0SG']
    return(tab)


def do_phot(dataId, refCat):
    """ Perform forced photometry on specified data and reference catalog

    Call the ForcedPhotometryExternalCatalogTask to measure forced photometry
    on the specified image, at the positions in the reference catalog. The
    data are fetched from the repository cached by the get_in_butler function
    (to speed up forced photometry on many images)

    Parameters:
    -----------
    dataId : dict
        Dictionary of key:value pairs specifying data
    refCat : lsst.afw.table.tableLib.SourceCatalog 
        Reference catalog containing positions for forced photometry

    Returns:
    --------
    measCat : lsst.afw.table table
        Forced photometry task output
    """
    in_butler = get_in_butler()
    exposure = in_butler.get('calexp', dataId=dataId)
    expWcs = exposure.getWcs()

    ftask = ForcedPhotExternalCatalogTask(out_butler)
    measCat = ftask.measurement.generateMeasCat(exposure, refCat, expWcs)

    ftask.measurement.attachTransformedFootprints(measCat, refCat, exposure, expWcs)
    ftask.measurement.run(measCat, exposure, refCat, expWcs)

    # Get magnitude information so it can be added to catalog metadata
    calib = exposure.getCalib()
    fluxMag0, fluxMag0Err = calib.getFluxMag0()

    meta = measCat.getTable().getMetadata()
    for (key, val) in dataId.iteritems():
        meta.add(key.upper(), val)
    meta.add('FLUXM0', fluxMag0)
    meta.add('FLUXM0SG', fluxMag0Err)
    measCat.getTable().setMetadata(meta)
    return(measCat)


@functools32.lru_cache()
def get_in_butler(repo_str='/datasets/gapon/data/DC_2013/calexps'):
    return(Butler(repo_str))

