import os
import tempfile
import functools32
from forcedPhotExternalCatalog import ForcedPhotExternalCatalogTask
from lsst.daf.persistence import Butler
from astropy.table import Table


if os.path.exists('/home/shupe/work/forcephot/output'):
    out_butler = Butler('/home/shupe/work/forcephot/output')
elif os.path.exists('/hydra/workarea/forcephot/output'):
    out_butler = Butler('/hydra/workarea/forcephot/output')
ftask = ForcedPhotExternalCatalogTask(out_butler)


def conv_afwtable_astropy(afwtable):
    with tempfile.NamedTemporaryFile() as tf:
        afwtable.writeFits(tf.name)
        tf.flush()
        tf.seek(0)
        atab = Table.read(tf.name, hdu=1)
    return(atab)


def doit(dataId, refCat):
    """ Perform forced photometry on dataId from repo_str at positions in refCat
    """
    in_butler = get_in_butler()
    exposure = in_butler.get('calexp', dataId=dataId)
    expWcs = exposure.getWcs()

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
    atab = conv_afwtable_astropy(measCat)
    return(atab)

@functools32.lru_cache()
def get_in_butler(repo_str='/datasets/gapon/data/DC_2013/calexps'):
    return(Butler(repo_str))
