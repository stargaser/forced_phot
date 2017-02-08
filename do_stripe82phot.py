import os
import functools32
from forcedPhotExternalCatalog import ForcedPhotExternalCatalogTask
from lsst.daf.persistence import Butler


if os.path.exists('/home/shupe/work/forcephot/output'):
    out_butler = Butler('/home/shupe/work/forcephot/output')
elif os.path.exists('/hydra/workarea/forcephot/output'):
    out_butler = Butler('/hydra/workarea/forcephot/output')


def doit(dataId, refCat):
    """ Perform forced photometry on dataId from repo_str at positions in refCat
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
