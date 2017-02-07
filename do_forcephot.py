from forcedPhotExternalCatalog import ForcedPhotExternalCatalogTask
from lsst.daf.persistence import Butler

def doit(repo_str, out_repo, dataId, refCat):
    """ Perform forced photometry on dataId from repo_str at positions in refCat
    """
    obut = Butler(out_repo)
    ftask = ForcedPhotExternalCatalogTask(obut)
    butler = Butler(repo_str)
    exposure = butler.get('calexp', dataId=dataId)
    expWcs = exposure.getWcs()

    measCat = ftask.measurement.generateMeasCat(exposure. refCat, expWcs)

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

