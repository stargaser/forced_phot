#!/usr/bin/env python

from __future__ import print_function, division

import numpy as np
import os

import lsst
import lsst.meas.base as measBase
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.pex.exceptions as pexExcept
import lsst.pipe.base as pipeBase
from lsst.afw.geom import Angle, degrees


from lsst.pipe.tasks.getRepositoryData import DataRefListRunner

# FrocedRaDecCentroidConfig & Plugin are lifted from Michael Wood-Vasey's
# branch of meas_base: https://github.com/lsst/meas_base/tree/u/wmwv/radeccentroid
#
# They're transplanted here for convenience and so that we don't have to set
# up an old version of meas_base.
class ForcedRaDecCentroidConfig(measBase.ForcedPluginConfig):
    pass


# Note the use of the measBase.register decorator to register the plugin for
# use with the framework and to give it a convenient name.
#
# Normally we use the prefix "base_" for plugins which are defined in
# meas_base itself, and "ext_" for plugins which are defined in extension
# packages (meas_extensions_simpleShape, meas_extensions_shapeHSM, etc). Here
# I'm using "suit_" just to keep things separate. The names are a matter of
# convention, and you can actually use whatever you like.
@measBase.register("suit_RaDecCentroid")
class ForcedRaDecCentroidPlugin(measBase.ForcedPlugin):
    """A centroid pseudo-algorithm for forced measurement that simply transforms the centroid
    from the reference catalog RA, Dec to the measurement coordinate system.  This is used as
    the slot centroid for external-catalog-based forced measurement, allowing subsequent measurements
    to simply refer to the slot value just as they would in single-frame measurement.
    """

    ConfigClass = ForcedRaDecCentroidConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.CENTROID_ORDER

    def __init__(self, config, name, schemaMapper, metadata):
        measBase.ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
        schema = schemaMapper.editOutputSchema()
        # Allocate x and y fields, join these into a single FunctorKey for ease-of-use.
        xKey = schema.addField(name + "_x", type="D", doc="transformed reference centroid column",
                               units="pixel")
        yKey = schema.addField(name + "_y", type="D", doc="transformed reference centroid row",
                               units="pixel")
        self.centroidKey = lsst.afw.table.Point2DKey(xKey, yKey)
        # Because we're taking the reference position as given, we don't bother transforming its
        # uncertainty and reporting that here, so there are no sigma or cov fields.  We do propagate
        # the flag field, if it exists.
        if "slot_Centroid_flag" in schemaMapper.getInputSchema():
            self.flagKey = schema.addField(name + "_flag", type="Flag",
                                           doc="whether the reference centroid is marked as bad")
        else:
            self.flagKey = None

    def measure(self, measRecord, exposure, refRecord, refWcs):
        targetWcs = exposure.getWcs()
        targetPos = targetWcs.skyToPixel(refRecord.getCoord())
        measRecord.set(self.centroidKey, targetPos)


# This TaskRunner does TWO special things:

# - It passes the dataset & coord_file arguments to the task's run() method;
# - It calls the task *once*, passing it all of the dataRefs as list, rather
#   than repeatedly, once per dataRef (which is the default behaviour).
class TaskRunnerWithArgs(DataRefListRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        result = [(parsedCmd.id.refList, {"dataset": parsedCmd.dataset, "coord_file": parsedCmd.coord_file})]
        return result

    def __call__(self, dataRefList):
        task = self.TaskClass(config=self.config, log=self.log)
        result = task.run(dataRefList[0], **dataRefList[1])

        if self.doReturnResults:
             return pipeBase.Struct(
                 dataRefList=dataRefList,
                 metadata=task.metadata,
                 result=result,
                 )
        return result


def load_external_catalog_info(coord_file):
    """Load name, ra, dec information from an external csv catalog.

    >>> tmpfile = 'foo.csv'
    >>> with open(tmpfile, 'w') as f:
    ...     f.write("Name,RA,Dec\\n")
    ...     f.write("SN2011gy,52.397163,40.867481")
    >>> foo = load_external_catalog_info(tmpfile)
    >>> print(foo)  # doctest: +NORMALIZE_WHITESPACE
      Name       RA       Dec
    -------- --------- ---------
    SN2011gy 52.397163 40.867481

    """
    from astropy.table import Table

    info = Table.read(coord_file, format='ascii.csv', names=('Name', 'RA', 'Dec'))
    return info


class ForcedExternalCatalogMeasurementTask(measBase.ForcedMeasurementTask):
    def attachTransformedFootprints(self, sources, refCat, exposure, expWcs):
        """Default implementation for attaching Footprints to blank sources prior to measurement
        Footprints for forced photometry must be in the pixel coordinate system of the image being
        measured, while we want to use RA, Dec position from an external catalog.
        This implementation takes RA, Dec from an external catalog, a fixed radius,
        and creates and sets footprints in the exposure's pixel system.
        Note that ForcedPhotImageTask delegates to this method in its own attachFootprints method.
        attachFootprints can then be overridden by its subclasses to define how their Footprints
        should be generated.
        See the documentation for run() for information about the relationships between run(),
        generateMeasCat(), and attachTransformedFootprints().
        """
        footprint_radius = 5  # pixels

        for srcRecord, refRecord in zip(sources, refCat):

            # Add footprints
            #  See https://community.lsst.org/t/how-do-i-do-forced-photometry-on-a-set-of-ra-dec/1074/9
            # From TallJimbo (Jim Bosch)
            # "There's a Footprint constructor that takes an integer position and a radius,
            #  so I think something like this should work:"
            # coord = lsst.afw.coord.IcrsCoord(lsst.afw.geom.Point2D(ra, dec), lsst.afw.geom.degrees)
            # fpCenter = lsst.afw.geom.Point2I(wcs.skyToPixel(coord))
            # footprint = lsst.afw.detection.Footprint(fpCenter, radius)

            # coord = afwCoord.IcrsCoord(afwGeom.Point2D(row['RA'], row['Dec']), degrees)
            coord = refRecord.getCoord()
            fpCenter = afwGeom.Point2I(expWcs.skyToPixel(coord))
            footprint = afwDetection.Footprint(fpCenter, footprint_radius)
            srcRecord.setFootprint(footprint)


class ForcedPhotExternalCatalogConfig(lsst.pex.config.Config):
    """Config class for forced measurement driver task."""

    measurement = lsst.pex.config.ConfigurableField(
        target=ForcedExternalCatalogMeasurementTask,
        doc="subtask to do forced measurement")

    def setDefaults(self):
        # RaDecCentroid takes the centroid from the reference catalog and uses it.
        # Note that we are using the "suit_RaDecCentroid" plugin registered
        # above, as well as PsfFlux (which is defined in meas_base, and is
        # automatically registered when that package is imported).
        self.measurement.plugins.names = ["suit_RaDecCentroid", "base_PsfFlux"]
        self.measurement.slots.shape = None
        self.measurement.slots.centroid = "suit_RaDecCentroid"


class ForcedPhotExternalCatalogTask(pipeBase.CmdLineTask):

    ConfigClass = ForcedPhotExternalCatalogConfig
    RunnerClass = TaskRunnerWithArgs
    _DefaultName = "ForcedPhotExternalCatalogTask"

    def __init__(self, butler=None, **kwargs):
        super(lsst.pipe.base.CmdLineTask, self).__init__(**kwargs)

        # We need to generate a minimal SourceTable schema which the
        # measurement subtask can work with.
        # In this context, "schema" means the list of fields in the table.
        # makeMinimalSchema() will provide the bare-bones that we need, with
        # fields for the source ID, position and any parents. Thus:
        #
        # >>> afwTable.SourceTable.makeMinimalSchema()
        # Schema(
        #     (Field['L'](name="id", doc="unique ID"), Key<L>(offset=0, nElements=1)),
        #     (Field['Angle'](name="coord_ra", doc="position in ra/dec"), Key<Angle>(offset=8, nElements=1)),
        #     (Field['Angle'](name="coord_dec", doc="position in ra/dec"), Key<Angle>(offset=16, nElements=1)),
        #     (Field['L'](name="parent", doc="unique ID of parent source"), Key<L>(offset=24, nElements=1)),
        # )
        self.refSchema = afwTable.SourceTable.makeMinimalSchema()
        self.makeSubtask("measurement", refSchema=self.refSchema)

    def create_source_catalog_from_external_catalog(self, dataRef, sources, debug=False):
        src_cat = afwTable.SourceCatalog(afwTable.SourceTable.makeMinimalSchema())
        for row in sources:
            record = src_cat.addNew()
            record.set('coord_ra', Angle(row['RA']*degrees))
            record.set('coord_dec', Angle(row['Dec']*degrees))

        if debug:
            print(src_cat['coord_ra'], src_cat['coord_dec'])
        return(src_cat)

    def __process_one_dataRef(self, dataRef, sources):
        butler = dataRef.getButler()
        exposure = butler.get(self.dataset, dataId=dataRef.dataId)
        expWcs = exposure.getWcs()
        refCat = self.create_source_catalog_from_external_catalog(dataRef, sources)
        measCat = self.measurement.generateMeasCat(exposure, refCat, expWcs)
        self.log.info("Performing forced measurement on science image %s" % (dataRef.dataId))

        try:
            self.measurement.attachTransformedFootprints(measCat, refCat, exposure, expWcs)
            self.measurement.run(measCat, exposure, refCat, expWcs)
        except pexExcept.InvalidParameterError as e:
            # This likely means that our Footprint fell off the edge of the
            # Exposure. We'll just warn and move on.
            # Note that we actully provide a ForcedPhotCcdTask which can take
            # care of this for you by deliberately only performing
            # measurements for which data is available. See the logic here:
            # https://github.com/lsst/meas_base/blob/master/python/lsst/meas/base/forcedPhotCcd.py#L105
            # I think trying to adopt it here makes this example more
            # complicated than it need to be for the moment.
            self.log.warn(str(e))

        # Our measurement catalog includes flux in raw units (ie, counts).
        # Let's convert them to magnitudes and store them in the output
        # catalog. Note that in production we'd actually use a separate task
        # to post-process the output of measurement to calibrated form (see
        # https://github.com/lsst/pipe_tasks/blob/master/python/lsst/pipe/tasks/transformMeasurement.py),
        # but we'll do this here for the sake of providing a nice example.
        #
        # Note that we can't add columns to a table that already exists (it
        # involves shifting things about in memory & would be very
        # inefficient), so we're going to make a copy of the existing measCat
        # with some new columns added.
        #
        # We can use a SchemaMapper to map fields in the existing table to the
        # new one.
        mapper = afwTable.SchemaMapper(measCat.schema)

        # Map through all the existing fields
        mapper.addMinimalSchema(measCat.schema)

        # Add two new ones containing floats
        mapper.editOutputSchema().addField("base_PsfFlux_mag", doc="magnitude derived from PsfFlux", type="F")
        mapper.editOutputSchema().addField("base_PsfFlux_magErr", doc="error on magnitude derived from PsfFlux", type="F")

        # Now create a new catalog and map over the fields from the old one
        newCat = afwTable.SourceCatalog(mapper.getOutputSchema())
        newCat.extend(measCat, mapper=mapper)

        # Get magnitude information so it can be added to catalog metadata
        calib = exposure.getCalib()
        fluxMag0, fluxMag0Err = calib.getFluxMag0()

        # Using getColumnView() we can operate on an entire column of the
        # table at once. Since our tables here have length 1 this isn't buying
        # us much, but with longer tables it could be much more efficient.
        #
        # Here, we use it calculate the magnotides for all sources in the
        # measCat...
        magColumn, magErrColumn = calib.getMagnitude(measCat.getColumnView()['base_PsfFlux_flux'],
                                                     measCat.getColumnView()['base_PsfFlux_fluxSigma'])

        # ...and then set those columns in our output table.
        newCat.getColumnView()['base_PsfFlux_mag'] = magColumn
        newCat.getColumnView()['base_PsfFlux_magErr'] = magErrColumn

        # It's usually idiomatic in Tasks to return Structs, which are
        # basically a home-brewed version of a named tuple.
        return pipeBase.Struct(measCat=newCat, fluxMag0=fluxMag0, fluxMag0Err=fluxMag0Err)


    def run(self, dataRefs, coord_file=None, dataset=None):
        """ Perform forced photometry on the dataRef exposure at the locations in coord_file.
        """
        # Specifying the dataset type as an argument to the task works fine in
        # this case.
        # Usually, though, we hard-code the dataset type a particular task can
        # operate on within the task itself, and have multiple different types
        # of tasks for the different datasets.
        #
        # This, for example, we have separate ForcedPhotCcd and
        # ForcedPhotCoadd tasks, despite the fact that fundamentally they do
        # the same thing, ie perform forced photometry on an image.
        #
        # The reason for this is reproducibility. When a task is run and
        # produces an output dataset, we store the task configuration in the
        # output repository as well as the output data. We then guarantee that
        # that is *the* result of running that task and storing the result in
        # that repository: there can be no ambiguity.
        #
        # Using your model, where datasets are specified on the command line,
        # you would have to change task configurations to process (say)
        # calexps and coadds, and store the results of two different
        # executions of the same task with different configurations to the
        # output repository. That makes provenance much harder to trace.
        self.dataset = dataset

        sources = load_external_catalog_info(coord_file)

        # We'll process each dataRef separately, generating a catalog for each.
        #
        # For simplicity, we're just using a list comprehension for this here.
        # But note that the framework actually provides primitives for
        # providing distributed processing of Tasks like this across a cluster
        # (or just using multiple CPUs on a given machine). Look at the
        # ctrl_pool package for the framework, or pipe_drivers for examples.
        # (Unfortunately, neither of them are well documented, and we expect
        # them to be replaced entirely over the next several months)
        catalogs = [self.__process_one_dataRef(dataRef, sources).measCat
                    for dataRef in dataRefs]

        # Next, we'll concatenate all those catalogs into a single master
        # catalog (by simply appending them to the first one in the list).
        # Note that in doing this we're actually losing ID information for
        # each measurement: there's no way to track which exposure the
        # measurement was made on. To track that, you could add another column
        # (or columns) to the tables generated in __process_one_dataRef, and
        # copy in the metadata you want to track from the Exposure.
        measCat = catalogs.pop(0)
        for catalog in catalogs:
            measCat.extend(catalog, deep=True)

        self.writeOutput(measCat)

    def writeOutput(self, sources):
        """Write the catalog to FITS"""
        # Note that we are explicitly NOT writing through the Butler (ie,
        # calling dataRef.put(). That's because the results we're collecting
        # here aren't semantically well-suited to being stored as a Butler
        # dataset. We normally expect a Butler dataset to represent a
        # reproducible result: for example, the results of performing
        # measurement on a particular image. However, in this case we have
        # performed measurement on an effectively arbitrary collection of
        # image data in response to user input. This doesn't fit naturally
        # into the Butler hierarchy as we currently conceive it.
        #
        # In future, there are a couple of approaches we could consider:
        #
        # - Creating a per-query output repository. Then you could have a
        #   "user query result" dataset within the repository which would have
        #   an unambiguous identity: the result of running the particular
        #   query associated with the repository.
        # - Dynamically create new Butler dataset types on the fly. Thus, you
        #   could create a new output dataset for each query received and
        #   retrieve it on demand. This isn't currently supported by the
        #   Butler, but it is planned (DM-4180).
        #
        # Neither of the above are concrete proposals for how you'd actually
        # want to do this in practice. :-)
        #
        # For now, we'll just write to a FITS file.
        sources.writeFits(os.path.join("output", "result.fits"))

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)

        # Can I make an argument which is a dataset type?
        parser.add_id_argument("--id", "src", help="data ID of the image")
        parser.add_argument("--dataset", default="calexp",
                            help="dataset to perform forced photometry on: 'calexp', 'deepDiff'")
        parser.add_argument("--coord_file",
                            help="File with coordinates to photometry. " +
                            "Each line should be Name,RA,Dec with RA, Dec as decimal degrees.")
        return parser

    # Overriding these two functions prevent the task from attempting to save the config.
    def _getConfigName(self):
        return None

    def _getMetadataName(self):
        return None


if __name__ == "__main__":

    ForcedPhotExternalCatalogTask.parseAndRun()
