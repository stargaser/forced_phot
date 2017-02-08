
import lsst.afw.table as afwTable
from lsst.daf.persistence import Butler
from lsst.afw.geom import Angle, degrees

def make_refcat(ra, dec, repo_str='/home/shupe/work/forcephot/output'):
    out_butler = Butler(repo_str)
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


