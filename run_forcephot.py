
import concurrent.futures
import lsst.log
from astropy.table import Table, vstack, join, Column
from astropy.time import Time
from dax_utils import imgserv_json_to_df, get_image_table
from stripe82phot import (do_phot, parse_phot_table, 
                          make_refcat, parse_phot_table)



def do_one(dataId):
    global src_cat
    try:
        rval = do_phot(dataId=dataId, refCat=src_cat)
    except:
        return
    return(rval)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="""
            Run forced photometry on specified images and output a time-series table
        """)
    parser.add_argument('ra', help='Right Ascension in decimal degees')
    parser.add_argument('dec', help='Right Ascension in decimal degees')
    parser.add_argument('filter_name', help='name of filter to subset input table [ugriz]')
    parser.add_argument('output_table', help='output table path in FITS format')
    args = parser.parse_args()
    ra = float(args.ra)
    dec = float(args.dec)
    filter_name = args.filter_name
    output_table = args.output_table

    lsst.log.setLevel('', lsst.log.INFO)

    df_f = get_image_table(ra, dec, filter_name, 'Science_Ccd_Exposure')
    ids = [{'run':row.run, 'field':row.field, 'camcol':row.camcol, 
            'filter':row.filterName.encode()} 
            for index, row in df_f.iterrows()]

    global src_cat
    src_cat = make_refcat([ra], [dec])


    with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
        results = executor.map(do_one, ids)

    afwtabs = [r for r in results if r is not None]

    tbl_list = [parse_phot_table(t) for t in afwtabs]

    alltabs = vstack(tbl_list)

    # merge
    intab = Table.from_pandas(df_f)
    outtab = join(alltabs, intab, keys=['run','camcol','field','filterName'],
                  join_type='left')
    t = Time(outtab['expMidpt'], format='isot', scale='utc')
    outtab['mjd'] = t.mjd
    outtab.sort('mjd')
    outtab['ra'] = outtab['coord_ra'].to('deg')
    outtab['dec'] = outtab['coord_dec'].to('deg')
    mycols = ['mjd','psfMag','psfMagErr','ra','dec',
           'expMidpt','run','field','camcol','filterName',
           'objectId','base_RaDecCentroid_x',
           'base_RaDecCentroid_y','base_PsfFlux_flux','base_PsfFlux_fluxSigma']
    outtab.keep_columns(mycols)
    newtab = outtab[mycols]
    newtab.write(output_table, format='ipac')

if __name__ == '__main__':
    main()
