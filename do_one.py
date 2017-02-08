

def do_one(dataId):
    try:
        rval = doit(dataId=dataId, refCat = src_cat)
    except:
        return
    return(rval)

