def json2python(json):
    from simplejson import loads
    try:
        return loads(json)
    except:
        pass
    return None


def niceprint(x, thisindent='', indent=''):
    subindent = indent + '  '
    if type(x) is dict:
        print thisindent + '{'
        for (k,v) in x.items():
            print subindent + k + ' = ',
            niceprint(v, '', subindent)
        print indent + '}'
    elif type(x) is list:
        print thisindent + '['
        for k in x:
            niceprint(k, subindent, subindent)
        print indent + ']'
    else:
        print thisindent + repr(x)
