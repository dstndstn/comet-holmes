import urllib2
import urllib
import sys
import threading
import time

from util import *

# See:
#  http://code.google.com/apis/ajaxsearch/documentation/reference.html#_restUrlBase
#
# 'referer' should be *your* site
# 'imgsize' can be: "icon", "small", "medium", "large", "xlarge",
#     "xxlarge", "huge", or any combo separate by "|".
#     default is no limit on size.
# 'imgtype' can be: "photo", "lineart", etc...
# 'rsz': result set size: 'small'=4, 'large'=8.
# 'filetype': can be 'jpg', 'png', 'gif', 'bmp'
#
def imagesearch(query, referer, imgsize=None, start=0, imgtype=None,
                rsz='large', key=None, filetype=None):
    url = 'http://ajax.googleapis.com/ajax/services/search/images?v=1.0'
    url += '&q=' + urllib.quote_plus(query)
    if imgsize:
        url += '&imgsz=' + imgsize
    if start:
        url += '&start=%i' % start
    if rsz:
        url += '&rsz=' + rsz
    if key:
        url += '&key=' + key
    if filetype:
        url += '&as_filetype=' + filetype

    req = urllib2.Request(url, None, {'Referer': referer})
    f = urllib2.urlopen(req)
    json = f.read()
    py = json2python(json)
    return py


def fetchall(*args1, **args2):
    mythreads = []

    start = 0
    while True:
        args2['start'] = start
        res = imagesearch(*args1, **args2)

        #niceprint(res)
        print
        print

        #print 'Starting at %i' % start
        rd = res['responseData']
        if rd is None:
            print 'No results'
            break
        cursor = rd['cursor']
        moreresultsurl = cursor['moreResultsUrl']
        estimatedResults = int(cursor['estimatedResultCount'])

        #print 'Estimated %i results.' % estimatedResults
        #print 'More results:', moreresultsurl

        results = rd['results']
        for i,r in enumerate(results):
            url = r['url']
            w = int(r['width'])
            h = int(r['height'])
            #print 'url: %ix%i' % (w,h), url
            print '# %i (%ix%i)' % (start+i, w, h)
            print url

            def download(url, fn):
                f = urllib2.urlopen(url)
                data = f.read()
                o = open(fn, 'w')
                o.write(data)
                o.close()

            t = threading.Thread(target=download, args=(url, 'gimage-%02i' % (start+i)))
            mythreads.append(t)
            t.start()


        start += len(results)

    while True:
        nalive = sum([t.isAlive() for t in mythreads])
        print 'Waiting for %i threads.' % nalive
        if nalive == 0:
            break
        time.sleep(10)



if __name__ == '__main__':

    args1 = ('Comet Holmes', 'http://astrometry.net/',)
    args2 = {
        imgtype='photo', filetype='jpg',
        key='ABQIAAAA7dWWcc9pB-GTzZE7CvT6SRTUmheot7Jkn4LXD-2luhklyOUjjRQsxt84r-N3HEJGAqD16UgFd5pp2g'
        }

    if False:
        args2['imgsize'] = 'huge'
        fetchall(*args1, **args2)

    if True:
        args2['imgsize'] = 'xxlarge'
        fetchall(*args1, **args2)

