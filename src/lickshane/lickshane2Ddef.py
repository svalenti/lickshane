import numpy as np
import pylab as pl
from astropy.io import fits as pyfits


def checkwavelength_arc(xx1, yy1, xx2, yy2, xmin, xmax, inter=True):
    # print "LOGX:: Entering `checkwavelength_arc` method/function in
    # %(__file__)s" % globals()
    minimo = max(min(xx1), min(xx2)) + 50
    massimo = min(max(xx1), max(xx2)) - 50
    yy1 = [0 if e < 0 else e for e in np.array(yy1)]
    yy2 = [0 if e < 0 else e for e in np.array(yy2)]
    _shift, integral = [], []
    for shift in range(-500, 500, 1):
        xxnew = xx1 + shift / 10.
        yy2interp = np.interp(xxnew, xx2, yy2)
        yy2timesyy = yy2interp * yy1
        xxcut = np.compress((np.array(xxnew) >= minimo) & (
            np.array(xxnew) <= massimo), np.array(xxnew))
        yycut = np.compress((np.array(xxnew) >= minimo) & (
            np.array(xxnew) <= massimo), np.array(yy2timesyy))
        integrale = np.trapz(yycut, xxcut)
        integral.append(integrale)
        _shift.append(shift / 10.)
    result = _shift[integral.index(max(integral))]
    if inter:
        # import matplotlib as mpl
        #   mpl.use("TKAgg")
        import pylab as pl
        pl.ion()
        pl.clf()
        ratio = np.trapz(yy1, xx1) / np.trapz(yy2, xx2)
        yy3 = np.array(yy2) * float(ratio)
        xx4 = xx1 + result
        pl.plot(xx1, yy1, label='spectrum')
        pl.plot(xx2, yy3, label='reference sky')
        pl.plot(xx4, yy1, label='shifted spectrum')
        pl.legend(numpoints=1, markerscale=1.5)
        if xmin != '' and xmax != '':
            pl.xlim(xmin, xmax)
    return result


def continumsub(imagefile, _order1, _order2):
    # print "LOGX:: Entering `continumsub` method/function in %(__file__)s" %
    # globals()
    import lickshane
    from pyraf import iraf
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.continuum']
    for t in toforget:
        iraf.unlearn(t)
    lickshane.util.delete('tsky.fits')
    iraf.specred.continuum(imagefile, output='tsky.fits', type='difference',
                           interact='no', function='legendre', niterat=300, low_rej=3, high_re=2, sample='*',
                           order=_order1, ask='YES')
    lickshane.util.delete(imagefile)
    iraf.continuum('tsky.fits', output=imagefile, type='difference',
                   interact='no', function='spline1', overrid='yes', niterat=10, low_rej=3, high_re=1, sample='*',
                   order=_order2, ask='YES')
    lickshane.util.delete('tsky.fits')
    return imagefile


def skyfrom2d(fitsfile, skyfile, interac=True):
    # print "LOGX:: Entering `skyfrom2d` method/function in %(__file__)s" %
    # globals()
    import lickshane
    tpe = pyfits.open(fitsfile)[0].header.get('version')

    if tpe in ['kastr']:
        yy1 = pyfits.open(fitsfile)[0].data[:, :].mean(1)
        crval2 = pyfits.open(fitsfile)[0].header.get('CRVAL2')
        cd2 = pyfits.open(fitsfile)[0].header.get('CD2_2')
        
        lickshane.util.delete('new3.fits')
        hdu = pyfits.PrimaryHDU(yy1)
        hdulist = pyfits.HDUList([hdu])
        hdulist[0].header.update('CRVAL1', crval2)
        hdulist[0].header.update('CD1_1', cd2)
        hdulist.writeto('new3.fits')
        hdulist.close()
        
        fitsfile = lickshane.lickshane2Ddef.continumsub('new3.fits', 6, 1)
        yy1 = pyfits.open(fitsfile)[0].data
        xx1 = np.arange(len(yy1))
        aa1 = crval2 + (xx1) * cd2
#        lickshane.util.delete('new3.fits')
        rangea,rangeb = 5000,10000
    else:
        # need to check
        yy1 = pyfits.open(fitsfile)[0].data[:, :].mean(0)
        crval1 = pyfits.open(fitsfile)[0].header.get('CRVAL1')
        cd1 = pyfits.open(fitsfile)[0].header.get('CD1_1')
        
        lickshane.util.delete('new3.fits')
        hdu = pyfits.PrimaryHDU(yy1)
        hdulist = pyfits.HDUList([hdu])
        hdulist[0].header.update('CRVAL1', crval1)
        hdulist[0].header.update('CD1_1', cd1)
        hdulist.writeto('new3.fits')
        hdulist.close()

        fitsfile = lickshane.lickshane2Ddef.continumsub('new3.fits', 6, 1)
        yy1 = pyfits.open(fitsfile)[0].data
        xx1 = np.arange(len(yy1))
        aa1 = crval1 + (xx1) * cd1
        lickshane.util.delete('new3.fits')
        rangea,rangeb = 3500,5000

    skyff = pyfits.open(skyfile)[0].data
    crval1 = pyfits.open(skyfile)[0].header.get('CRVAL1')
    cd1 = pyfits.open(skyfile)[0].header.get('CD1_1')
    skyxx = np.arange(len(skyff))
    skyaa = crval1 + (skyxx) * cd1
    shift = lickshane.lickshane2Ddef.checkwavelength_arc(aa1, yy1, skyaa, skyff, rangea, rangeb, interac)
    return shift
