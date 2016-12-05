#  cosmic correction
#  lacosmic  iraf modules rewritten in pyraf
#
#

import os
import numpy as np
import math

try:       from astropy.io import fits as pyfits
except:    import pyfits

# We define the laplacian kernel to be used
laplkernel = np.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])
# Other kernels :
growkernel = np.ones((3, 3))

# FITS import - export


def fromfits(infilename, hdu=0, verbose=True):
    # print "LOGX:: Entering `fromfits` method/function in %(__file__)s" %
    # globals()
    """
        Reads a FITS file and returns a 2D numpy array of the data.
        Use hdu to specify which HDU you want (default = primary = 0)
        """

    pixelarray, hdr = pyfits.getdata(infilename, hdu, header=True)
    pixelarray = np.asarray(pixelarray).transpose()

    pixelarrayshape = pixelarray.shape
    if verbose:
        print "FITS import shape : (%i, %i)" % (pixelarrayshape[0], pixelarrayshape[1])
        print "FITS file BITPIX : %s" % (hdr["BITPIX"])
        print "Internal array type :", pixelarray.dtype.name

    return pixelarray, hdr


def tofits(outfilename, pixelarray, hdr=None, verbose=True):
    # print "LOGX:: Entering `tofits` method/function in %(__file__)s" %
    # globals()
    """
        Takes a 2D numpy array and write it into a FITS file.
        If you specify a header (pyfits format, as returned by fromfits()) it will be used for the image.
        You can give me boolean numpy arrays, I will convert them into 8 bit integers.
        """
    pixelarrayshape = pixelarray.shape
    if verbose:
        print "FITS export shape : (%i, %i)" % (pixelarrayshape[0], pixelarrayshape[1])

    if pixelarray.dtype.name == "bool":
        pixelarray = np.cast["uint8"](pixelarray)

    if os.path.isfile(outfilename):
        os.remove(outfilename)

    if hdr == None:  # then a minimal header will be created
        hdu = pyfits.PrimaryHDU(pixelarray.transpose())
    else:  # this if else is probably not needed but anyway ...
        hdu = pyfits.PrimaryHDU(pixelarray.transpose(), hdr)

    hdu.writeto(outfilename, output_verify='ignore')

    if verbose:
        print "Wrote %s" % outfilename

###################################################


def lacos_spec(_input, output='clean.fits', outmask='mask.fits', gain=1.3, readn=9,\
               xorder=9, yorder=3, sigclip=4.5, sigfrac=0.5, objlim=1, niter=4, instrument='kastr', verbose=True, interactive=False):
    # print "LOGX:: Entering `lacos` method/function in %(__file__)s" %
    # globals()
    import lickshane
    import sys
    import re
    import os
    import string
    from pyraf import iraf
    import numpy as np

    oldoutput, galaxy, skymod, med5 = 'oldoutput.fits', 'galaxy.fits', 'skymod.fits', 'med5.fits'
    blk, lapla, med3, med7, sub5, sigima, finalsel = 'blk.fits', 'lapla.fits', 'med3.fits', 'med7.fits', 'sub5.fits', 'sigima.fits', 'finalsel.fits'
    deriv2, noise, sigmap, firstsel, starreject = 'deriv2.fits', 'noise.fits', 'sigmap.fits', 'firstsel.fits', 'starreject.fits'
    inputmask = 'inputmask.fits'
    # set some parameters in standard IRAF tasks
    iraf.convolve.bilinear = 'no'
    iraf.convolve.radsym = 'no'
    # create Laplacian kernel
    # laplkernel = np.array([[0.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 0.0]])
    f = open('_kernel', 'w')
    f.write('0 -1 0;\n-1 4 -1;\n0 -1 0')
    f.close()
    # create growth kernel
    f = open('_gkernel', 'w')
    f.write('1 1 1;\n1 1 1;\n1 1 1')
    f.close()
    gkernel = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
    lickshane.util.delete(galaxy)
    lickshane.util.delete(skymod)
    lickshane.util.delete(oldoutput)


    if not output:
        output = _input
    else:
        os.system('cp '+_input+' '+output)
    os.system('cp '+_input+' '+oldoutput)

    arrayinput, headerinput = lickshane.cosmics.fromfits(oldoutput, verbose=False)
    lickshane.cosmics.tofits(outmask, np.float32(
        arrayinput - arrayinput), headerinput, verbose=False)



    if instrument in ['kastr']:
        axis1 = 1
        axis2 = 2
    elif instrument in ['kastb']:
        axis1 = 2
        axis2 = 1

    # subtract object spectra if desired
    if xorder > 0:
        iraf.fit1d(oldoutput, galaxy, "fit", axis= axis1, order=xorder, func="leg", low=4.,
                   high=4., nav=1, inter='no', sample="*", niter=3, grow=0, cursor="")
        iraf.imarith(oldoutput, "-", galaxy, oldoutput)
    else:
        lickshane.cosmics.tofits(galaxy, np.float32(
            arrayinput - arrayinput), headerinput, verbose=False)


    # Subtract sky lines
    if yorder >0:
        iraf.fit1d(oldoutput, skymod, "fit", axis= axis2, order=yorder, func="leg", low=4., high=4.,
                   inter='no', sample="*", nav=1, niter=3, grow=0, cursor="")
        iraf.imarith(oldoutput, "-", skymod, oldoutput)
    else:
        lickshane.cosmics.tofits(skymod, np.float32(
            arrayinput - arrayinput), headerinput, verbose=False)

    arrayoldoutput, headeroldoutput = lickshane.cosmics.fromfits(
        oldoutput, verbose=False)

    # add object spectra to sky model
    iraf.imarith(skymod, "+", galaxy, skymod)

###########
##  start iteration
###########
    ii = 0
    while ii < niter:
          print ii
          # add median of residuals to sky model
          lickshane.util.delete(med5)
          iraf.median(oldoutput, med5, 5, 5, zlor='INDEF',
                      zhir='INDEF', verbose='no')
#          m5 = ndimage.filters.median_filter(_inputarray, size=5, mode='mirror')
          iraf.imarith(skymod, "+", med5, med5)
       
          # take second-order derivative (Laplacian) of input image
          # kernel is convolved with subsampled image, in order to remove negative
          # pattern around high pixels
          lickshane.util.delete(blk)
          lickshane.util.delete(lapla)
          lickshane.util.delete(deriv2)
          lickshane.util.delete(noise)
          lickshane.util.delete(sigmap)
       
          iraf.blkrep(oldoutput, blk, 2, 2)
          iraf.convolve(blk, lapla, '_kernel')
          iraf.imreplace(lapla, 0, upper=0, lower='INDEF')
          iraf.blkavg(lapla, deriv2, 2, 2, option="average")
       
          # create noise model
          iraf.imutil.imexpr(expr='sqrt(a*' + str(gain) + '+' + str(readn) +
                             '**2)/' + str(gain), a=med5, output=noise, verbose='no')
          iraf.imreplace(med5, 0.00001, upper=0, lower='INDEF')
       
          # divide Laplacian by noise model
          iraf.imutil.imexpr(expr='(a/b)/2', a=deriv2, b=noise,
                             output=sigmap, verbose='no')
       
          # removal of large structure (bright, extended objects)
          lickshane.util.delete(med5)
          iraf.median(sigmap, med5, 5, 5, zlo='INDEF', zhi='INDEF', verbose='no')
          iraf.imarith(sigmap, "-", med5, sigmap)
       
          # find all candidate cosmic rays
          # this selection includes sharp features such as stars and HII regions
          arraysigmap, headersigmap = lickshane.cosmics.fromfits(sigmap, verbose=False)
          arrayf = np.where(arraysigmap < sigclip, 0, arraysigmap)
          arrayf = np.where(arrayf > 0.1, 1, arrayf)
          lickshane.cosmics.tofits(firstsel, np.float32(
              arrayf), headersigmap, verbose=False)
       
          # compare candidate CRs to median filtered image
          # this step rejects bright, compact sources from the initial CR list
          # subtract background and smooth component of objects
          lickshane.util.delete(med3)
       
          iraf.median(oldoutput, med3, 3, 3, zlo='INDEF', zhi='INDEF', verbose='no')
          lickshane.util.delete(med7)
          lickshane.util.delete('_' + med3)
          iraf.median(med3, med7, 7, 7, zlo='INDEF', zhi='INDEF', verbose='no')
          iraf.imutil.imexpr(expr='(a-b)/c', a=med3, b=med7,
                             c=noise, output='_' + med3, verbose='no')
          iraf.imreplace('_' + med3, 0.01, upper=0.01, lower='INDEF')
       
       
          # compare CR flux to object flux
          lickshane.util.delete(starreject)
          iraf.imutil.imexpr(expr='(a*b)/c', a=firstsel, b=sigmap,
                             c="_" + med3, output=starreject, verbose='no')
#   ######         #####        ######          #####       ######  FOUND A BUG ? 
#          iraf.imutil.imexpr(expr='a+b+c', a=firstsel, b=sigmap,
#                             c="_" + med3, output=starreject, verbose='no')
       
          # discard if CR flux <= objlim * object flux
          iraf.imreplace(starreject, 0, upper=objlim, lower='INDEF')
          iraf.imreplace(starreject, 1, lower=0.5, upper='INDEF')
          iraf.imarith(firstsel, "*", starreject, firstsel)
       
          # grow CRs by one pixel and check in original sigma map
          arrayfirst, headerfirst = lickshane.cosmics.fromfits(firstsel, verbose=False)
          arraygfirst = lickshane.cosmics.my_convolve_with_FFT2(arrayfirst, gkernel)
       
          arraygfirst = np.where(arraygfirst > 0.5, 1, arraygfirst)
          arraygfirst = arraygfirst * arraysigmap
          arraygfirst = np.where(arraygfirst < sigclip, 0, arraygfirst)
          arraygfirst = np.where(arraygfirst > 0.1, 1, arraygfirst)
       
       
          # grow CRs by one pixel and lower detection limit
          sigcliplow = sigfrac * sigclip
       
          # Finding neighbouring pixels affected by cosmic rays
          arrayfinal = lickshane.cosmics.my_convolve_with_FFT2(arraygfirst, gkernel)
          arrayfinal = np.where(arrayfinal > 0.5, 1, arrayfinal)
          arrayfinal = arrayfinal * arraysigmap
          arrayfinal = np.where(arrayfinal < sigcliplow, 0, arrayfinal)
          arrayfinal = np.where(arrayfinal > 0.1, 1, arrayfinal)
       
          # determine number of CRs found in this iteration
          arraygfirst = (1 - (arrayfinal - arrayfinal)) * arrayfinal
          npix = [str(int(np.size(np.where(arraygfirst > 0.5)) / 2.))]
       
          # create cleaned output image; use 3x3 median with CRs excluded
          arrayoutmask = np.where(arrayfinal > 1, 1, arrayfinal)
          lickshane.cosmics.tofits(outmask, np.float32(
              arrayoutmask), headerfirst, verbose=False)
          lickshane.util.delete(inputmask)
          arrayinputmask = (1 - (10000 * arrayoutmask)) * arrayoldoutput
          lickshane.cosmics.tofits(inputmask, np.float32(
              arrayinputmask), headerfirst, verbose=False)
          lickshane.util.delete(med5)
          iraf.median(inputmask, med5, 5, 5, zloreject=-
                      9999, zhi='INDEF', verbose='no')
          iraf.imarith(outmask, "*", med5, med5)
          lickshane.util.delete(output)
          iraf.imutil.imexpr(expr='(1-a)*b+c', a=outmask, b=oldoutput,
                             c=med5, output=output, verbose='no')

          lickshane.util.delete(oldoutput)
          os.system('cp '+output+' '+oldoutput)

          # add sky and object spectra back in
          iraf.imarith(output, "+", skymod, output)
          # cleanup and get ready for next iteration
          ii = ii + 1
          if npix == 0:
              ii = niter
         # delete temp files

    lickshane.util.delete(blk + "," + lapla + "," + deriv2 + "," + med5)
    lickshane.util.delete(med3 + "," + med7 + "," + noise + "," + sigmap)
    lickshane.util.delete(firstsel + "," + starreject)
    lickshane.util.delete(finalsel + "," + inputmask)
    lickshane.util.delete(oldoutput + "," + skymod + "," + galaxy)
    lickshane.util.delete("_" + med3 + ",_" + sigmap)
    lickshane.util.delete('_kernel' + "," + '_gkernel')
    lickshane.util.delete(outmask)

#################################################


def clean_image(img, cleanimg):
    # print "LOGX:: Entering `clean_image` method/function in %(__file__)s" %
    # globals()
    import lickshane
    from lickshane.util import readkey3, readhdr, delete
    array, header = lickshane.cosmics.fromfits(img, verbose=False)
    import warnings

    def fxn():
        # print "LOGX:: Entering `fxn` method/function in %(__file__)s" %
        # globals()
        warnings.warn(" ", DeprecationWarning)

    original_filters = warnings.filters[:]
    # Ignore warnings.
    warnings.simplefilter("ignore")
    try:
        c = lickshane.cosmics.cosmicsimage(array, gain=readkey3(header, 'gain'), readnoise=readkey3(
            header, 'ron'), sigclip=5.0, sigfrac=0.3, objlim=5.0, verbose=False)
        c.run(maxiter=4, verbose=False)
        fxn()
    finally:
        warnings.filters = original_filters

    if not cleanimg:
        delete(img)
        cleanimg = img
    lickshane.cosmics.tofits(cleanimg, c.cleanarray, header, verbose=False)
    return cleanimg

##################################################


def my_convolve_with_FFT2(image1, kernel):
    # print "LOGX:: Entering `my_convolve_with_FFT2` method/function in
    # %(__file__)s" % globals()
    import numpy as np
    r1, c1 = image1.shape
    r2, c2 = kernel.shape
    r = r1 + r2 - 1
    c = c1 + c2 - 1
    padr = (r - r1) / 2  # add to each side
    padc = (c - c1) / 2  # add to each of top and bottom

    # pad the edges, with the same values found on the edges
    lside = np.empty((image1.shape[0], padc))
    rside = np.empty((image1.shape[0], padc))
    for i in range(padc):
        lside[:, i] = np.copy(image1[:, 0]).T
        rside[:, i] = np.copy(image1[:, image1.shape[1] - 1]).T
    # end for loop

    image1 = np.hstack((lside, image1, rside))
    top = np.empty((padr, image1.shape[1]))
    bot = np.empty((padr, image1.shape[1]))
    # pad the top and bottom
    for i in range(padr):
        top[i] = np.copy(image1[0, :]).T
        bot[i] = np.copy(image1[image1.shape[0] - 1, :]).T
    # end for loop
    image1 = np.vstack((top, image1, bot))
    rOrig = r
    cOrig = c
    pr2 = int(np.log(r) / np.log(2.0) + 1.0)
    pc2 = int(np.log(c) / np.log(2.0) + 1.0)
    r = 2**pr2
    c = 2**pc2
    fftimage = np.fft.fft2(image1, s=(r, c)) * \
        np.fft.fft2(kernel[::-1, ::-1], s=(r, c))
    ret = np.fft.ifft2(fftimage).real
    return ret[(rOrig - r1):rOrig, (cOrig - c1):cOrig]

###################################################
