from astropy.io import fits as pyfits
import numpy as np

def fluxcalib2d(img2d, sensfun):  # flux calibrate 2d images
    # print "LOGX:: Entering `fluxcalib2d` method/function in %(__file__)s" %
    # globals()
    import re
    import string
    import lickshane

    data2d, hdr2d = pyfits.getdata(img2d, 0, header=True)
    xxd = np.arange(len(data2d[:, 0]))
    # aad=crvald+(xxd)*cdd
    # crvald=readkey3(hdr2d,'CRVAL2')
    # cdd=readkey3(hdr2d,'CD2_2')
    crvald = pyfits.open(img2d)[0].header.get('CRVAL2')
    cdd = pyfits.open(img2d)[0].header.get('CD2_2')
    _exptime = lickshane.util.readkey3(lickshane.util.readhdr(img2d), 'exptime')
    _airmass = lickshane.util.readkey3(lickshane.util.readhdr(img2d), 'airmass')
    #  read sensfunction and interpole pixel of 2D image
    yys = pyfits.open(sensfun)[0].data
    crvals = pyfits.open(sensfun)[0].header.get('CRVAL1')
    cds = pyfits.open(sensfun)[0].header.get('CD1_1')
    yys = (10 ** (yys / 2.5)) * cds  # from sens iraf in sens flux
    xxs = np.arange(len(yys))
    aasens = crvals + (xxs) * cds
    xxs2 = (aasens - crvald) / cdd
    aasens2 = np.interp(xxd, xxs2, yys)
    #  read atmosferic function and interpole pixel of 2D image
    aae, yye = lickshane.util.ReadAscii2(
        lickshane.__path__[0] + '/standard/extinction/lick.dat')
    aae, yye = np.array(aae, float), np.array(yye, float)
    xxe = (aae - crvald) / cdd
    atm_xx = np.interp(xxd, xxe, yye)
    aircorr = 10 ** (0.4 * np.array(atm_xx) * _airmass)
    img2df = re.sub('.fits', '_2df.fits', img2d)
    for i in range(len(data2d[0, :])):
        data2d[:, i] = ((np.array(data2d[:, i] / _exptime) *
                         np.array(aircorr)) / aasens2) * 1e20
    lickshane.util.delete(img2df)
    pyfits.writeto(img2df, np.float32(data2d), hdr2d)
    lickshane.util.updateheader(
        img2df, 0, {'SENSFUN': [string.split(sensfun, '/')[-1], '']})
    lickshane.util.updateheader(img2df, 0, {
                          'BUNIT': ['10^20 erg/cm2/s/Angstrom', 'Physical unit of array values']})
    return img2df



def telluric_atmo(imgstd):
    # print "LOGX:: Entering `telluric_atmo` method/function in %(__file__)s"
    # % globals()

    import lickshane
    from pyraf import iraf

    iraf.images(_doprint=0)
    iraf.noao(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.onedspec(_doprint=0)
    toforget = ['imfilter.gauss', 'specred.apall', 'longslit.identify', 'longslit.reidentify', 'specred.standard',
                'onedspec.wspectext']
    for t in toforget:
        iraf.unlearn(t)

    _grism = lickshane.util.readkey3(lickshane.util.readhdr(imgstd), 'grism')

    imgout = 'invers_atmo_' + imgstd
    lickshane.util.delete(imgout)
    iraf.set(direc=lickshane.__path__[0] + '/')
    _cursor = 'direc$standard/ident/cursor_sky_0'
    iraf.noao.onedspec.bplot(imgstd, cursor=_cursor,
                             spec2=imgstd, new_ima=imgout, overwri='yes')
    xxstd, ffstd = lickshane.util.readspectrum(imgout)
    if _grism in ['300_7500_3000','300_7500_5000','300_7500_6000']:
        llo2 = np.compress((np.array(xxstd) >= 7550) & (
            np.array(xxstd) <= 7750), np.array(xxstd))
        llh2o = np.compress((np.array(xxstd) >= 7100) & (
            np.array(xxstd) <= 7500), np.array(xxstd))
        ffo2 = np.compress((np.array(xxstd) >= 7550) & (
            np.array(xxstd) <= 7750), np.array(ffstd))
        ffh2o = np.compress((np.array(xxstd) >= 7100) & (
            np.array(xxstd) <= 7500), np.array(ffstd))
    elif _grism in ['600_4310']:
        llo2 = np.compress((np.array(xxstd) >= 6830) & (
            np.array(xxstd) <= 7100), np.array(xxstd))
        llh2o = np.compress((np.array(xxstd) >= 7100) & (
            np.array(xxstd) <= 7500), np.array(xxstd))
        ffo2 = np.compress((np.array(xxstd) >= 6830) & (
            np.array(xxstd) <= 7100), np.array(ffstd))
        ffh2o = np.compress((np.array(xxstd) >= 7100) & (
            np.array(xxstd) <= 7500), np.array(ffstd))
    if _grism in ['600_4310', '300_7500_3000','300_7500_5000','300_7500_6000']:
        _skyfileh2o = 'direc$standard/ident/ATLAS_H2O.fits'
        _skyfileo2 = 'direc$standard/ident/ATLAS_O2.fits'
        atlas_smooto2 = '_atlas_smoot_o2.fits'
        atlas_smooth2o = '_atlas_smoot_h2o.fits'
        _sigma = 200
        lickshane.util.delete(atlas_smooto2)
        lickshane.util.delete(atlas_smooth2o)
        iraf.imfilter.gauss(_skyfileh2o, output=atlas_smooth2o, sigma=_sigma)
        iraf.imfilter.gauss(_skyfileo2, output=atlas_smooto2, sigma=_sigma)
        llskyh2o, ffskyh2o = lickshane.util.readspectrum(atlas_smooth2o)
        llskyo2, ffskyo2 = lickshane.util.readspectrum(atlas_smooto2)
        ffskyo2cut = np.interp(llo2, llskyo2, ffskyo2)
        ffskyh2ocut = np.interp(llh2o, llskyh2o, ffskyh2o)
        _scaleh2o = []
        integral_h2o = []
        for i in range(1, 21):
            j = 0.6 + i * 0.04
            _ffskyh2ocut = list((np.array(ffskyh2ocut) * j) + 1 - j)
            diff_h2o = abs(_ffskyh2ocut - ffh2o)
            integraleh2o = np.trapz(diff_h2o, llh2o)
            integral_h2o.append(integraleh2o)
            _scaleh2o.append(j)
        _scaleo2 = []
        integral_o2 = []
        for i in range(1, 21):
            j = 0.6 + i * 0.04
            _ffskyo2cut = list((np.array(ffskyo2cut) * j) + 1 - j)
            diff_o2 = abs(_ffskyo2cut - ffo2)
            integraleo2 = np.trapz(diff_o2, llo2)
            integral_o2.append(integraleo2)
            _scaleo2.append(j)
        sh2o = _scaleh2o[np.argmin(integral_h2o)]
        so2 = _scaleo2[np.argmin(integral_o2)]
        telluric_features = ((np.array(ffskyh2o) * sh2o) +
                             1 - sh2o) + ((np.array(ffskyo2) * so2) + 1 - so2) - 1
        telluric_features = np.array([1] + list(telluric_features) + [1])
        llskyo2 = np.array([1000] + list(llskyo2) + [15000])
        telluric_features_cut = np.interp(xxstd, llskyo2, telluric_features)
        _imgout = 'atmo_' + imgstd

        data1, hdr = pyfits.getdata(imgstd, 0, header=True)
        data1[0] = np.array(telluric_features_cut)
        data1[1] = data1[1] / data1[1]
        data1[2] = data1[2] / data1[2]
        data1[3] = data1[3] / data1[3]
        lickshane.util.delete(_imgout)
        pyfits.writeto(_imgout, np.float32(data1), hdr)
        lickshane.util.delete(atlas_smooto2)
        lickshane.util.delete(atlas_smooth2o)
        lickshane.util.delete(imgout)
    else:
        _imgout = ''
        print '### telluric correction with model not possible '
    return _imgout


def checkwavestd(imgex, _interactive):
    # print "LOGX:: Entering `checkwavestd` method/function in %(__file__)s" %
    # globals()
    import lickshane
    import numpy as np

    print '\n### Warning: check in wavelenght with sky lines not performed\n'
    if _interactive in ['yes', 'YES', 'Yes', 'Y', 'y']:
        answ = raw_input(
            '\n### Do you want to check the wavelengh calibration with tellurich lines [[y]/n]? ')
        if not answ:
            answ = 'y'
    else:
        answ = 'y'
    if answ in ['y', 'yes']:
        print '\n### check wavelength calibration with tellurich lines \n'
        _skyfile = lickshane.__path__[0] + '/standard/ident/sky_new_0.fits'
        skyff = 1 - (pyfits.open(_skyfile)[0].data)
        crval1 = pyfits.open(_skyfile)[0].header.get('CRVAL1')
        cd1 = pyfits.open(_skyfile)[0].header.get('CD1_1')
        skyxx = np.arange(len(skyff))
        skyaa = crval1 + (skyxx) * cd1
        atmofile = lickshane.lickshane1Ddef.atmofile(imgex, 'atmo2_' + imgex)
        atmoff = 1 - (pyfits.open(atmofile)[0].data[0][0])
        crval1 = pyfits.open(atmofile)[0].header.get('CRVAL1')
        cd1 = pyfits.open(atmofile)[0].header.get('CD1_1')
        atmoxx = np.arange(len(atmoff))
        atmoaa = crval1 + (atmoxx) * cd1
        shift = lickshane.lickshane2Ddef.checkwavelength_arc(
            atmoaa, atmoff, skyaa, skyff, '', '')
    else:
        shift = 0
    zro = pyfits.open(imgex)[0].header.get('CRVAL1')
    if _interactive in ['yes', 'YES', 'Yes', 'Y', 'y']:
        answ = raw_input(
            '\n### do you want to correct the wavelengh calibration with this shift: ' + str(shift) + ' [[y]/n] ? ')
        if not answ:
            answ = 'y'
        if answ.lower() in ['y', 'yes']:
            lickshane.util.updateheader(imgex, 0, {'CRVAL1': [zro + int(shift), '']})
            lickshane.util.updateheader(imgex, 0, {'shift': [float(shift), '']})
    else:
        lickshane.util.updateheader(imgex, 0, {'CRVAL1': [zro + int(shift), '']})
        lickshane.util.updateheader(imgex, 0, {'shift': [float(shift), '']})


# ###################################

def atmofile(imgstd, imgout=''):
    # print "LOGX:: Entering `atmofile` method/function in %(__file__)s" %
    # globals()
    from pyraf import iraf
    import os
    import lickshane

    iraf.noao(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.set(direc=lickshane.__path__[0] + '/')
    _cursor = 'direc$standard/ident/cursor_sky_0'
    if not imgout:
        imgout = 'atmo_' + imgstd
    os.system('rm -rf ' + imgout)
    iraf.noao.onedspec.bplot(imgstd, cursor=_cursor,
                             spec2=imgstd, new_ima=imgout, overwri='yes')
    return imgout


def sensfunction(standardfile, _function, _order, _interactive):
    # print "LOGX:: Entering `sensfunction` method/function in %(__file__)s" %
    # globals()
    import re
    import os
    import sys
    import lickshane
    import datetime
    from pyraf import iraf

    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.scopy', 'specred.sensfunc', 'specred.standard']
    for t in toforget:
        iraf.unlearn(t)
    iraf.specred.scopy.format = 'multispec'
    iraf.specred.verbose = 'no'
    hdrs = lickshane.util.readhdr(standardfile)
    try:
        _outputsens = 'sens_' + str(lickshane.util.readkey3(hdrs, 'date-night')) + '_' + \
                      str(lickshane.util.readkey3(hdrs, 'grism')) + '_' + str(lickshane.util.readkey3(hdrs, 'filter')) + '_' + \
                      re.sub('.dat', '', lickshane.util.readkey3(
                          hdrs, 'stdname')) + '_' + str(MJDtoday)
    except:
        sys.exit('Error: missing header -stdname- in standard ' +
                 str(standardfile) + '  ')

    _outputsens = lickshane.util.name_duplicate(standardfile, _outputsens, '')
    if os.path.isfile(_outputsens):
        if _interactive.lower() != 'yes':
            lickshane.util.delete(_outputsens)
        else:
            answ = raw_input(
                'sensitivity function already computed, do you want to do it again [[y]/n] ? ')
            if not answ:
                answ = 'y'
            if answ.lower() in ['y', 'yes']:
                lickshane.util.delete(_outputsens)

    if not os.path.isfile(_outputsens):
        iraf.set(direc=lickshane.__path__[0] + '/')
        _caldir = 'direc$standard/MAB/'
        _extinctdir = 'direc$standard/extinction/'
        _observatory = 'lick'
        _extinction = 'lick.dat'
        refstar = 'm' + \
            re.sub('.dat', '', pyfits.open(standardfile)
                   [0].header.get('stdname'))
        _airmass = lickshane.util.readkey3(hdrs, 'airmass')
        _exptime = lickshane.util.readkey3(hdrs, 'exptime')
        _outputstd = 'std_' + str(lickshane.util.readkey3(hdrs, 'grism')) + '.fits'
        lickshane.util.delete(_outputstd)
        lickshane.util.delete(_outputsens)
        iraf.specred.standard(input=standardfile, output=_outputstd, extinct=_extinctdir + _extinction,
                              caldir=_caldir, observa=_observatory, star_nam=refstar, airmass=_airmass,
                              exptime=_exptime, interac=_interactive)
        iraf.specred.sensfunc(standard=_outputstd, sensitiv=_outputsens, extinct=_extinctdir + _extinction,
                              ignorea='yes', observa=_observatory, functio=_function, order=_order,
                              interac=_interactive)

        data, hdr = pyfits.getdata(standardfile, 0, header=True)  # added later
        data1, hdr1 = pyfits.getdata(
            _outputsens, 0, header=True)  # added later
        lickshane.util.delete(_outputsens)  # added later
        pyfits.writeto(_outputsens, np.float32(data1), hdr)  # added later
    return _outputsens

###########################################################################################################

##################################################

def extractspectrum(img, dv, _ext_trace, _dispersionline, _interactive, _type, automaticex=False):
    # print "LOGX:: Entering `extractspectrum` method/function in
    # %(__file__)s" % globals()
    import glob
    import os
    import string
    import sys
    import re
    import lickshane
    import datetime

    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    from pyraf import iraf
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.apall', 'specred.transform']
    for t in toforget:
        iraf.unlearn(t)

    dv = lickshane.util.dvex()
    hdr = lickshane.util.readhdr(img)
    _gain = lickshane.util.readkey3(hdr, 'gain')
    _rdnoise = lickshane.util.readkey3(hdr, 'ron')
    _grism = lickshane.util.readkey3(hdr, 'grism')

    _instrument = lickshane.util.readkey3(hdr, 'version')

    imgex = re.sub('.fits', '_ex.fits', img)
    imgfast = re.sub(string.split(img, '_')[-2] + '_', '', img)
    # imgfast=re.sub(str(MJDtoday)+'_','',img)
    if not os.path.isfile(imgex) and not os.path.isfile(
            'database/ap' + re.sub('.fits', '', img)) and not os.path.isfile(
            'database/ap' + re.sub('.fits', '', imgfast)):
        _new = 'yes'
        _extract = 'yes'
    else:
        if automaticex:
            if _interactive in ['Yes', 'yes', 'YES', 'y', 'Y']:
                answ = 'x'
                while answ not in ['o', 'n', 's']:
                    answ = raw_input(
                        '\n### New extraction [n], extraction with old parameters [o], skip extraction [s] ? [o]')
                    if not answ:
                        answ = 'o'
                if answ == 'o':
                    _new, _extract = 'no', 'yes'
                elif answ == 'n':
                    _new, _extract = 'yes', 'yes'
                else:
                    _new, _extract = 'yes', 'no'
            else:
                _new, _extract = 'no', 'yes'
        else:
            if _interactive in ['Yes', 'yes', 'YES', 'y', 'Y']:
                answ = 'x'
                while answ not in ['y', 'n']:
                    answ = raw_input(
                        '\n### do you want to extract again [[y]/n] ? ')
                    if not answ:
                        answ = 'y'
                if answ == 'y':
                    _new, _extract = 'yes', 'yes'
                else:
                    _new, _extract = 'yes', 'no'
            else:
                _new, _extract = 'yes', 'yes'
    if _extract == 'yes':
        lickshane.util.delete(imgex)
        if _dispersionline:
            question = 'yes'
            while question == 'yes':
                _z1, _z2, goon = lickshane.util.display_image(img, 1, '', '', False)
                dist = raw_input(
                    '\n### At which line do you want to extract the spectrum [' + str(dv['line'][_grism]) + '] ? ')
                if not dist:
                    dist = 400
                try:
                    dist = int(dist)
                    question = 'no'
                except:
                    print '\n### input not valid, try again:'
        else:
            dist = dv['line'][_grism]
        if _ext_trace in ['yes', 'Yes', 'YES', True]:
            lista = glob.glob('*ex.fits')
            if lista:
                for ii in lista:
                    print ii
                _reference = raw_input(
                    '\### which object do you want to use for the trace [' + str(lista[0]) + '] ? ')
                if not _reference:
                    _reference = lista[0]
                _reference = re.sub('_ex', '', _reference)
                _fittrac = 'no'
                _trace = 'no'
            else:
                sys.exit('\n### error: no extracted spectra in the directory')
        else:
            _reference = ''
            _fittrac = 'yes'
            _trace = 'yes'
        if _new == 'no':
            if not os.path.isfile('database/ap' + re.sub('.fits', '', img)):
                lickshane.util.repstringinfile('database/ap' + re.sub('.fits', '', imgfast),
                                         'database/ap' +
                                         re.sub('.fits', '', img), re.sub(
                                             '.fits', '', imgfast),
                                         re.sub('.fits', '', img))
            _find = 'no'
            _recenter = 'no'
            _edit = 'no'
            _trace = 'no'
            _fittrac = 'no'
            _mode = 'h'
            _resize = 'no'
            _review = 'no'
            iraf.specred.mode = 'h'
            _interactive = 'no'
        else:
            iraf.specred.mode = 'q'
            _mode = 'q'
            _find = 'yes'
            _recenter = 'yes'
            _edit = 'yes'
            _review = 'yes'
            _resize = dv[_type]['_resize']

        if _instrument == 'kastb':
            iraf.specred.dispaxi = 1
        elif _instrument == 'kastr':
            iraf.specred.dispaxi = 2

        iraf.specred.apall(img, output=imgex, referen=_reference, trace=_trace, fittrac=_fittrac, find=_find,
                           recenter=_recenter, edit=_edit,
                           nfind=1, extract='yes', backgro='fit', gain=_gain, readnoi=_rdnoise, lsigma=4, usigma=4,
                           format='multispec',
                           b_function='legendre', b_sample=dv[_type]['_b_sample'], clean='yes', pfit='fit1d',
                           lower=dv[_type]['_lower'], upper=dv[_type][
                               '_upper'], t_niter=dv[_type]['_t_niter'],
                           width=dv[_type]['_width'],
                           radius=dv[_type]['_radius'], line=dist, nsum=dv[
                               _type]['_nsum'], t_step=dv[_type]['_t_step'],
                           t_nsum=dv[_type]['_t_nsum'],
                           t_nlost=dv[_type]['_t_nlost'], t_sample=dv[
                               _type]['_t_sample'], resize=_resize,
                           t_order=dv[_type]['_t_order'],
                           weights=dv[_type]['_weights'], interactive=_interactive, review=_review, mode=_mode)

        lickshane.util.repstringinfile('database/ap' + re.sub('.fits', '', img), 'database/ap' + re.sub('.fits', '', imgfast),
                                 re.sub('.fits', '', img), re.sub('.fits', '', imgfast))
    else:
        print '\n### skipping new extraction'
    return imgex



###########################################################################################################

def lickshane1Dredu(files, _interactive,  _ext_trace, _dispersionline, _automaticex, _verbose=False):
    import lickshane
    import datetime
    import os
    import re
    import string
    import sys

    liststandard = ''
    listatmo0 = ''

#    os.environ["PYRAF_BETA_STATUS"] = "1"
    _extinctdir = 'direc$standard/extinction/'
    _extinction = 'lick.dat'
    _observatory = 'lick'

    now = datetime.datetime.now()
    datenow = now.strftime('20%y%m%d%H%M')
    MJDtoday = 55927 + (datetime.date.today() - datetime.date(2012, 01, 01)).days
    scal = np.pi / 180.
    dv = lickshane.util.dvex()
    std, rastd, decstd, magstd = lickshane.util.readstandard('standard_lick_mab.txt')
    objectlist = {}
    for img in files:
        hdr = lickshane.util.readhdr(img)
        img = re.sub('\n', '', img)
        lickshane.util.correctcard(img)
        _ra = lickshane.util.readkey3(hdr, 'RA')
        _dec = lickshane.util.readkey3(hdr, 'DEC')
        _object = lickshane.util.readkey3(hdr, 'object')
        _grism = lickshane.util.readkey3(hdr, 'grism')
        _slit = lickshane.util.readkey3(hdr, 'slit')
        dd = np.arccos(np.sin(_dec * scal) * np.sin(decstd * scal) + np.cos(_dec * scal) *
                       np.cos(decstd * scal) * np.cos((_ra - rastd) * scal)) * ((180 / np.pi) * 3600)
        if min(dd) < 100:
            _type = 'stdsens'
        else:
            _type = 'obj'
        print img,_type
        if min(dd) < 100:
            lickshane.util.updateheader(
                img, 0, {'stdname': [std[np.argmin(dd)], '']})
            lickshane.util.updateheader(
                img, 0, {'magstd': [float(magstd[np.argmin(dd)]), '']})
        if _type not in objectlist:
            objectlist[_type] = {}
        if (_grism,  _slit) not in objectlist[_type]:
            objectlist[_type][_grism, _slit] = [img]
        else:
            objectlist[_type][_grism, _slit].append(img)

    from pyraf import iraf
    iraf.set(stdimage='imt2048')
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    iraf.imutil(_doprint=0)
    toforget = ['imutil.imcopy', 'specred.sarith', 'specred.standard']
    for t in toforget:
        iraf.unlearn(t)
    iraf.specred.verbose = 'no'
    iraf.specred.dispaxi = 2
    iraf.set(direc=lickshane.__path__[0] + '/')
    sens = {}
    print objectlist
    outputfile = []
    if 'obj' in objectlist.keys():
        tpe = 'obj'
    elif 'stdsens' in objectlist.keys():
        tpe = 'stdsens'
    else:
        sys.exit('error: no objects and no standards in the list')

    for setup in objectlist[tpe]:
        extracted = []
        listatmo = []
        if setup not in sens:
            sens[setup] = []
        if tpe == 'obj':
            print '\n### setup= ', setup, '\n### objects= ', objectlist['obj'][setup], '\n'
            for img in objectlist['obj'][setup]:
                #              hdr=readhdr(img)
                print '\n\n### next object= ', img, ' ', lickshane.util.readkey3(lickshane.util.readhdr(img), 'object'), '\n'
                
                #_automaticex = ''

                imgex = lickshane.lickshane1Ddef.extractspectrum(img, dv, _ext_trace, _dispersionline, _interactive, 'obj',
                                                 automaticex=_automaticex)
                if not os.path.isfile(imgex):
                    sys.exit('### error, extraction not computed')
                if not lickshane.util.readkey3(lickshane.util.readhdr(imgex), 'shift') and \
                        lickshane.util.readkey3(lickshane.util.readhdr(imgex), 'shift') != 0.0:
                    if setup in ['300_7500']:
                        lickshane.lickshane1Ddef.checkwavestd(imgex, _interactive)
                    else:
                        print 'wave check using teluric not possible'

                extracted.append(imgex)
                if imgex not in outputfile:
                    outputfile.append(imgex)
                lickshane.util.updateheader(
                    imgex, 0, {'TRACE1': [img, 'Originating file']})
                if os.path.isfile('database/ap' + re.sub('_ex.fits', '', imgex)):
                    if 'database/ap' + re.sub('_ex.fits', '', imgex) not in outputfile:
                        outputfile.append(
                            'database/ap' + re.sub('_ex.fits', '', imgex))
            print '\n### all object with this setup extracted\n'

        if liststandard:
            standardlist = liststandard
            _type = 'stdfromdreducer'
        else:
            try:
                standardlist = objectlist['stdsens'][setup]
                _type = 'stdsens'
            except:
                standardlist = ''
                _type = ''
        if _type == 'stdfromdreducer' and len(extracted) >= 1:
            _outputsens2 = lichshane.util.searchsens(extracted[0], standardlist)[0]
            print '\n### using standard from reducer ' + str(_outputsens2)
        elif _type not in ['stdsens', 'stdfromdreducer'] and len(extracted) >= 1:
            _outputsens2 = lickshane.util.searchsens(extracted[0], '')[0]
            os.system('cp ' + _outputsens2 + ' .')
            _outputsens2 = string.split(_outputsens2, '/')[-1]
            print '\n### no standard in the list, using standard from archive'
        else:
            for simg in standardlist:
                print '\n###  standard for setup ' + \
                      str(setup) + ' = ', simg, ' ', lickshane.util.readkey3(
                          lickshane.util.readhdr(simg), 'object'), '\n'
                simgex = lickshane.lickshane1Ddef.extractspectrum(
                    simg, dv, False, False, _interactive, 'std', automaticex=_automaticex)
                lickshane.util.updateheader(
                    simgex, 0, {'TRACE1': [simg, 'Originating file']})
                if not lickshane.util.readkey3(lickshane.util.readhdr(simgex), 'shift') and lickshane.util.readkey3(lickshane.util.readhdr(simgex), 'shift') != 0.0:
                    lickshane.lickshane1Ddef.checkwavestd(simgex, _interactive)
                print simgex
                atmofile = lickshane.lickshane1Ddef.telluric_atmo(simgex)  # atmo file2
                print atmofile
                lickshane.util.updateheader(atmofile, 0, {'TRACE1': [simgex, 'Originating file']})
                if tpe != 'obj' and atmofile not in outputfile:
                    outputfile.append(atmofile)
                if not listatmo0:
                    listatmo.append(atmofile)
                sens[setup].append(simgex)
                if simgex not in outputfile:
                    outputfile.append(simgex)
            print '\n### standard available: ', sens[setup]

            if tpe == 'obj':
                if len(sens[setup]) > 1:
                    goon = 'no'
                    while goon != 'yes':
                        stdused = raw_input(
                            '\n### more than one standard for this setup, which one do you want to use [' + sens[setup][
                                0] + '] ?')
                        if not stdused:
                            stdused = sens[setup][0]
                        if os.path.isfile(stdused):
                            goon = 'yes'
                else:
                    stdused = sens[setup][0]
                stdvec = [stdused]
            else:
                stdvec = sens[setup]
            for stdused in stdvec:
                stdusedclean = re.sub('_ex', '_clean', stdused)
                lickshane.util.delete(stdusedclean)
                iraf.specred.sarith(
                    input1=stdused, op='/', input2=atmofile, output=stdusedclean, format='multispec')
                _outputsens2 = lickshane.lickshane1Ddef.sensfunction(stdusedclean, 'spline3', 16, _interactive)
                lickshane.util.updateheader(_outputsens2, 0, {'TRACE1': [stdused, 'Originating file']})

                if _outputsens2 not in outputfile:
                    outputfile.append(_outputsens2)
        if _outputsens2 and tpe == 'obj':
            ####################################################
            for img in objectlist['obj'][setup]:  # flux calibrate 2d images
                imgd = fluxcalib2d(img, _outputsens2)
                lickshane.util.updateheader(
                    imgd, 0, {'TRACE1': [img, 'Originating files']})
                if imgd not in outputfile:
                    outputfile.append(imgd)
            ####################################################
            #    flux calib in the standard way
            if not listatmo and listatmo0:
                listatmo = listatmo0[:]
            for _imgex in extracted:
                _airmass = lickshane.util.readkey3(
                    lickshane.util.readhdr(_imgex), 'airmass')
                _exptime = lickshane.util.readkey3(
                    lickshane.util.readhdr(_imgex), 'exptime')
                _imgf = re.sub('_ex.fits', '_f.fits', _imgex)
                lickshane.util.delete(_imgf)
                qqq = iraf.specred.calibrate(input=_imgex, output=_imgf, sensiti=_outputsens2, extinct='yes',
                                             flux='yes',
                                             extinction=_extinctdir + _extinction, observatory=_observatory,
                                             airmass=_airmass, ignorea='yes', exptime=_exptime, fnu='no')
                hedvec = {'SENSFUN': [_outputsens2, ''],
#                          'SNR': [lickshane.util.StoN2(_imgf, False), 'Average signal to noise ratio per pixel'],
                          'BUNIT': ['erg/cm2/s/Angstrom', 'Physical unit of array values'],
                          'TRACE1': [_imgex, 'Originating file'],
                          'ASSON1': [re.sub('_f.fits', '_2df.fits', _imgf), 'Name of associated file'],
                          'ASSOC1': ['ANCILLARY.2DSPECTRUM', 'Category of associated file']}
                lickshane.util.updateheader(_imgf, 0, hedvec)
                if _imgf not in outputfile:
                    outputfile.append(_imgf)
                if listatmo:
                    atmofile = lickshane.util.searcharc(_imgex, listatmo)[0]
                    if atmofile:
                        _imge = re.sub('_f.fits', '_e.fits', _imgf)
                        lickshane.util.delete(_imge)
                        iraf.specred.sarith(input1=_imgf, op='/', input2=atmofile, output=_imge, w1='INDEF', w2='INDEF',
                                            format='multispec')
                        try:
                            iraf.imutil.imcopy(
                                input=_imgf + '[*,1,2]', output=_imge + '[*,1,2]', verbose='no')
                        except:
                            pass
                        try:
                            iraf.imutil.imcopy(
                                input=_imgf + '[*,1,3]', output=_imge + '[*,1,3]', verbose='no')
                        except:
                            pass
                        try:
                            iraf.imutil.imcopy(
                                input=_imgf + '[*,1,4]', output=_imge + '[*,1,4]', verbose='no')
                        except:
                            pass
                        if _imge not in outputfile:
                            outputfile.append(_imge)
                        if atmofile not in outputfile:
                            outputfile.append(atmofile)
                        lickshane.util.updateheader(
                            _imge, 0, {'ATMOFILE': [atmofile, '']})
                        lickshane.util.updateheader(
                            _imge, 0, {'TRACE1': [_imgf, 'Originating file']})
                        imgin = _imge
                    else:
                        imgin = _imgf
                else:
                    imgin = _imgf
                imgasci = re.sub('.fits', '.asci', imgin)

                lickshane.util.delete(imgasci)
                iraf.onedspec(_doprint=0)
                iraf.onedspec.wspectext(
                    imgin + '[*,1,1]', imgasci, header='no')
                if imgasci not in outputfile:
                    outputfile.append(imgasci)
            


    return objectlist,'ddd'

