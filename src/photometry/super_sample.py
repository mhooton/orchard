# -*- coding: utf-8 -*-
"""
Usage:
    super_sample.py [options] (<FILE_LIST>) (<OUT_NAME>)

Options:
    -h --help  Show help text
    --factor=FACTOR  What oversampling factor to use [default: 5]
    --size=SIZE  how large a region around each star to stack [default: 11]
    --stars=STARS  How many stars to stack in each quadrant [default: 100]
    --binning=BINNING  use a binning factor for long time series? [default: 1]

"""

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import pickle
import threading
import scipy.optimize as opt
import os
import numpy as np
import multiprocessing.dummy as multithreading
import multiprocessing
from functools import partial
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in divide')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in sqrt')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='Number of calls to function has reached maxfev.*')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='divide by zero encountered in divide')

matplotlib.rc('text', usetex=False)

def super_sample(filelist,inputcat,factor,size,stars,binning,tag,nproc=4):

    p = multiprocessing.Pool(nproc)

    # make curry

    files = []
    files_model = ''
    files_residuals = ''
    files_psf = ''

    for line in open(filelist,'r'):
        files += [line.rstrip('\n')]
        files_psf += files[-1].rstrip('.fits')+'_psf.png,'
        files_model += files[-1].rstrip('.fits')+'_model.png,'
        files_residuals += files[-1].rstrip('.fits')+'_residuals.png,'
    
    data_points = len(files)

    f_5x = []
    f_5y = []
    theta = []
    mjd = []

    curry = []
    for i in range(0,data_points):
        curry += [[files[i],inputcat,factor,size,stars,tag]]
        with pf.open(files[i]) as imdata:
            mjd +=[imdata[0].header['MJD']]

    output = p.map(uncurry_call_find_fwhm,curry)

    for o in output:
        f_5x += o[0]['f_5']
        f_5y += o[1]['f_5']
        theta += o[2]['f_5']

    plt.plot(mjd,f_5x,'bo')
    plt.plot(mjd,f_5y,'ro')

    plt.savefig(tag+'_xy.png', bbox_inches=0)

    plt.close()

    plt.plot(mjd,theta,'ro')

    plt.savefig(tag+'_theta.png', bbox_inches=0)

    plt.close()

    fps = 5

    outputf = tag+'_psf.avi'
    command = 'mencoder mf://'+files_psf.rstrip(',')+' -mf w=800:h=600:fps='+str(fps)+':type=png -ovc raw -oac copy -o '+outputf
    print(command)
    os.system(command)

    outputf = tag+'_residuals.avi'
    command = 'mencoder mf://'+files_residuals.rstrip(',')+' -mf w=800:h=600:fps='+str(fps)+':type=png -ovc raw -oac copy -o '+outputf
    print(command)
    os.system(command)

    outputf = tag+'_model.avi'
    command = 'mencoder mf://'+files_model.rstrip(',')+' -mf w=800:h=600:fps='+str(fps)+':type=png -ovc raw -oac copy -o '+outputf
    print(command)
    os.system(command)

    quit()

    i = 0
    for file in files:
        fwhm_a, fwhm_b, t = call_find_fwhm(file,inputcat,factor,size,stars,tag=tag)
        f_5x += [fwhm_a['f_5']]
        f_5y += [fwhm_b['f_5']]
        theta += [t['f_5']]
        i +=1

 # labels = ['f_1','f_3','f_5','f_7','f_9']

#  for label in labels:
#        condense_data(label,tag=tag)

#  plot_everything(files[0:10],labels,binning,tag=tag)

def uncurry_call_find_fwhm(c):
    return call_find_fwhm(c[0],c[1],c[2],c[3],c[4],tag=c[5])

def condense_data(label,tag=''):

    fwhm = []
    ellipse = []
    theta = []
    x = []
    y = []


    theta += pickle.load(open(tag+'theta'+label+'.p'))
    x += pickle.load(open(tag+'fwhm_a'+label+'.p'))
    y += pickle.load(open(tag+'fwhm_b'+label+'.p'))

    for i in range(0,len(x)):
        rad = theta[i]*np.pi/180.0
        theta[i] = theta[i]*np.pi/180.0
        theta[i] = theta[i] - 2*np.pi*(int(theta[i]/(2*np.pi)))
        # print theta[i]

        if x[i] < y[i]:
            theta[i] -= np.pi/2.0

        if theta[i] < 0:
            theta[i] += 2*np.pi
        if theta[i] > np.pi:
            theta[i] -= np.pi

        theta[i] = theta[i]*180.0/np.pi

    #  print 'major?',abs(y[i]/np.cos(rad))

        ellipse_frame, average = find_minor_axis(rad,x[i],y[i])

#        A = 1.0

#        xp = arange(0,220)
#        yp = arange(0,220)

#        x0 = 110
#        y0 = 110

#        f = gaussian2d(A,theta[i],x[i]*10,y[i]*10,x0,y0,xp,yp)
#        fig = plt.figure()
#        a = fig.add_subplot(1,1,1)
#        imshow(f)
#        a.set_ylim(0,220)

#        show()

    files = []
    for line in open(filelist, 'r'):
        files += [line.strip()]

        fwhm += [average]

        ellipse += [ellipse_frame]

    print('the mean orientation of this star is:',mean(theta))

    pickle.dump(fwhm, open(tag+'_fwhm_'+str(label)+'.p','wb'))
    pickle.dump(ellipse, open(tag+'_ellipse_'+str(label)+'.p','wb'))
    pickle.dump(theta, open(tag+'_theta_'+str(label)+'.p','wb'))

    theta += pickle.load(open(tag + 'theta' + label + '.p'))
    x += pickle.load(open(tag + 'fwhm_x' + label + '.p'))
    y += pickle.load(open(tag + 'fwhm_y' + label + '.p'))

    for i in range(0, len(x)):
        rad = theta[i] * np.pi / 180.0
        theta[i] = theta[i] * np.pi / 180.0
        theta[i] = theta[i] - 2 * np.pi * (int(theta[i] / (2 * np.pi)))
        # print theta[i]

    i = 0

    fwhm = []
    ellipse = []
    theta = []

    for label in labels:
        fwhm += [pickle.load(open(tag+'_fwhm_'+str(label)+'.p'))]
        ellipse += [pickle.load(open(tag+'_ellipse_'+str(label)+'.p'))]
        theta += [pickle.load(open(tag+'_theta_'+str(label)+'.p'))]


    focus_position = []
    tel_alt = []
    mjd = []

    for line in files:
        with pf.open(line) as imdata:
            if 1.0*i/binning == int(i/binning):
                focus_position +=[imdata[0].header['FCSR_PHY']]
                tel_alt +=[imdata[0].header['TEL_ALT']]
                mjd +=[imdata[0].header['MJD']]
        i += 1
        if i > binning*(len(ellipse[0])-1): break


    subplot(2, 1, 1)
    for fwhms in fwhm:
        plot(mjd,fwhms,'o')
    ylabel('fwhm (pix)')
    xlabel('MJD')

    subplot(2, 1, 2)
    plot(mjd,tel_alt,'o')
    ylabel('Altitude')
    xlabel('MJD')
    ylim(0,100)

    savefig(tag+'_timeseries.png', bbox_inches=0)
    close()

    for fwhms in fwhm:
        plot(tel_alt,fwhms,'o')
    ylabel('fwhm (pix)')
    xlabel('Altitude')
    xlim(0,100)
    savefig(tag+'_altitude.png', bbox_inches=0)
    close()

    for thetas in theta:
        plot(mjd,thetas,'o') 
    ylabel('theta')
    xlabel('MJD')
    savefig(tag+'_theta.png', bbox_inches=0)
    close()

    for ellipses in ellipse:
        plot(mjd,ellipses,'o') 
    ylabel('ellipticity')
    xlabel('MJD')
    savefig(tag+'_ellipse.png', bbox_inches=0)
    close()

#for i in cross_sec_x:
#  plot(i)
#  show()

#im_ani.save('im.mp4',metadata={'artist':'Tom'})
#  i = 0
#  for fwhms in fwhm:
#        i += 1
#        plot(focus_position,array(fwhms) + i,'o') 
#  ylabel('fwhm (pix)')
#  xlabel('focus position (mm)')
#  savefig('plot_bin/'+tag+'_focuspos.png', bbox_inches=0)
#  close()

def find_minor_axis(theta,sx,sy):

    # use math instead of numpy?

    a = ((np.cos(theta))**2.0)/(2*sx**2) + ((np.sin(theta))**2.0)/(2*sy**2)

    b = -np.sin(2.0*theta)/(4*sx**2.0) + sin(2.0*theta)/(4*sy**2.0)

    c = ((np.sin(theta))**2.0)/(2*sx**2) + ((np.cos(theta))**2.0)/(2*sy**2)

    w_1 = np.sqrt(np.log(2)) / sqrt(a + 2*b*tan(theta+np.pi/2.0) + c*(tan(theta+np.pi/2.0))**2.0)

    w_2 = np.sqrt(np.log(2)) / sqrt(a + 2*b*tan(theta) + c*(tan(theta))**2.0)

    if abs(theta) > 2.0*np.pi:
        theta = theta - 2.0*np.pi*int((theta / (2.0*np.pi)))

    if theta <0 :
        theta = theta + 2*np.pi


    w_1 = abs(w_1/np.sin(theta)) 
    w_2 = abs(w_2/np.cos(theta))

    if abs(w_1) > abs(w_2):
        w_major = w_1
        w_minor = w_2
    else:
        w_major = w_2
        w_minor = w_1

    ellipse = np.sqrt(1.0 - (w_minor / w_major)**2.0)

    average = (w_minor + w_major)/2.0

    xx = 0.1*np.arange(0,220)
    func = tan(theta+3*np.pi/2)
    yy = xx*func

    x0 = 11

    # print ellipse, average

    return ellipse, average

def call_find_fwhm(file,inputcat,factor,size,stars,tag=''):

    tag = file.rstrip('.fits')

    #labels = ['f_1','f_3','f_5','f_7','f_9']
    labels = ['f_1','f_2','f_3','f_4','f_5','f_6','f_7','f_8','f_9']

    fwhm_a = [[]]*len(labels)
    fwhm_b = [[]]*len(labels)
    fwhm = [[]]*len(labels)
    theta = [[]]*len(labels)

    # create array of zeros of length 55x55
    zero_array = np.zeros((factor*size,factor*size))

    fwhm = {'f_1':[],'f_2':[],'f_3':[],'f_4':[],'f_5':[],'f_6':[],'f_7':[],'f_8':[],'f_9':[]}
    fwhm_a = {'f_1':[],'f_2':[],'f_3':[],'f_4':[],'f_5':[],'f_6':[],'f_7':[],'f_8':[],'f_9':[]}
    fwhm_b = {'f_1':[],'f_2':[],'f_3':[],'f_4':[],'f_5':[],'f_6':[],'f_7':[],'f_8':[],'f_9':[]}
    theta = {'f_1':[],'f_2':[],'f_3':[],'f_4':[],'f_5':[],'f_6':[],'f_7':[],'f_8':[],'f_9':[]}
    data = {'f_1':zero_array.copy(),'f_2':zero_array.copy(),'f_3':zero_array.copy(),'f_4':zero_array.copy(),'f_5':zero_array.copy(),'f_6':zero_array.copy(),'f_7':zero_array.copy(),'f_8':zero_array.copy(),'f_9':zero_array.copy()}
    lengths = {'f_1':True,'f_2':True,'f_3':True,'f_4':True,'f_5':True,'f_6':True,'f_7':True,'f_8':True,'f_9':True}

    f1 = plt.figure()
    f2 = plt.figure()
    f3 = plt.figure()

    fwhm_extract(file,inputcat,factor,size,stars,labels,tag)

    for label in labels:

        # load in the data from the fwhm_extract function - array of sample squares
        dat = pickle.load(open(tag+label+'.p','rb'))

        data[label] += dat
        # normalise (haven't we already normalised before storing in fwhm_extract?)
        data[label] = data[label]/data[label].max()
        fwhm_a_frame, fwhm_b_frame, theta_frame, residuals, model = find_2dfwhm(data[label],factor,size)
        # print fwhm_a_frame, fwhm_b_frame, theta_frame, label
        
        fwhm_a[label] += [fwhm_a_frame]
        fwhm_b[label] += [fwhm_b_frame]
        fwhm[label] += [(fwhm_a_frame + fwhm_b_frame)/2.0]
        theta[label] += [theta_frame]

        # plot the results
        label_no = int(label.split('_')[-1])
        a1 = f1.add_subplot(3,3,label_no)
        ticks = factor*np.arange(size)
        a1.set_yticks(ticks)
        a1.set_yticklabels(np.arange(size))
        a1.set_xticks(ticks)
        a1.set_xticklabels(np.arange(size))
        reverse = (0,a1.get_ylim()[1] + factor)
        a1.set_ylim(reverse)
        cax = a1.imshow(data[label], interpolation='none',cmap='afmhot')
        a1.grid(True)
        center = factor*((size)/2.0)
        a1.plot(center,center,'gx')
        lengths[label] = False
        a1.get_xaxis().set_ticklabels([])
        a1.get_yaxis().set_ticklabels([])

        # ONLY COMPUTE THIS ONCE...
        if label_no == 5:
            stub = os.path.basename(file)
            av_fwhm_a = round(fwhm_a_frame,2)
            av_fwhm_b = round(fwhm_b_frame,2)

            try:
                eccentricity = round(np.sqrt(1.0-((av_fwhm_b/av_fwhm_a)**2)),2)
            except:
                print("Eccentricity is NaN")
                print(av_fwhm_a, av_fwhm_b)

            av_theta = round(theta_frame,0)



        a2 = f2.add_subplot(3,3,label_no)
        ticks = factor*np.arange(size)
        a2.set_yticks(ticks)
        a2.set_yticklabels(np.arange(size))
        a2.set_xticks(ticks)
        a2.set_xticklabels(np.arange(size))
        reverse = (0,a2.get_ylim()[1] + factor)
        a2.set_ylim(reverse)
        cax = a2.imshow(abs(residuals), interpolation='none',cmap='afmhot',vmin=min(data[label].flatten()),vmax=max(data[label].flatten()))
        a2.grid(True)
        center = factor*((size)/2.0)
        a2.plot(center,center,'gx')
        lengths[label] = False
        a2.get_xaxis().set_ticklabels([])
        a2.get_yaxis().set_ticklabels([])


        a3 = f3.add_subplot(3,3,label_no)
        ticks = factor*np.arange(size)
        a3.set_yticks(ticks)
        a3.set_yticklabels(np.arange(size))
        a3.set_xticks(ticks)
        a3.set_xticklabels(np.arange(size))
        reverse = (0,a3.get_ylim()[1] + factor)
        a3.set_ylim(reverse)
        cax = a3.imshow(model, interpolation='none',cmap='afmhot',vmin=min(data[label].flatten()),vmax=max(data[label].flatten()))
        a3.grid(True)
        center = factor*((size)/2.0)
        a3.plot(center,center,'gx')
        data[label] = zero_array.copy()
        lengths[label] = False
        a3.get_xaxis().set_ticklabels([])
        a3.get_yaxis().set_ticklabels([])

    f1.suptitle('PSF_a_5: {0:.2f}        PSF_b_5: {1:.2f}        e: {2:.2f}      theta: {3:.0f}'.format(av_fwhm_a, av_fwhm_b, eccentricity, av_theta),size=15)
    f2.suptitle('PSF_a_5: {0:.2f}        PSF_b_5: {1:.2f}        e: {2:.2f}      theta: {3:.0f}'.format(av_fwhm_a, av_fwhm_b, eccentricity, av_theta),size=15)
    f3.suptitle('PSF_a_5: {0:.2f}        PSF_b_5: {1:.2f}        e: {2:.2f}      theta: {3:.0f}'.format(av_fwhm_a, av_fwhm_b, eccentricity, av_theta),size=15)

    f1.text(0.5, 0.05, '{0}'.format(stub), ha="center", va="center",size=15)
    f2.text(0.5, 0.05, '{0}'.format(stub), ha="center", va="center",size=15)
    f3.text(0.5, 0.05, '{0}'.format(stub), ha="center", va="center",size=15)

    f1.savefig(file.rstrip('.fits')+'_psf.png', bbox_inches=0)
    f2.savefig(file.rstrip('.fits')+'_residuals.png', bbox_inches=0)
    f3.savefig(file.rstrip('.fits')+'_model.png', bbox_inches=0)

    plt.close()
    plt.close()
    plt.close()

    for label in labels:
        os.system('rm '+tag+label+'.p')

    cache = False
    if cache == True:
        for label in labels:
            pickle.dump(fwhm[label], open(tag+'fwhm'+label+'.p','wb'))
            pickle.dump(fwhm_a[label], open(tag+'fwhm_a'+label+'.p','wb'))
            pickle.dump(fwhm_b[label], open(tag+'fwhm_b'+label+'.p','wb'))
            pickle.dump(theta[label], open(tag+'theta'+label+'.p','wb'))

    return fwhm_a, fwhm_b, theta

def find_2dfwhm(data,factor,size):
    # this function
    #x1 = [1.0,np.pi,0.6,0.4,size/2.0,size/2.0]

    x1 = [1.0,np.pi,0.6,0.1,size/2.0,size/2.0]

    # np.arange - Return evenly spaced values within a given interval (in this case 0 to (factor*size) in steps of 1)
    # /factor to normalise the array from 0 to size with steps of (1/factor)
    x = np.arange(1.0*factor*size)/factor
    y = np.arange(1.0*factor*size)/factor

    # perform a least squares optimisation using scipy.optimize
    # Minimize the sum of squares of a set of equations. The parameters are as follows:
    # func : (callable) should take at least one (possibly length N vector) argument and returns M floating
    # point numbers. It must not return NaNs or fitting might fail.
    # x0 : (ndarray) The starting estimate for the minimization.
    # args : (tuple, optional) Any extra arguments to func are placed in this tuple.

    # gaussian2d_fit returns residuals from fitting a 2D gaussian to the data therefore this least squares procedure
    # will minimise the residuals
    p, success = opt.leastsq(gaussian2d_fit, x1, args=(x,y,data))

    # use the optimised values for the gaussian to create a model:
    model = gaussian2d(p[0],p[1],p[2],p[2]+abs(p[3]),p[4],p[5],x,y)

    # find the full width at half max (fwhm) values in
    # fwhm = 2*sqrt(2ln(2))*std_dev therefore fwhm_b is along short axis and fwhm_a is along long axis?
    fwhm_b = 2.0*(np.sqrt(2.0*np.log(2.0)))*p[2]
    fwhm_a = 2.0*(np.sqrt(2.0*np.log(2.0)))*(p[2] + abs(p[3]))
    # convert theta to degrees?
    theta = p[1]*180.0/np.pi

    #while theta > 360:
        #theta -= 180

    #while theta < 0:
        #theta += 180

    #plt.imshow(model, interpolation='none')
    #plt.show()

    # calculate the residuals as the model - actual data
    residuals = model -data

    return fwhm_a, fwhm_b, theta, residuals, model

def gaussian2d_fit(p,x,y,data):
    #x1 = [1.0, np.pi, 0.6, 0.1, size / 2.0, size / 2.0]
    #x = np.arange(1.0 * factor * size) / factor
    #y = np.arange(1.0 * factor * size) / factor
    # function which takes in (x1, x, y, data) and fits a 2D gaussian to this data

    f = gaussian2d(p[0],p[1],p[2],p[2]+abs(p[3]),p[4],p[5],x,y)

    # residuals is the difference between this gaussian function and the actual data
    # .flatten will convert the array given into one single array
    try:
        residuals = (data - f).flatten()/np.sqrt(abs(data.flatten()))
    except:
        print("Data is zero")
        print(abs(residuals))
    # replace any residuals which are infinite with 0
    residuals[abs(residuals) == np.inf] = 0.0

    return residuals

def gaussian2d(A,theta,sx,sy,x0,y0,x,y):
    # arguments when calling this function:
    #p = [1.0, np.pi, 0.6, 0.1, size / 2.0, size / 2.0]
    #x = np.arange(1.0 * factor * size) / factor
    #y = np.arange(1.0 * factor * size) / factor
    # therefore A = 1, theta = pi, sx = 0.6, sy = 0.7, x0 = size/2 and y0 = size/2
    # sx and sy are the spread of the blob and x0 and y0 are the central coords (which is centre of (x,y))

    # http://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function

    # create 2D grid of size x by y
    xx, yy = np.meshgrid(x,y)

    # these coeffs a, b and c are for the general form of f (see wiki)
    a = ((np.cos(theta))**2.0)/(2*sx**2) + ((np.sin(theta))**2.0)/(2*sy**2)

    b = -np.sin(2.0*theta)/(4*sx**2.0) + np.sin(2.0*theta)/(4*sy**2.0)

    c = ((np.sin(theta))**2.0)/(2*sx**2) + ((np.cos(theta))**2.0)/(2*sy**2)

    # 2D elliptical gaussian function:
    f = A*np.exp(-(a*(x0-xx)**2 + 2*b*(x0- xx)*(y0 - yy) + c*(y0 - yy)**2.0))
    return f

def fwhm_extract(image_name,inputcat,factor,size,stars,condition_name_list,tag=''):

    # fits open science image
    with pf.open(image_name) as imdata:
        image = imdata[0].data

        # increase size by 2 - why?
        size += 2

        # fits open the imcore list output
        with pf.open(image_name + '.phot') as photdata:

            # fits open stack image catalogue
            with pf.open(inputcat) as incat:
                mean_fluxes = incat[1].data['isophotal_flux']
                # subset of fluxes which are less than the median of fluxes greater than median of all fluxes (top half)
                # and greater than the median of fluxes less than the median of all fluxes (bottom half)
                # this means you end up with values between 1/4 and 3/4  of the sorted fluxes
                IQR = (mean_fluxes < (np.median(mean_fluxes[mean_fluxes > np.median(mean_fluxes)]))) & \
                      (mean_fluxes > (np.median(mean_fluxes[mean_fluxes < np.median(mean_fluxes)])))
                selection = np.argsort(mean_fluxes)[-(stars + 100):-100]

            # find the x and y coords of the imcore list catalogue
            xpos = photdata[1].data['X_coordinate']
            ypos = photdata[1].data['Y_coordinate']
            # extract the x and y positions of the selection
            # xpos = xpos[selection]
            # ypos = ypos[selection]
            xpos = xpos[IQR]
            ypos = ypos[IQR]

#        with pf.open('../focus_test/test.cat') as photdata:
#            mean_fluxes = photdata[2].data['FLUX_APER']
#            IQR = [(mean_fluxes < (np.median(mean_fluxes[mean_fluxes > median(mean_fluxes)]))) & (mean_fluxes > (median(mean_fluxes[mean_fluxes < median(mean_fluxes)])))]
#            xpos = photdata[2].data['X_IMAGE'][mean_fluxes > 3*mean(mean_fluxes)]
#            ypos = photdata[2].data['Y_IMAGE'][mean_fluxes > 3*mean(mean_fluxes)]

        imdata.close()

        # calculate phase by extracting the fractional part of the coordinates
        x_phase = np.array([abs(x-int(x)) for x in xpos])
        y_phase = np.array([abs(y-int(y)) for y in ypos])

        imdata.close()

#        n,bins,patches = hist(x_phase,10,normed=1)
#        show()
#        n,bins,patches = hist(y_phase,10,normed=1)
#        show()

        # mx = 2048
        # my = 2048

#        top_left = [(xpos < mean(xpos)) & (ypos > mean(ypos))],xpos[(xpos < mean(xpos)) & (ypos > mean(ypos))]
#        top_right = [(xpos > mean(xpos)) & (ypos > mean(ypos))],xpos[(xpos > mean(xpos)) & (ypos > mean(ypos))]
#        bottom_left = [(xpos < mean(xpos)) & (ypos < mean(ypos))],xpos[(xpos < mean(xpos)) & (ypos < mean(ypos))]
#        bottom_right = [(xpos > mean(xpos)) & (ypos < mean(ypos))],xpos[(xpos > mean(xpos)) & (ypos < mean(ypos))]
#        condition_list = [top_left,top_right,bottom_left,bottom_right]
#        condition_name_list = ['top_left','top_right','bottom_left','bottom_right']


        # mx = 2048
        # my = 2048
        mx = np.shape(image)[1]
        my = np.shape(image)[0]

        # (xpos < mx/3) is a boolean list for each xpos
        # therefore each f is a boolean array where all conditions are met
        #split pixel grid into thirds for both x and y - 1,2,3

        f_1 = (xpos < mx / 3) & (ypos > 2 * my / 3)
        f_2 = (xpos > mx / 3) & (xpos < 2 * mx / 3) & (ypos > 2 * my / 3)
        f_3 = (xpos > 2 * mx / 3) & (ypos > 2 * my / 3)
        f_4 = (xpos < mx / 3) & (ypos < 2 * my / 3) & (ypos > my / 3)
        f_5 = (xpos > mx / 3) & (xpos < 2 * mx / 3) & (ypos < 2 * my / 3) & (ypos > my / 3)
        f_6 = (xpos > 2 * mx / 3) & (ypos < 2 * my / 3) & (ypos > my / 3)
        f_7 = (xpos < mx / 3) & (ypos < my / 3)
        f_8 = (xpos > mx / 3) & (xpos < 2 * mx / 3) & (ypos < my / 3)
        f_9 = (xpos > 2 * mx / 3) & (ypos < my / 3)

        condition_list = [f_1,f_2,f_3,f_4,f_5,f_6,f_7,f_8,f_9]

        for x in condition_name_list:
            i = int(x.split('_')[-1])-1
            get_psf(ypos[condition_list[i]], xpos[condition_list[i]], image, size, factor, x, tag)

#        threads = []
#        for i in range(0,len(condition_list)):
#            t = threading.Thread(target=get_psf, args = (ypos[condition_list[i][0]],xpos[condition_list[i][0]],image,size,factor,condition_name_list[i],tag))
#            threads.append(t)
#            threads = np.append(threads,t)
#        [x.start() for x in threads]
#        [x.join() for x in threads]

def get_psf(ypos,xpos,image,size,factor,condition_name,tag=''):

    a_size = size*factor
    stack = np.zeros((a_size,a_size))
    mini = 2*size
    # maxi = 2048 -(2*size)
    # maxix = 2046-(2*size)
    # maxiy = 2044-(2*size)
    maxix = np.shape(image)[1] - (2*size)
    maxiy = np.shape(image)[0] - (2*size)
    # condition checks which positions lie between min and max i
    # we exclude regions where we lie within 2*size of any side
    condition = (ypos > mini) & (ypos < maxiy) & (xpos > mini) & (xpos < maxix)

    xpos = np.array(xpos[condition])
    ypos = np.array(ypos[condition])

    for n in range(0, len(xpos)):
        coords = [ypos[n], xpos[n]]
        # print 'n: ' + str(n)
        # print 'coords: ' + str(coords)

        # oversampled is an array which represents a sample square centred on coords
        oversampled = return_sample_square(coords, image, size, factor)
        # repeat for all coords and store in stack:
        stack += oversampled
    # cleanup

    #stack_copy = stack
    # reduce stack to elements between (5,5) and (59,59) - why?
    stack  = stack[factor:-factor,factor:-factor]

    # for (x,y), item in np.ndenumerate(stack_copy):
    #     if item == stack[-1,-1]:
    #         print (x,y)

    # normalise
    stack = stack - np.median(stack)
    stack = stack / stack.max()

    outname = tag+condition_name+'.p'

    # write variable stack to file
    pickle.dump(stack, open(outname,'wb'))


def return_sample_square(c,image,size,factor):
    # c - coordinates
    # r = 1/2 (n-1)
    r = (size-1.0)/2.0
    # round will round to nearest integer but return a float
    # look at square formed between (coordinates - (r+1)) and (coordinates + r) =  -(size/2 + 1/2) to +(size/2 - 1/2)
    raw = image[int(round(c[0])-r-1):int(round(c[0])+r),int(round(c[1])-r-1):int(round(c[1])+r)]

    # I changed these round because the coords are stored (y,x)
    # offsets are the fractional parts of the coordinates * factor
    yoffset = (round(c[0]) - c[0])*factor
    xoffset = (round(c[1]) - c[1])*factor

    oversampled = []
    for row in raw:
        oversampled_row = []
        # for every pixel in this square
        for pixel in row:
            # multiply by factor
            oversampled_row += [pixel]*factor
        # multiply each row by factor (why?)
        # oversampled will be a list of rows where each row is a list of pixels
        oversampled += [oversampled_row]*factor
    # convert to array
    oversampled = np.array(oversampled)
 
    # now recenter on the central pixel

    oversampled = recenter(oversampled,xoffset,yoffset)

    # now normalise the output

    maxval = oversampled.max() + 1.0

    oversampled = oversampled / maxval

    return oversampled

def recenter(oversampled,xshift,yshift):

    xshift = -int(round(xshift))
    yshift = -int(round(yshift))

    blanks = np.zeros((len(oversampled[:,]),abs(xshift))).T
    if xshift < 0:
        oversampled=np.append(blanks,oversampled[:xshift],axis=0)
    else:
        xshift = xshift
        blanks = np.zeros((len(oversampled[:,]),abs(xshift))).T
        oversampled=np.append(oversampled[xshift:],blanks,axis=0)

    blanks = np.zeros((len(oversampled[:,]),abs(yshift)))

    if yshift < 0:
        oversampled=np.append(blanks,oversampled[:,:yshift],axis=1)
    else:
        oversampled=np.append(oversampled[:,yshift:],blanks,axis=1)
    return oversampled

if __name__ == '__main__':

    description='''
    Produces diagnostic plot of the PSF by super-sampling and stacking the star images
    on the chip. Relies on good astrom to function correctly.

    Currently requires a list of files and a casutools imcore output file. (note, not a imcore_list
    file)
    Searches for a 'catcache' directory in the same path as the image files, this is used for the
    comparison.
    Should modify so a direct path to the catanp.log file is given.
    '''

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('filelist')
    parser.add_argument('inputcat')
    parser.add_argument('outname')
    parser.add_argument('--factor', type=int, required=False, default=5,
                                            help='What oversampling factor to use [default: 5]')
    parser.add_argument('--size', type=int, required=False, default=11,
                                            help='how large a region around each star to stack [default: 11]')
    parser.add_argument('--stars', type=int, required=False, default=100,
                                            help='How many stars to stack in each quadrant [default: 100]')
    parser.add_argument('--binning', type=int, required=False, default=1,
                                            help='use a binning factor for long time series? [default: 1]')
    parser.add_argument('--nproc', type=int, required=False, default=4,
                                            help='use multiprocessing? [default: 4]')

    args = parser.parse_args()

    super_sample(args.filelist,args.inputcat,args.factor,args.size,args.stars,args.binning,args.outname,args.nproc)

# vim: sw=2
