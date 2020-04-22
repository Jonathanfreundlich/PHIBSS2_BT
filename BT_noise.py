# GET NOISE FROM IMAGE THROUGH A GAUSSIAN FIT

print 'Running BT_noise.py'

from astropy.modeling import models, fitting

def get_noise(data,fthreshold=5,savefigure=False):
    # RETURNS THE SKY OFFSET AND THE NOISE STANDARD DEVIATION
    # offset, noise = get_noise(data)
    threshold=fthreshold*median(data)
    data_values=data.reshape(1,size(data))
    data_values=data_values[where(data_values<threshold)]
    nbins=100
    #
    figure('noise')
    clf()
    counts, bins, patches=hist(data_values,bins=nbins,color='b')
    meanbins=zeros(size(counts))
    for i in range(size(counts)):
        meanbins[i]=(bins[i]+bins[i+1])/2.
    p_init = models.Gaussian1D(amplitude=counts.max(),mean=median(data),stddev=std(bins)/2.)
    fit_p = fitting.LevMarLSQFitter()
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        p = fit_p(p_init, meanbins, counts)
    #
    plot(meanbins,p(meanbins),color='r',linewidth=2)
    title(r'offset = %g, noise = %g'%(p.mean.value, p.stddev.value))
    print 'offset =', p.mean.value, 'noise =', p.stddev.value
    #
    if savefigure:
        savefig('noise.pdf')
    #
    return p.mean.value, p.stddev.value


    
    