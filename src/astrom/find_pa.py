import math
import astropy.io.fits as fits
import argparse

class ccd:

    def load_astrometry(self, fname):
        hdulist = fits.open(fname)
        hdr = hdulist[0].header

        self.cd11 = hdr['CD1_1']
        self.cd12 = hdr['CD1_2']
        self.cd21 = hdr['CD2_1']
        self.cd22 = hdr['CD2_2']
        self.jd = hdr['JD']

        if 'ASTSOLVE' in hdr:
            self.solved = hdr['ASTSOLVE']
        else:
            self.solved = 'F'

    def determinant(self):
        self.det = (self.cd11*self.cd22) - (self.cd12*self.cd21)


def pa(fname):

    a = ccd()
    a.load_astrometry(fname)
    a.determinant()

    if a.det >= 0:
        parity = 1.
    else:
        parity = -1.

    T = parity * a.cd11 + a.cd22
    A = parity * a.cd21 - a.cd12

    # arctan(A/T), added the 180 degrees
    # orient = 180-math.degrees(math.atan2(A, T))
    orient = math.degrees(math.atan(A/T))
    # print math.degrees(math.atan2(A, T))
    # print math.degrees(math.atan(A/T))
    # print math.degrees(math.atan2(-A, T))
    # print math.degrees(math.atan(-A/T))
    # print math.degrees(math.atan2(-A, -T))
    # print math.degrees(math.atan(-A/-T))
    # print math.degrees(math.atan2(A, -T))
    # print math.degrees(math.atan(A/(-T)))


    return orient, a.jd, a.solved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fname')
    args = parser.parse_args()
    fname=args.fname
    pa(fname)

# def write_header(hdr,val,fname):
#
#     hdulist = fits.open(fname)
#     hdr_old = hdulist[0].header
#     data = hdulist_old[0].data
#     hdr_old.set(hdr, val)
#     fits.writeto(fname, data, hdr_old, overwrite=True)