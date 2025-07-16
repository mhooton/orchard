from astropy.coordinates import SkyCoord
from astropy import units as u
from calibration.pipeutils import open_fits_file
import sys
import glob


def fix_target_name(targetname, ra, dec):
    targetname = targetname.split('BIS')[0]
    targetname = targetname.replace("SPC", "SP")
    targetname = targetname.replace("SPEC", "SP")
    targetname = targetname.replace("_", "-")
    targetname = targetname.split(".")[0]
    if '+' in targetname:
        targsplit = targetname.split("+")

    else:
        targsplit = targetname.split("-")

    # depending on the formatting of RA/Dec (from headers) get the 'corrected' targetname
    # this name may have issues if the PM is high, account for that in next step below
    if " " not in str(ra):
        c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
        ra = c.ra.hms
        ra = str(int(ra[0])).zfill(2) + str(int(ra[1])).zfill(2)
        dec = c.dec.dms
        if int(dec[0]) > 0:
            sign = '+'
            hr = str(int(dec[0])).zfill(2)
            min = str(int(dec[1])).zfill(2)
        else:
            sign = "-"
            hr = str(int(dec[0]))[1:].zfill(2)
            min = str(int(dec[1]))[1:].zfill(2)
        dec = sign + hr + min
        expected_tname = "SP" + ra + dec
    else:
        ra = str(ra).replace(" ", "")
        dec = str(dec).replace(" ", "")
        expected_tname = "SP" + ra[:4] + dec[:5]
        sign = dec[0]

    # if the formatting and the sign was already correct, use the name in the header to avoid issues with PM & stacks
    if len(targsplit) > 1:
        if (len(targsplit[0]) == 6) and (len(targsplit[1]) == 4):
            if sign == targetname[6]:
                expected_tname = targetname
        else:
            # if the formatting or sign is wrong then correct the targetname but use the new name
            if len(targsplit[0]) >= 6 and len(targsplit[1]) >= 4:
                # if "-" in targetname:
                if sign == "-":
                    targetname = targsplit[0][:6] + "-" + targsplit[1][:4]
                elif sign == "+":  # +" in targetname:
                    targetname = targsplit[0][:6] + "+" + targsplit[1][:4]

                # expected_tname = targetname

    return expected_tname, targetname


def find_existing_stack(targetname, stackdir):
    exists = False
    tname = targetname
    flag = False

    if "-" in targetname:
        t_split = targetname.split("-")
        splitter = "-"
    elif "+" in targetname:
        t_split = targetname.split("+")
        splitter = "+"

    t_split[0] = t_split[0][2:]
    t_split = [int(t) for t in t_split]

    p1 = "SP" + str(t_split[0] - 2) + splitter + str(t_split[1])
    p2 = "SP" + str(t_split[0] - 1) + splitter + str(t_split[1])
    p3 = "SP" + str(t_split[0] + 1) + splitter + str(t_split[1])
    p4 = "SP" + str(t_split[0] + 2) + splitter + str(t_split[1])
    p5 = "SP" + str(t_split[0]) + splitter + str(t_split[1] - 2)
    p6 = "SP" + str(t_split[0]) + splitter + str(t_split[1] - 1)
    p7 = "SP" + str(t_split[0]) + splitter + str(t_split[1] + 1)
    p8 = "SP" + str(t_split[0]) + splitter + str(t_split[1] + 2)
    p9 = "SP" + str(t_split[0] - 2) + splitter + str(t_split[1] - 2)
    p10 = "SP" + str(t_split[0] - 2) + splitter + str(t_split[1] - 1)
    p11 = "SP" + str(t_split[0] - 2) + splitter + str(t_split[1] + 1)
    p12 = "SP" + str(t_split[0] - 2) + splitter + str(t_split[1] + 2)
    p13 = "SP" + str(t_split[0] - 1) + splitter + str(t_split[1] - 2)
    p14 = "SP" + str(t_split[0] - 1) + splitter + str(t_split[1] - 1)
    p15 = "SP" + str(t_split[0] - 1) + splitter + str(t_split[1] + 1)
    p16 = "SP" + str(t_split[0] - 1) + splitter + str(t_split[1] + 2)
    p17 = "SP" + str(t_split[0] + 1) + splitter + str(t_split[1] - 2)
    p18 = "SP" + str(t_split[0] + 1) + splitter + str(t_split[1] - 1)
    p19 = "SP" + str(t_split[0] + 1) + splitter + str(t_split[1] + 1)
    p20 = "SP" + str(t_split[0] + 1) + splitter + str(t_split[1] + 2)
    p21 = "SP" + str(t_split[0] + 2) + splitter + str(t_split[1] - 2)
    p22 = "SP" + str(t_split[0] + 2) + splitter + str(t_split[1] - 1)
    p23 = "SP" + str(t_split[0] + 2) + splitter + str(t_split[1] + 1)
    p24 = "SP" + str(t_split[0] + 2) + splitter + str(t_split[1] + 2)

    possible_dups = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21,
                     p22, p23, p24]
    # print(possible_dups)

    for p in possible_dups:
        stacks = glob.glob("%s/%s_outstack_*.fits" % (stackdir, p))
        if len(stacks) > 0:
            if exists == True:
                flag = True
            exists = True
            tname = p
            # print("This target already exists under a different name! " + p)

    return exists, tname


def correct_target_names(imlist, targ, outdir):
    # Convert double-hyphens back to spaces for processing
    targ = targ.replace('--', ' ')

    stackdir = outdir + "/StackImages"

    with open(imlist) as infile:
        # loop through all science fits filenames
        fnames = [line.strip() for line in infile]

    if len(fnames) > 1:
        try:
            with open_fits_file(fnames[1]) as hdulist:
                # tname = hdulist[0].header['object']
                ra = hdulist[0].header['ra']
                dec = hdulist[0].header['dec']

            targ = targ.upper()

            if "SP" in targ and "WASP" not in targ and "TESS" not in targ and "NGTS" not in targ:
                # print targ, ra, dec
                exp_targ, targetname = fix_target_name(targ, ra, dec)
                dup, tname_dup = find_existing_stack(exp_targ, stackdir)
                # targ = exp_targ
                if dup:
                    targ = tname_dup
                else:
                    targ = exp_targ  # targetname
            elif "TRAPPIST" in targ:
                # targ = "TRAPPIST-1"
                targ = "SP2306-0502"
            print(targ.replace(' ', '--'))
        except:
            print(targ.upper().replace(' ', '--'))
    else:
        print(targ.upper().replace(' ', '--'))


def main():
    imlist = sys.argv[1]
    targ = sys.argv[2]
    outdir = sys.argv[3]

    correct_target_names(imlist, targ, outdir)


if __name__ == '__main__':
    main()
