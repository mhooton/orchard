import fitsio
import os
import glob
from astropy.io import fits, ascii
import numpy as np
import sys


def import_gaia_ids_40pc(fname):
    # 40 PC TARGET LIST
    data = ascii.read(fname, delimiter=" ", header_start=0, data_start=1)
    gaia_ids = [str(x) for x in data['Gaia_ID,'][1:]]
    sp_ids = [str(d).upper() for d in data['Sp_ID,'][1:]]
    return gaia_ids, sp_ids


def create_targetlist(runtlist, stackdir, outdir, fname_tlist):
    accept_filts = ['Clear', 'i', 'I+z', 'z', 'g', 'r', 'B', 'Exo', 'V', 'R']
    run_gaia, run_sp = [], []

    # get list of all output.fits files
    outfits = glob.glob("%s/*_*_output.fits" % (outdir))
    targs = [os.path.basename(o).split("_")[0] for o in outfits if os.path.basename(o).split("_")[1] in accept_filts]
    filts = [os.path.basename(o).split("_")[1] for o in outfits if os.path.basename(o).split("_")[1] in accept_filts]

    # remove duplicate targs
    utargs = list(set(targs))

    # get all Gaia IDs and their SPECULOOS IDs
    gaia_ids_tlist, spec_ids = import_gaia_ids_40pc(fname_tlist)
    spec_ids = [s.upper() for s in spec_ids]

    stacks = glob.glob("%s/*_outstack_*.fits" % (stackdir))
    stack_targs = [os.path.basename(o).split("_")[0] for o in stacks if os.path.basename(o).split("_")[0][:2] == "SP"]
    ustack_targs = list(set(stack_targs))

    in_stack = []
    # for every run target add to list with names
    for t in range(len(utargs)):
        if utargs[t] in gaia_ids_tlist:
            # if the target is in the 40pc targetlist get SPECULOOS name
            ind = np.where(np.array(gaia_ids_tlist) == utargs[t])[0][0]
            # print(targs[t],spec_ids[ind])
            in_stack.append((utargs[t], spec_ids[ind]))
            run_gaia.append(utargs[t].upper())
            run_sp.append(spec_ids[ind])

        else:
            # the target isn't in the 40pc targetlist
            if utargs[t].upper() in ustack_targs:
                # in_stack.append(targs[t].upper())
                # print(targs[t],'IN STACK')
                # run_gaia.append("N")
                # run_sp.append(utargs[t].upper())
                if utargs[t].upper() in spec_ids:
                    ind = np.where(np.array(spec_ids) == utargs[t].upper())[0][0]
                    # print(utargs[t].upper(), gaia_ids_tlist[ind])
                    if gaia_ids_tlist[ind] != "0":
                        run_gaia.append(gaia_ids_tlist[ind])
                        run_sp.append(utargs[t].upper())
                    else:
                        run_gaia.append("N")
                        run_sp.append(utargs[t].upper())
                else:
                    run_gaia.append("N")
                    run_sp.append(utargs[t].upper())
            else:
                # print(utargs[t],"NOT IN STACK")
                # if it has a gaia ID
                if len(utargs[t]) > 17:
                    gaia_files = glob.glob("%s/20*/*/%s_*_output.fits" % (outdir, utargs[t]))
                    # print(gaia_files)
                    if len(gaia_files) == 1:
                        run_gaia.append(utargs[t].upper())
                        run_sp.append(gaia_files[0].split("/")[-2].upper())
                        # print(utargs[t], gaia_files[0].split("/")[-2],"\n")
                    else:
                        tmp_names = []
                        for g in gaia_files:
                            # print(g.split("/")[-2].upper())
                            tmp_names.append(g.split("/")[-2].upper())

                        # print(tmp_names)
                        if len(list(set(tmp_names))) == 1:
                            # only one unique name
                            run_gaia.append(utargs[t].upper())
                            run_sp.append(tmp_names[0])
                        else:
                            # print("MULTIPLE NAMES FOR SAME TARGET")
                            run_gaia.append(utargs[t].upper())
                            run_sp.append(tmp_names[0])

                else:
                    if utargs[t].upper() in spec_ids:
                        ind = np.where(np.array(spec_ids) == utargs[t].upper())[0][0]
                        if gaia_ids_tlist[ind] != "0":
                            # print(utargs[t].upper(),gaia_ids_tlist[ind])
                            run_gaia.append(gaia_ids_tlist[ind])
                            run_sp.append(utargs[t].upper())
                        else:
                            run_gaia.append("N")
                            run_sp.append(utargs[t].upper())
                    else:
                        run_gaia.append("N")
                        run_sp.append(utargs[t].upper())

    # what's the case where there's a stack but no outfits?
    # for t in range(len(ustack_targs)):
    #     if ustack_targs[t] not in run_sp:
    #         print("CHECK ",ustack_targs[t])
    #         for i in range(len(run_sp)):
    #             if run_sp[i] in ustack_targs[t]:
    #                 print("Can remove ",ustack_targs[t],run_sp[i],i,t)

    run_gaia = [j for i, j in sorted(zip(run_sp, run_gaia))]
    run_sp = sorted(run_sp)

    if not os.path.exists(runtlist):
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(runtlist), exist_ok=True)
        # Create the file with the header content
        with open(runtlist, 'w') as f:
            f.write("Gaia_ID SPECULOOS_ID\n\n")

    ascii.write([run_gaia, run_sp], runtlist, names=['Gaia_ID', 'SPECULOOS_ID'], overwrite=True)
    return


def get_runtargetlist(runtlist):
    data = ascii.read(runtlist, header_start=0, data_start=1)
    return data


def add_target_to_list(runtargs, targ, gaia, runtlist):
    run_sp = [str(s) for s in runtargs['SPECULOOS_ID'][1:]]
    run_gaia = [str(s) for s in runtargs['Gaia_ID'][1:]]

    run_sp.append(targ)
    run_gaia.append(gaia)

    run_gaia = [j for i, j in sorted(zip(run_sp, run_gaia))]
    run_sp = sorted(run_sp)

    ascii.write([run_gaia, run_sp], runtlist, names=['Gaia_ID', 'SPECULOOS_ID'], overwrite=True)

    return


def find_target_in_list(runtargs, targ, gaia):
    if targ in runtargs['SPECULOOS_ID'] or gaia in runtargs['Gaia_ID']:
        return True, targ
    else:
        for i in runtargs['SPECULOOS_ID']:
            if i in targ:
                return True, i

        return False, targ


def already_run_targs(fname_tlist, outdir, targ, gaia, add_targ=True):
    # Convert double-hyphens back to spaces for processing
    targ = targ.replace('--', ' ')

    # print(basedir,tel,v,targ,gaia)

    # outdir = basedir + "/output/" + v + "/"
    stackdir = outdir + "/StackImages/"
    runtarglist = stackdir + "/already_run.txt"
    # fname_tlist = os.path.dirname(basedir) + "/ml_40pc.txt"

    # if os.path.exists(targlist):
    #     os.remove(targlist)
    alreadyrun = False
    newtarg = targ

    if os.path.exists(runtarglist):

        # import already run targets
        runtargs = get_runtargetlist(runtarglist)

        # check if the current target is in that list
        inlist, newtarg = find_target_in_list(runtargs, targ.upper(), gaia)
        if add_targ == True and inlist == False:
            add_target_to_list(runtargs, targ.upper(), gaia, runtarglist)

        if inlist:
            alreadyrun = True

    else:
        # if the targetlist doesn't exist then create it
        create_targetlist(runtarglist, stackdir, outdir, fname_tlist)

    print(newtarg.replace(' ', '--'))


def main():
    fname_tlist = sys.argv[1]
    outdir = sys.argv[2]
    targ = sys.argv[3]
    gaia = sys.argv[4]

    already_run_targs(fname_tlist, outdir, targ, gaia)

if __name__ == '__main__':
    main()
