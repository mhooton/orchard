import sys
import datetime as dt
import glob
import argparse

def read_file(fname, silent=False):
    f = open(fname, "r")
    content = f.read()

    splitlines = content.split("\n")

    try:
        if splitlines[3][0] == ";":
            spl_gaia = splitlines[3].split(" ")[1]
            if len(spl_gaia) > 16:
                if silent:
                    print(spl_gaia)
                else:
                    print("SUCCESSFULLY EXTRACTED GAIA ID FROM PLAN: " + spl_gaia)
                return spl_gaia
            else:
                print("N")
        else:
            print("N")
    except:
        print("N")
        pass

    return None


def find_plan(dir, date, targname):
    date_dt = dt.datetime.strptime(date, "%Y%m%d")
    date_formatted = date_dt.strftime("%Y-%m-%d")
    search_pattern = "%s/schedule/Plans_by_date/%s/Obj_%s.txt" % (dir, date_formatted, targname)
    flist = glob.glob(search_pattern)
    return flist


def gaia_id_from_schedule(dir, date, targname, silent=False):
    # Convert double-hyphens back to spaces for searching plan files
    targname = targname.replace('--', ' ')

    fplan = find_plan(dir, date, targname)

    if len(fplan) >= 1:
        gaia = read_file(fplan[0], silent=silent)
        if gaia is None:
            print("N")
    else:
        print("N")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dir', help='Observatory directory')
    parser.add_argument('date', help='Date in YYYYMMDD format')
    parser.add_argument('targname', help='Target name')
    parser.add_argument('--silent', action='store_true',
                        help='Only output the Gaia ID, no success messages')

    args = parser.parse_args()

    gaia_id_from_schedule(args.dir, args.date, args.targname, silent=args.silent)


if __name__ == '__main__':
    main()
