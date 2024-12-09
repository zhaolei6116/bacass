import argparse
import csv


__version__ = "v0.0.1"
__date__    = "2024.11.25"
__author__  = "zhaolei"


def gc_parse(infile):

    gc_ratio = {}
    contig_list = []
    with open(infile, 'r', encoding='utf-8') as inF:
        reader = csv.reader(inF, delimiter='\t')
        for row in reader:
            gc_ratio[row[0]] = (int(row[4]) + int(row[5]))/int(row[1]) * 100
            gc_ratio[row[0]] = round(gc_ratio[row[0]], 2)
            contig_list.append(row[0])
    
    return gc_ratio, contig_list


def cov_dep_parse(infile):

    cov_dep = {}
    with open(infile, 'r', encoding='utf-8') as inF:
        reader = csv.reader(inF, delimiter='\t')
        next(reader)
        for row in reader:
            cov_dep[row[0]] = {}
            cov_dep[row[0]]['length'] = row[2]
            cov_dep[row[0]]['cov'] = row[5]
            cov_dep[row[0]]['dep'] = row[6]
    
    return cov_dep


def get_result(gc_ratio, contig_list, cov_dep, outfile):

    head_line = ["Contig", "Length", "GC content(%)", "Coverage", "Sequence depth"]

    with open(outfile, "w") as outF:
        outF.write("\t".join(head_line) + "\n")
        for contig in contig_list:
            temp_list = [contig, cov_dep[contig]['length'], str(gc_ratio[contig]), cov_dep[contig]['cov'], cov_dep[contig]['dep']]
            outF.write("\t".join(temp_list) + "\n")
    

def setup_args():
    parser = argparse.ArgumentParser(description = 'get contig stat',
            formatter_class=argparse.RawTextHelpFormatter,
            epilog = f"Version: {__version__}\nDate: {__date__}\nAuthor: {__author__}")
    parser.add_argument('-v','--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('-gc', type=str, help='The path to gc stat', required = True)
    parser.add_argument('-cov', type=str, help='The path to cov and depth stat', required = True)
    parser.add_argument('-o', type=str, help='output file', default="contig_stat.xls")
    args = parser.parse_args()

    return args


def main():
    args = setup_args()
    gc_ratio, contig_list = gc_parse(args.gc)
    cov_dep = cov_dep_parse(args.cov)
    get_result(gc_ratio, contig_list, cov_dep, args.o)



if __name__ == "__main__":
    main()

