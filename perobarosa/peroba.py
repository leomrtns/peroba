import os, logging, argparse, sys, pathlib

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

current_working_dir = os.getcwd()

def run_align (args):
    from perobarosa import task_seq
    task_seq.align (args.fasta, args.alignment, args.reference)


class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    parser = ParserWithErrorHelp(
    description="""
    peroba top level
    """) 

    parser.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING, help="Print debugging statements")
    parser.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
    parser.add_argument('-o', '--outdir', action="store", help="Output database directory. Default: working directory")
    subp= parser.add_subparsers(dest='command')

    up_aln = subp.add_parser('align', help="add sequences to an alignment")
    #up_aln.add_argument('-s', '--fasta', metavar='fas', nargs='+', required=True, help="unaligned sequences")
    up_aln.add_argument('fasta', help="unaligned sequences")
    up_aln.add_argument('-a', '--alignment', metavar='aln', nargs="+", help="optional file with aligned sequences")
    up_aln.add_argument('-r', '--reference', metavar='fas', help="optional file with reference genome (default=MN908947.3)")
    up_aln.set_defaults(func = run_align)

    args = parser.parse_args()
    args.func(args)

    logging.basicConfig(level=args.loglevel)
    if args.outdir: 
        args.outdir = os.path.join(current_working_dir, args.outdir)
        common.pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist
    else: 
        args.outdir = current_working_dir

    #args.func(args) # calls task 

if __name__ == '__main__':
    main()
