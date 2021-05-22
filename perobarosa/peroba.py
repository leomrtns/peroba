import os, logging, argparse, sys, pathlib, multiprocessing, datetime

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba_TOP  %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M%S")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

defaults = {
    "current_dir": os.getcwd() + "/",
    "timestamp": datetime.datetime.now().strftime("%y%m%d_%H%M%S"),
    "n_threads": multiprocessing.cpu_count(),
    "reference": os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/MN908947.3.fas") 
    }

def run_align (args):
    from perobarosa import task_seq
    if args.reference: defaults["reference"] = args.reference

    if args.length is None: args.length = 28000 
    if args.length < 10000: 
        logger.warning (f"Length {args.length} is way too short for genome alignment; changing to default 28k");
        args.length = 28000
    task_seq.align (args.fasta, defaults, args.alignments, args.output, args.length, args.ambiguous)

class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    parser = ParserWithErrorHelp()
    #description="""
    #peroba top level
    #""") 

    parser.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING, help="Print debugging statements")
    parser.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
    parser.add_argument('-t', '--nthreads', metavar='int', help="Number of threads requested (default = maximum available)")
    parser.add_argument('--outdir', action="store", help="Output database directory. Default: working directory")
    subp= parser.add_subparsers(dest='command')

    up_aln = subp.add_parser('align', help="add sequences to an alignment")
    #up_aln.add_argument('-s', '--fasta', metavar='fas', nargs='+', required=True, help="unaligned sequences")
    up_aln.add_argument('fasta', help="unaligned sequences")
    up_aln.add_argument('-a', '--alignments', metavar='aln', nargs="+", help="optional files with aligned sequences")
    up_aln.add_argument('-r', '--reference', metavar='fas', help="optional file with reference genome (default=MN908947.3)")
    up_aln.add_argument('-o', '--output', metavar='aln', help="incremental output alignment (i.e. only new sequences)")
    up_aln.add_argument('-A', '--ambiguous', metavar='float', help="maximum allowed ambiguity for uvaia (default = 0.1)")
    up_aln.add_argument('-l', '--length', metavar='int', help="exclude sequences shorter than this (default = 25k)")
    up_aln.set_defaults(func = run_align)

    args = parser.parse_args()

    logging.basicConfig(level=args.loglevel)
    if args.outdir: 
        defaults["current_dir"] = args.outdir = os.path.join(defaults["current_dir"], args.outdir)
        common.pathlib.Path(defaults["current_dir"]).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist
    if args.nthreads:
        if args.nthreads < 1: args.nthreads = 1
        defaults["n_threads"] = args.nthreads

    args.func(args) # calls task 

if __name__ == '__main__':
    main()
