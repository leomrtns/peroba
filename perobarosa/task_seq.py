#!/usr/bin/env python
from perobarosa.utils_seq import *

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='task_%(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

default_ref_seq = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/MN908947.3.fas") 

def align (fasta, alignment = None, reference = None):
    if reference is None:
        reference = default_ref_seq
    if alignment is None:
