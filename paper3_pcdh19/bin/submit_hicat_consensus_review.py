#!/usr/bin/env python3
"""Thin entry point for the configured HiCAT stage-03 saved-consensus review."""
from pathlib import Path
import sys
bundle=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(bundle/'scripts'))
from hicat.consensus_review_submission import main
if __name__=='__main__':main(bundle)
