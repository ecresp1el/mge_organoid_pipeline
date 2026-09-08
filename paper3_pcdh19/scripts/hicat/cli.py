"""Execute the frozen HiCAT technical pilot; never choose a latest input."""
from pathlib import Path
import argparse
from .workflow import PilotWorkflow


def main():
    """Parse the explicit run directory, run the pilot, and print its asset path."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    args=parser.parse_args()
    print('Published pilot:',PilotWorkflow(args.run_dir).run(),flush=True)


if __name__=='__main__':
    main()
