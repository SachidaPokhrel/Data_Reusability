#!/usr/bin/env python3
"""Project CLI: stage execution, source acquisition, validation and audit."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pe_upstream.cli import main
if __name__ == '__main__':
    main()
